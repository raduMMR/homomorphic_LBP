#include <iostream>
#include "basic_primitives.h"
using namespace std;

double clock_diff(const clock_t &t1, const clock_t &t2){
    return double(t2 - t1) / CLOCKS_PER_SEC;
}

vector<Ctxt*> hom_LBP(vector<Ctxt*> enc_pixeli, vector<vector<Ctxt*>> vecini, int t_bits);


void test_LBP() {
    assert(NSLOTS != 0);

    int t_bits=8;

    vector<long> pixeli;
    for(int i=0; i<NSLOTS; i++) {
        pixeli.push_back(rand() % 256);
    }

    vector<vector<long>> vecini(t_bits);
    for(int i=0; i<t_bits; i++) {
        for(int j=0; j<NSLOTS; j++) {
            // toti vecinii de pe pozitia i ai pixelilor.
            vecini[i].push_back(rand() % 256);
        }
    }

    vector<long> lbp_codes(NSLOTS);
    for(int i=0; i<NSLOTS; i++) {
        lbp_codes[i] = 0;
        for(int j=0; j<t_bits; j++) {
            lbp_codes[i] |= (vecini[j][i] >= pixeli[i]) << j;
        }
    }

    // encrypting pixels.
    cout << "Encrypting pixels ...\n";
    vector<vector<Ctxt*> > enc_vecini(t_bits);
    for(int i=0; i<t_bits; i++) {
        enc_vecini[i] = encryptIntVal(vecini[i], t_bits);
    }
    vector<Ctxt*> enc_pixeli = encryptIntVal(pixeli, t_bits);
    cout << "Done pixels encryption.\n";

    // computing LBP codes.
    cout << "Homomorphic LBP computation ...\n";
    vector<Ctxt*> hom_lbp = hom_LBP(enc_pixeli, enc_vecini, t_bits);
    cout << "Done LBP computation.\n";

    // comparison.
    bool success = true;
    vector<long> dec_lbp = decryptIntVal(hom_lbp);
    for(int i=0; i<dec_lbp.size(); i++) {
        if(dec_lbp[i] != lbp_codes[i]) {
            cout << "ESEC\n";
            success = false;
            break;
        }
    }
    if( success == true) {
        cout << "Succes!!!!!!!!\n";
    }

    // cleaning up.
    cout << "Cleaning up ctxts...\n";
    for(int i=0; i<t_bits; i++) {
        for(int j=0; j<t_bits; j++) {
            delete enc_vecini[i][j];
        }
    }

    for(int i=0; i<t_bits; i++) {
        delete enc_pixeli[i];
    }
    cout << "Done cleaning up ctxts.\n";

}


void test_Compute_s() {

    int t_bits = 8;
    int val1, val2;
    vector<Ctxt*> nr1;
    vector<Ctxt*> nr2;

    for(int i=0; i<128; i++) {
        val1 = rand() % 256;
        val2 = rand() % 256;

        vector<long> batch1(NSLOTS, val1);
        vector<long> batch2(NSLOTS, val2);

        nr1 = encryptIntVal(batch1, t_bits);
        nr2 = encryptIntVal(batch2, t_bits);

        Ctxt *gte = compute_s(0, t_bits, nr1, nr2);
        // cout << val1 << " >= " << val2 << " => " << decryptBitVal(gte) << endl;

        bool success = true;
        vector<long> plain_gte = decryptBitVal(gte);
        for(int j=0; j<NSLOTS; j++) {
            if(plain_gte[j] != (val1 >= val2)) {
                cout << "Rezultat gresit, slot " << j << endl;
                success = false;
                break;
            }
        }

        // cleaning.
        delete gte;
        for(int j=0; j<nr1.size(); j++) {
            delete nr1[j];
        }
        for(int j=0; j<nr2.size(); j++) {
            delete nr2[j];
        }

        if(success == true) {
            cout << "Test " << i << " incheiat cu succes.\n";
        }
        else {
            cout << "Test " << i << " esuat\n";
            break;
        }
    }

}

int main(int argc, char **argv) {

    ArgMapping amap;

    bool dry=false;
    amap.arg("dry", dry, "dry=1 for a dry-run");

    long R=1;
    amap.arg("R", R, "number of rounds");

    long p=2;
    amap.arg("p", p, "plaintext base");

    long r=1;
    amap.arg("r", r,  "lifting");

    long d=1;
    amap.arg("d", d, "degree of the field extension");
    amap.note("d == 0 => factors[0] defines extension");

    long c=2;
    amap.arg("c", c, "number of columns in the key-switching matrices");

    
    long k=80;
    amap.arg("k", k, "security parameter");

    long L=0;
    amap.arg("L", L, "# of levels in the modulus chain",  "heuristic");

    long s=0;
    amap.arg("s", s, "minimum number of slots");

    long repeat=1;
    amap.arg("repeat", repeat,  "number of times to repeat the test");

    long chosen_m=0;
    amap.arg("m", chosen_m, "use specified value as modulus", NULL);

    Vec<long> mvec;
    amap.arg("mvec", mvec, "use product of the integers as  modulus", NULL);
    amap.note("e.g., mvec='[5 3 187]' (this overwrite the m argument)");

    Vec<long> gens;
    amap.arg("gens", gens, "use specified vector of generators", NULL);
    amap.note("e.g., gens='[562 1871 751]'");

    Vec<long> ords;
    amap.arg("ords", ords, "use specified vector of orders", NULL);
    amap.note("e.g., ords='[4 2 -4]', negative means 'bad'");

    long seed=0;
    amap.arg("seed", seed, "PRG seed");

    long nt=1;
    amap.arg("nt", nt, "num threads");

    amap.arg("noPrint", noPrint, "suppress printouts");

    amap.parse(argc, argv);

    SetSeed(ZZ(seed));
    SetNumThreads(nt);
    
    if (L==0) { // determine L based on R,r
        L = 3*R+3;
        if (p>2 || r>1) { // add some more primes for each round
        long addPerRound = 2*ceil(log((double)p)*r*3)/(log(2.0)*FHE_p2Size) +1;
        L += R * addPerRound;
        }
    }

    long w = 64; // Hamming weight of secret key
    //  long L = z*R; // number of levels

    if (mvec.length()>0)
        chosen_m = computeProd(mvec);
    std::cout << argv[0] << ": ";
    long m = FindM(k, L, c, p, d, s, chosen_m, !noPrint);

    setDryRun(dry);
    cout << "Generare context ...\n";
    setGlobalVariables(p, r, d, c, k, w, L, m, gens, ords);
    cout << "Terminat de generat context.\n";

    // test_Compute_s();
    clock_t begin = clock();
    test_LBP();
    clock_t end = clock();

    cout << "TIMP: " << clock_diff(begin, end) << " secunde.\n";

    cout << "Cleaning up ...\n";
    cleanGlobalVariables();
    cout << "Terminat cleaning up.\n";

    cout << "Program terminat.\n";
    return 0;
}

/*************************************************************************************/
vector<Ctxt*> hom_LBP(vector<Ctxt*> enc_pixeli, vector<vector<Ctxt*>> vecini, int t_bits) {
    vector<Ctxt*> lbp_codes;
    for(int i=0; i<t_bits; i++) {
        lbp_codes.push_back(compute_s(0, t_bits, vecini[i], enc_pixeli));
    }
    return lbp_codes;
}

/*************************************************************************************/
