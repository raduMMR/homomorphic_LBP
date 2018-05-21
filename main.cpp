#include <iostream>
#include <tuple>
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include "basic_primitives.h"
using namespace std;

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

void refreshCtxt(Ctxt *&ctxt) {
    vector<long> bits;
    bits =  decryptBitVal(ctxt);
    delete ctxt;
    ctxt = encryptBitVal(bits);
}

class HE_INT {

public:

    HE_INT() {
        initialize_HE_INT();
    }

    HE_INT(const HE_INT& he_int) {
       clone(he_int);
    }

    HE_INT& operator=(const HE_INT& he_int) {
        clone(he_int);
        return (*this);
    }

    void clone(const HE_INT& he_int) {
        initialize_HE_INT();
        for(int i=0; i<enc_bits.size(); i++) {
            (*this->enc_bits[i]) = *(he_int.enc_bits[i]);
        }
    }

    vector<Ctxt*> enc_bits;

    // make HE_INT an encrypted representation of zero,
    // with enc_bits[0] = LSB and enc_bits[T_BITS-1] = MSB.
    void initialize_HE_INT() {
        vector<long> zeroes(NSLOTS, 0);
        this->enc_bits = encryptIntVal(zeroes, T_BITS);
    }
    

    // TO DO: adauga metoda copy unui Ctxt.
    void add1Bit(Ctxt* ctxt) {
        Ctxt *carry = encryptBitVal(vector<long>(NSLOTS, 0));
        *carry = (*enc_bits[0]);
        carry->multiplyBy((*ctxt));

        enc_bits[0]->addCtxt(*ctxt);

        Ctxt *new_carry = encryptBitVal(vector<long>(NSLOTS, 0));
        for(int i=1; i<enc_bits.size(); i++) {
            *new_carry = (*enc_bits[i]);
            new_carry->multiplyBy(*carry);
            enc_bits[i]->addCtxt(*carry);
            // TO DO: use this costly refresh only when mandatory.
            refreshCtxt(enc_bits[i]);
            *carry = *new_carry;
        }

        delete carry;
        delete new_carry;
    }

    void setSkipping(const Ctxt &isSkipped) {
        for(int i=0; i<enc_bits.size(); i++) {
            enc_bits[i]->multiplyBy(isSkipped);
        }
    }

    ~HE_INT() {
        for(int i=0; i<enc_bits.size(); i++) {
            delete enc_bits[i];
        }
    }


};

typedef vector<Ctxt*> ENC_INT;

/**
 * @param enc_nums - vector of encrypted numbers on bits representation
 * @return - the number of occurences for each number in the vector
 */
vector<tuple<vector<Ctxt*>, HE_INT>> hom_counter(vector<ENC_INT> &enc_nums) {


    // this vector will be used to count the occurence for each ENC_INT.
    // an entry of frequencies vector will hold the enc number ENC_INT and 
    // its number of hits HE_INT.
    vector<tuple<ENC_INT, HE_INT>> frequencies;

    vector<long> zeroes(NSLOTS, 0);

    // for each number keep a parallel vector 
    // with a Ctxt corresponding to each number
    // the slot is initialized with encryption of zero
    // when the number is hit, mark the slot as "seen" => Enc(1) 
    ENC_INT isCounted = encryptIntVal(zeroes, enc_nums.size());

    Ctxt *isSkipped = encryptBitVal(zeroes);

    for(int i=0; i<enc_nums.size(); i++) {
        HE_INT frequency;

        (*isSkipped) = (*isCounted[i]);

        isSkipped->negate(); // for multiplication purpose, see below.
        // if the number has been previously hit, then multiply with 0
        // the frequency of the counting.

        for(int j=i; j<enc_nums.size()-1; j++) {
            Ctxt* areEqual = compute_z(0, T_BITS, enc_nums[i], enc_nums[j]);

            refreshCtxt(areEqual);

            frequency.add1Bit(areEqual);

            // mark number as counted.
            (*isCounted[j]) = (*areEqual); 

            delete areEqual;
        }

        // skip number if already counted.
        // that means to set the frequency to zero for another hit
        // of a previous number. When doing chi-square computation
        // simply do not take into account elements with zero frequency.
        frequency.setSkipping(*isSkipped);

        frequencies.push_back(make_tuple(enc_nums[i], frequency));
    }

    // cleaning up.
    delete isSkipped;

    for(int i=0; i<isCounted.size(); i++) {
        delete isCounted[i];
    }

    return frequencies;
}

int myrandom (int i) { return std::rand()%i;}

void test_hom_counter() {
    // generate a random vector cu NSLOTS values.
    vector<long> myvector(NSLOTS, 0);
    for(int i=0; i<NSLOTS; i++) {
        myvector[i] = rand() % 3;
    }
    // myvector is going to be encrypted in one vector<Ctxt*>
    // so we'll shuffle this vector to have different values
    // in the corresponding slots.

    int VEC_SIZE = 10;

    vector<vector<Ctxt*>> enc_nums(VEC_SIZE);
    for(int i=0; i<VEC_SIZE; i++) {
        cout << myvector[0] << " ";
        enc_nums[i] = encryptIntVal(myvector, T_BITS);
        random_shuffle ( myvector.begin(), myvector.end(), myrandom);
    }
    cout << endl;

    vector<tuple<vector<Ctxt*>, HE_INT>> freqs = hom_counter(enc_nums);

    cout << "freqs.size() = " << freqs.size() << endl;

    for(int i=0; i<freqs.size(); i++) {
        vector<long> number = decryptIntVal(std::get<0>(freqs[i]));
        vector<long> frecventa = decryptIntVal(std::get<1>(freqs[i]).enc_bits);
        cout << "Numar = " << number[0] << ", Frecventa =  " << frecventa[0] << endl;
        // the memory of the number is deleted below, here it is just a pointer,
        // for(int j=0; j<std::get<0>(freqs[i]).size(); j++) {
        //     delete std::get<0>(freqs[i])[j];
        // }
    }

    // cleaning up.
    for(int i=0; i<VEC_SIZE; i++) {
        for(int j=0; j<T_BITS; j++) {
            delete enc_nums[i][j];
        }
    }

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
    setGlobalVariables(p, r, d, c, w, L, m, gens, ords);
    cout << "Terminat de generat context.\n";

    clock_t begin = clock();
    test_hom_counter();
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
