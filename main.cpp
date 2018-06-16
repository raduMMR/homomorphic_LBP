#include <iostream>
#include "basic_primitives.h"
#include "he_lbp.h"
#include "HE_FR_LBP.h"
#include <vector>
#include <tuple>
using namespace std;

class DivisionLookupTable {
    public:
    // < <divident, divisor>, quotient>, where quotient = divident/divisor
    vector<tuple<tuple<long, long>, long>> table;
};


void allocCtxt(int ns){

cout << "ns = " << ns << endl;
	int i=0;
	while(true){
		Ctxt* ctxt = encryptBitVal(vector<long>(ns, 1));
		cout << i << endl;
		i++;
	}
}

void test_Ctxt_Add(){
    vector<long> a(NSLOTS, 0);
    vector<long> b(NSLOTS, 1);
    vector<long> result;

    Ctxt* enc_zero = encryptBitVal(a);
    Ctxt* enc_one = encryptBitVal(b);

    Ctxt* sum = encryptBitVal(a);

    sum->addCtxt(*enc_zero);
    result = decryptBitVal(sum);
    cout << " 0 + 0 = " << result[0] << endl;

    (*sum) = (*enc_zero);
    sum->addCtxt(*enc_one);
    result = decryptBitVal(sum);
    cout << " 0 + 1 = " << result[0] << endl;

    (*sum) = (*enc_one);
    sum->addCtxt(*enc_zero);
    result = decryptBitVal(sum);
    cout << " 1 + 0 = " << result[0] << endl;

    (*sum) = (*enc_one);
    sum->addCtxt(*enc_one);
    result = decryptBitVal(sum);
    cout << " 1 + 1 = " << result[0] << endl;

    (*sum) = (*enc_zero);
    sum->addCtxt(*enc_zero, /*negative=*/ true);
    result = decryptBitVal(sum);
    cout << " 0 - 0 = " << result[0] << endl;

    (*sum) = (*enc_zero);
    sum->addCtxt(*enc_one, true);
    result = decryptBitVal(sum);
    cout << " 0 - 1 = " << result[0] << endl;

    (*sum) = (*enc_one);
    sum->addCtxt(*enc_one, true);
    result = decryptBitVal(sum);
    cout << " 1 - 1 = " << result[0] << endl;

    (*sum) = (*enc_one);
    sum->addCtxt(*enc_zero, true);
    result = decryptBitVal(sum);
    cout << " 1 - 0 = " << result[0] << endl;

    delete enc_zero;
    delete enc_one;
    delete sum;
}

int _xor(int a, int b){
    return (a+b)%2;
}

int special_xor(int a, int b){
    return _xor(a,0)*_xor(b,0);
}

int make_number(vector<int> vec){
    int num = 0;
    for(int i=0; i<vec.size(); i++){
        num |= vec[i] << i;
    }
    return num;
}

void test_Difference(){
    vector<int> a(8);
    vector<int> b(8);
    vector<int> result(8, 0);

    bool success = true;

    for(int k=0; k<10000; k++) {
        int num_a = rand() % 256;
        int num_b = rand() % 256;
        if(num_a < num_b) {
            int aux = num_a; num_a = num_b; num_b = aux;
        }
    
        for(int i=0;i<8; i++){
            a[i] = (num_a >> i) & 1;
            b[i] = (num_b >> i) & 1;
        }

        // cout << a << endl << b << endl;
        // cout << num_a << " - " << num_b << " = ";

        int borrowed = 0;
        for(int i=0; i<8; i++){
            int old_a = a[i];
            a[i] = abs(a[i]-borrowed);
            result[i] = abs(a[i] - b[i]);
            borrowed = _xor( 
                _xor( 
                    _xor(a[i], 1) * b[i], 
                    borrowed * _xor(old_a, 1)
                ), 
                special_xor(
                    _xor(a[i], 1)*b[i], 
                    borrowed * _xor(old_a, 1)
                    ) 
                );
        } 
        // cout << make_number(result) << endl;

        if(make_number(result) != (num_a - num_b)){
            cout << "FAIL\n";
            cout << result << endl;
            success = false;
            break;
        }
    }

    if(success == true) {
        cout << "TOTAL SUCCESS\n";
    }
}

void test_mult_time(){
    Ctxt* c1 = encryptBitVal(vector<long>(NSLOTS, 1));
    Ctxt* c2 = encryptBitVal(vector<long>(NSLOTS, 1));
    Ctxt* c3 = encryptBitVal(vector<long>(NSLOTS, 0));

    Ctxt* areEqual = compute_z(0, 1, vector<Ctxt*>(1, c1), vector<Ctxt*>(1, c2));

    cout << decryptBitVal(areEqual)[0] << endl;

    delete areEqual;
    areEqual = compute_z(0, 1, vector<Ctxt*>(1, c1), vector<Ctxt*>(1, c3));
    cout << decryptBitVal(areEqual)[0] << endl;

    delete areEqual;
    delete c3;

    // int i = 0;
    // do {
    //     c1->multiplyBy(*c2);
    //     i++;
    // }
    // while(decryptBitVal(c1)[0] == 1);
    // cout << "Numar de inmultiri = " << i << endl;

    delete c1;
    delete c2;
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

    /************* TESTING SPACE *********************************/
    clock_t begin = clock();
    // test_HE_FR_LBP();
	allocCtxt(s);
    clock_t end = clock();
    cout << "TIMP: " << clock_diff(begin, end) << " secunde.\n";
    /************* TESTING SPACE *********************************/


    cout << "Cleaning up ...\n";
    cleanGlobalVariables();
    cout << "Terminat cleaning up.\n";

    cout << "Program terminat.\n";
    return 0;
}


