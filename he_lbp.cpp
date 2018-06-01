#include <tuple>
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include "./he_lbp.h"

HE_INT::HE_INT() {
    this->initialize_HE_INT();
}

HE_INT::HE_INT(const HE_INT& he_int) {
    this->clone(he_int);
}

HE_INT& HE_INT::operator=(const HE_INT& he_int) {
    this->clone(he_int);
    return (*this);
}

void HE_INT::clone(const HE_INT& he_int) {
    this->initialize_HE_INT();
    for(int i=0; i<enc_bits.size(); i++) {
        (*this->enc_bits[i]) = *(he_int.enc_bits[i]);
    }
}

void HE_INT::initialize_HE_INT() {
    vector<long> zeroes(NSLOTS, 0);
    this->enc_bits = encryptIntVal(zeroes, T_BITS);
}
    
void HE_INT::add1Bit(Ctxt* ctxt) {
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
        pseudoRefreshCtxt(enc_bits[i]);
        *carry = *new_carry;
    }

    delete carry;
    delete new_carry;
}

void HE_INT::setSkipping(const Ctxt &isSkipped) {
    for(int i=0; i<enc_bits.size(); i++) {
        enc_bits[i]->multiplyBy(isSkipped);
    }
}

HE_INT::~HE_INT() {
    for(int i=0; i<enc_bits.size(); i++) {
        delete enc_bits[i];
    }
}

/*************************************************************************************/
void pseudoRefreshCtxt(Ctxt *&ctxt) {
    vector<long> bits;
    bits =  decryptBitVal(ctxt);
    delete ctxt;
    ctxt = encryptBitVal(bits);
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

    // encryption of one, used later.
    Ctxt *ENC_ONE = encryptBitVal(vector<long>(NSLOTS, 1));
    Ctxt *aux = encryptBitVal(vector<long>(NSLOTS, 1));

    // cout << "enc_nums.size() = " << enc_nums.size() << endl;
    ZZ z_one(1);
    for(int i=0; i<enc_nums.size(); i++) {
        HE_INT frequency;

        (*isSkipped) = (*isCounted[i]);

        // isSkipped->negate(); // for multiplication purpose, see below.
        isSkipped->addConstant(z_one);
        // cout << "~isSkipped = " << decryptBitVal(isSkipped)[0] << endl;
        // if the number has been previously hit, then multiply with 0
        // the frequency of the counting.

        for(int j=i; j<enc_nums.size(); j++) {

            Ctxt* areEqual = compute_z(0, T_BITS, enc_nums[i], enc_nums[j]);

            pseudoRefreshCtxt(areEqual);

            frequency.add1Bit(areEqual);

            *aux = *ENC_ONE;

            // mark number as counted, only if it has been hit.
            // (*isCounted[j]) = (*areEqual); 
            ENC_ONE->addCtxt(*areEqual, /*negative=*/true); 
            isCounted[j]->multiplyBy(*ENC_ONE);
            isCounted[j]->addCtxt(*areEqual);

            *ENC_ONE = *aux;

            delete areEqual;
        }
        // cout << "before Skipping, f = " << decryptIntVal(frequency.enc_bits)[0] << endl;

        // skip number if already counted.
        // that means to set the frequency to zero for another hit
        // of a previous number. When doing chi-square computation
        // simply do not take into account elements with zero frequency.
        frequency.setSkipping(*isSkipped);

        // cout << "after Skipping, f = " << decryptIntVal(frequency.enc_bits)[0] << endl;

        frequencies.push_back(make_tuple(enc_nums[i], frequency));
    }

    // cleaning up.
    delete isSkipped;
    delete ENC_ONE;
    delete aux;

    for(int i=0; i<isCounted.size(); i++) {
        delete isCounted[i];
    }

    return frequencies;
}

/*************************************************************************************/
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

/*************************************************************************************/
int myrandom (int i) { return std::rand()%i;}

/*************************************************************************************/
void test_hom_counter() {
    // generate a random vector cu NSLOTS values.
    vector<long> myvector(NSLOTS, 0);
    for(int i=0; i<NSLOTS; i++) {
        myvector[i] = rand() % 3;
    }
    // myvector is going to be encrypted in one vector<Ctxt*>
    // so we'll shuffle this vector to have different values
    // in the corresponding slots.

    int VEC_SIZE = 5;

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

/*************************************************************************************/
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

/*************************************************************************************/
Ctxt* hom_xor(Ctxt* a, Ctxt* b){
    Ctxt* result = encryptBitVal(vector<long>(NSLOTS, 0));
    *result = *a;
    result->addCtxt(*b);
    return result;
}

/*************************************************************************************/
Ctxt* special_hom_or(Ctxt* a, Ctxt* b){
    Ctxt* enc_zero = encryptBitVal(vector<long>(NSLOTS, 0));
    Ctxt* result = hom_xor(a, enc_zero);
    Ctxt* op2 = hom_xor(b, enc_zero);
    result->multiplyBy(*op2);
    delete op2;
    return result;
}

/*************************************************************************************/
void orderNumbers(const vector<Ctxt*> &a, const vector<Ctxt*> &b, 
                    vector<Ctxt*> &gt, vector<Ctxt*> &lt)
{
    Ctxt *t = compute_s(0, a.size(), a, b);
    pseudoRefreshCtxt(t);

    // this will be the greater than number between a and b.
    gt = encryptIntVal(vector<long>(NSLOTS, 0), a.size());

    // this will be the less than number between a and b
    lt = encryptIntVal(vector<long>(NSLOTS, 0), a.size());

    Ctxt* old_t = encryptBitVal(vector<long>(NSLOTS, 0));
    *old_t = *t;
    Ctxt* enc_one = encryptBitVal(vector<long>(NSLOTS, 1));

    for(int i=0; i<a.size(); i++){
        *gt[i] = (*enc_one);
        gt[i]->addCtxt(*t, /*negative=*/true);
        gt[i]->multiplyBy(*b[i]);
        t->multiplyBy(*a[i]);
        gt[i]->addCtxt(*t);
        (*t) = *old_t;

        *lt[i] = (*enc_one);
        lt[i]->addCtxt(*t, /*negative=*/true);
        lt[i]->multiplyBy(*a[i]);
        t->multiplyBy(*b[i]);
        lt[i]->addCtxt(*t);
        (*t) = *old_t;
    }

    // cleanup.
    delete old_t;
    delete t;
    delete enc_one;
}

/*************************************************************************************/
ENC_INT absDifference(const ENC_INT &a, const ENC_INT &b)
{
    assert(a.size() == b.size());
    int num_size = a.size();

    ENC_INT absDiff = encryptIntVal(vector<long>(NSLOTS, 0), num_size);
    
    vector<Ctxt*> gt, lt;
    orderNumbers(a, b, gt, lt);

    Ctxt* old_gt_i = encryptBitVal(vector<long>(NSLOTS, 0));
    Ctxt* enc_one = encryptBitVal(vector<long>(NSLOTS, 1));
    Ctxt* borrowed = encryptBitVal(vector<long>(NSLOTS, 0));

    vector<Ctxt*> aux(5);

    for(int i=0; i<num_size; i++){
        *old_gt_i = *gt[i];
        // a[i] = abs(a[i]-borrowed);
        gt[i]->addCtxt(*borrowed, /*negative=*/true); 
        // result[i] = abs(a[i] - b[i]);
        *absDiff[i] = *gt[i];
        absDiff[i]->addCtxt(*lt[i], /*negative=*/true); 
        // _xor(a[i], 1) * b[i]
        aux[0] = hom_xor(gt[i], enc_one);
        aux[0]->multiplyBy(*lt[i]);
        // borrowed * _xor(old_a, 1)
        aux[1] = hom_xor(old_gt_i, enc_one);
        aux[1]->multiplyBy(*borrowed);

        aux[2] = hom_xor(aux[0], aux[1]);
        aux[3] = special_hom_or(aux[0], aux[1]);
        aux[4] = hom_xor(aux[2], aux[3]);
        *borrowed = *aux[4];

        for(int i=0; i<5; i++){ delete aux[i]; }

        // borrowed = _xor( 
        //     _xor( 
        //         _xor(a[i], 1) * b[i], 
        //         borrowed * _xor(old_a, 1)
        //     ), 
        //     special_xor(
        //         _xor(a[i], 1)*b[i], 
        //         borrowed * _xor(old_a, 1)
        //      ) 
        //     );
    } 

    delete old_gt_i; delete enc_one; delete borrowed;
    assert(gt.size() == lt.size() );
    for(int i=0; i<gt.size(); i++){
        delete gt[i];
        delete lt[i];
    }

    return absDiff;
}

/*************************************************************************************/
void test_absDifference(){
    vector<long> a(NSLOTS, 0);
    vector<long> b(NSLOTS, 0);


    for(int i=0; i<10; i++){
        for(int j=0; j<NSLOTS; j++){
            a[j] = random()%256;
            b[j] = random()%256;
        }

        ENC_INT enc_a = encryptIntVal(a, 8);
        ENC_INT enc_b = encryptIntVal(b, 8);

        ENC_INT diff = absDifference(enc_a, enc_b);

        vector<long> result = decryptIntVal(diff);

        bool success = true;

        for(int j=0; j<NSLOTS; j++) {
            if(result[j] != abs(a[j]-b[j])){
                cout << "FAIL\n";
                cout << "|" << a[j] << " - " << b[j] << "| = " << result[j] << endl;
                success = false;
                break;
            }
        }

        // cleanup.
        for(int j=0; j<8; j++){
            delete enc_a[j];
            delete enc_b[j];
            delete diff[j];
        }

        if(success == false){
            break;
        }
    }
}
