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
vector<tuple<vector<Ctxt*>, HE_INT>> hom_counter(const vector<ENC_INT> &enc_nums) {

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
Ctxt* hom_xor(const Ctxt* a, const Ctxt* b){
    Ctxt* result = encryptBitVal(vector<long>(NSLOTS, 0));
    *result = *a;
    result->addCtxt(*b);
    return result;
}

/*************************************************************************************/
Ctxt* special_hom_or(const Ctxt* a, const Ctxt* b){
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


    for(int i=0; i<20; i++){
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
            } else {
		cout << result[j] << ", " << abs(a[j]-b[j]) << endl;
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

/*************************************************************************************/
// for debugging purposes.
void watchCtxt(const Ctxt* ctxt, const char * name, int line){
    cout << name << " : " << decryptBitVal(ctxt)[0] << " on line " << line << endl;
}

/*************************************************************************************/
Ctxt* hom_OR(const Ctxt* a, const Ctxt *b) {
    Ctxt *tmp1 = special_hom_or(a, b);
    Ctxt *tmp2 = hom_xor(a, b);
    Ctxt *result = hom_xor(tmp1, tmp2);
    delete tmp1; delete tmp2;
    return result;
}

/*************************************************************************************/
void hom_binarySum(ENC_INT &accumulator, const ENC_INT &term)
{
    Ctxt* carry = encryptBitVal(vector<long>(NSLOTS, 0));
    Ctxt* old_acc_i = encryptBitVal(vector<long>(NSLOTS, 0));
    // vector of Ctxts aux will be used to hold the following products:
    // aux[0] = carry*a, aux[1]=carry*b, aux[2]=a*b
    // which will be further used to compute the next carry using the formula
    // new_carry = (carry*a)+(carry*b)+(a*b)
    vector<Ctxt*> aux = encryptIntVal(vector<long>(NSLOTS, 0), 3);

    for(int i=0; i<term.size(); i++){
        *old_acc_i = *accumulator[i];
        // watchCtxt(accumulator[i], "accumulator[i]", 508);

        accumulator[i]->addCtxt(*term[i]);
        // watchCtxt(accumulator[i], "accumulator[i]", 511);
        accumulator[i]->addCtxt(*carry);
        // watchCtxt(accumulator[i], "accumulator[i]", 513);
        
        *aux[0] = *carry;
        aux[0]->multiplyBy(*old_acc_i);
        // watchCtxt(aux[0], "aux[0]", 517);
        *aux[1] = *carry;
        aux[1]->multiplyBy(*term[i]);
        // watchCtxt(aux[1], "aux[1]", 520);
        *aux[2] = *old_acc_i;
        aux[2]->multiplyBy(*term[i]);
        // watchCtxt(aux[2], "aux[2]", 523);

        delete carry;
        carry = hom_OR(aux[0], aux[1]);
        // watchCtxt(carry, "carry", 527);
        delete aux[0];
        aux[0] = hom_OR(carry, aux[2]);
        *carry = *aux[0];
        // watchCtxt(carry, "carry", 531);
    }

    pseudoRefreshCtxt(carry);

    for(int i=term.size(); i<accumulator.size(); i++){
        *old_acc_i = *accumulator[i];
        // watchCtxt(accumulator[i], "accumulator[i]", 538);
        accumulator[i]->addCtxt(*carry);
        // watchCtxt(accumulator[i], "accumulator[i]", 540);
        carry->multiplyBy(*old_acc_i);
        // watchCtxt(carry, "carry", 542);
    }

    // cleanup.
    delete carry;
    delete old_acc_i;
    for(int i=0; i<aux.size(); i++){
        delete aux[i];
    }
}

/*************************************************************************************/
void test_homBinarySum(){
    ENC_INT acc = encryptIntVal(vector<long>(NSLOTS, 0), 2*T_BITS);
    int sum = 0;
    int modulo = (int)pow(2, T_BITS)-1;

    for(int i=0; i<10; i++){
        int nr = rand() % modulo;
        cout << nr << " ";
        sum += nr;
        ENC_INT num = encryptIntVal(vector<long>(NSLOTS, nr), T_BITS);
        // cout << "Before hom_binarySum: " << decryptIntVal(acc)[0] << endl;
        hom_binarySum(acc, num);
        for(int j=0; j<acc.size(); j++){
            pseudoRefreshCtxt(acc[j]);
        }
        // cout << "After hom_binarySum: " << decryptIntVal(acc)[0] << endl;
        // cleanup
        for(int j=0; j<num.size(); j++) {
            delete num[j];
        }
    }
    cout << endl << "REGULAR SUM = " << sum << endl;
    cout << "SUM = " << decryptIntVal(acc)[0] << endl;

    // cleanup.
    for(int i=0; i<acc.size(); i++){
        delete acc[i];
    }
}

void refresh_ENC_INT(ENC_INT &enc_int) {
    for(int i=0; i<enc_int.size(); i++){
        pseudoRefreshCtxt(enc_int[i]);
    }
}

/*************************************************************************************/
ENC_INT absoluteValueMetric(ENC_HIST &h1, ENC_HIST &h2){
    ENC_INT sum = encryptIntVal(vector<long>(NSLOTS, 0), 2*T_BITS);
    ENC_INT tmp_diff;
    ENC_INT enc_zero = encryptIntVal(vector<long>(NSLOTS, 0), T_BITS);

    // ?? assert(h1.size() == h2.size());
    for(int i=0; i<h1.size(); i++) {
        for(int j=0; j<h2.size(); j++) {

            // refresh the ctxts before proceeding with further multiplications.
            refresh_ENC_INT(std::get<1>(h1[i]).enc_bits);
            refresh_ENC_INT(std::get<1>(h2[i]).enc_bits);

            tmp_diff = absDifference(std::get<1>(h1[i]).getNumber(), std::get<1>(h2[j]).getNumber());

            refresh_ENC_INT(tmp_diff);

            refresh_ENC_INT(std::get<0>(h1[i]));
            refresh_ENC_INT(std::get<0>(h2[i]));

            Ctxt *areEqual = compute_z(0, T_BITS, std::get<0>(h1[i]), std::get<0>(h2[j]));
		pseudoRefreshCtxt(areEqual);
            for(int k=0; k<tmp_diff.size(); k++){
                tmp_diff[k]->multiplyBy(*areEqual);
            }
// cout << "Dupa mult cu equal, tmp = " << decryptIntVal(tmp_diff)[0] << endl;
// cout << "Pentru nr: " << decryptIntVal(std::get<0>(h1[i]))[0] << ", " << decryptIntVal(std::get<0>(h2[j]))[0] << endl;

// and if the frequency is not 0
Ctxt* f1_is_zero = compute_z(0, T_BITS, std::get<1>(h1[i]).getNumber(), enc_zero);
Ctxt* f2_is_zero = compute_z(0, T_BITS, std::get<1>(h2[j]).getNumber(), enc_zero);
// cout << "inainte de f1*f2\n";
pseudoRefreshCtxt(f1_is_zero);
pseudoRefreshCtxt(f2_is_zero);
f1_is_zero->multiplyBy(*f2_is_zero);
// cout << "dupa f1*f2\n";
for(int k=0; k<tmp_diff.size(); k++){
	tmp_diff[k]->multiplyBy(*f1_is_zero);
}
// cout << "dupa tmp_diff * f1 *f2\n";
            hom_binarySum(sum, tmp_diff);
//	    cout << "sum = " << decryptIntVal(sum)[0] << endl;
            // cleanup.
delete f1_is_zero; delete f2_is_zero;
            delete areEqual;
            for(int k=0; k<tmp_diff.size(); k++){
                delete tmp_diff[k];
            }

            // cout << "Step i = " << i << ", j = " << j << endl;
        }
    }

// another cleanup
for(int i=0; i<enc_zero.size(); i++){
	delete enc_zero[i];
}

//	cout << "Before return sum\n";

    return sum;
}

/*************************************************************************************/
void test_absoluteValueMetric()
{
    vector<long> myvector(NSLOTS, 0);
    for(int i=0; i<NSLOTS; i++) {
        myvector[i] = rand() % 3;
    }

    int VEC_SIZE = 5;

    vector<ENC_INT> enc_nums(VEC_SIZE);
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
    }

    // compare the two identical frequencies.
    ENC_INT decision = absoluteValueMetric(freqs, freqs);
    cout << "Sum of abs diffs = " << decryptIntVal(decision)[0] << endl;

    // cleanup.
cout << "Stergere decision ...\n";
    for(int i=0; i<decision.size(); i++) {
        delete decision[i];
    }

    // Now, compare with a different histogram.
/*    for(int i=0; i<NSLOTS; i++) {
        myvector[i] = rand() % 3;
    }

    vector<ENC_INT> enc_nums2(VEC_SIZE);
    for(int i=0; i<VEC_SIZE; i++) {
        cout << myvector[0] << " ";
        enc_nums2[i] = encryptIntVal(myvector, T_BITS);
        random_shuffle ( myvector.begin(), myvector.end(), myrandom);
    }
    cout << endl;

    vector<tuple<vector<Ctxt*>, HE_INT>> freqs2 = hom_counter(enc_nums2);

    cout << "freqs2.size() = " << freqs2.size() << endl;
    for(int i=0; i<freqs2.size(); i++) {
        vector<long> number = decryptIntVal(std::get<0>(freqs2[i]));
        vector<long> frecventa = decryptIntVal(std::get<1>(freqs2[i]).enc_bits);
        cout << "Numar = " << number[0] << ", Frecventa =  " << frecventa[0] << endl;
    }

    // compare the two identical frequencies.
    ENC_INT decision2 = absoluteValueMetric(freqs, freqs2);
    cout << "Sum of abs diffs = " << decryptIntVal(decision2)[0] << endl;

    // cleaning up.
cout << "Stergere decision2\n";
    for(int i=0; i<decision.size(); i++) {
        delete decision2[i];
    }*/
// cout << "Stergere deciosion 2 completex\n";
cout << "Stergere enc_nums ...\n";
    for(int i=0; i<enc_nums.size(); i++) {
        for(int j=0; j<enc_nums[i].size(); j++) {
	cout<< "deleting enc_nums[" << i << "][" << j << "]\n";     
       delete enc_nums[i][j];
        }
    }
cout << "Stergere enc_nums completed.\n";
/*cout << "Stergere enc_nums2 ...\n";
    for(int i=0; i<VEC_SIZE; i++) {
        for(int j=0; j<T_BITS; j++) {
            delete enc_nums2[i][j];
        }
    }
cout << "Stergere enc_nums2 completed.\n";*/
    // BUG!! cleaning this memory results in corruption, because it was cleared previously
    // frequency contain only pointer to enc_nums !!!!
    // for(int i=0; i<freqs.size(); i++){
    //     for(int j=0; j<std::get<0>(freqs[i]).size(); j++){
    //         delete std::get<0>(freqs[i])[j];
    //     }
    // }
    // for(int i=0; i<freqs2.size(); i++){
    //     for(int j=0; j<std::get<0>(freqs2[i]).size(); j++){
    //         delete std::get<0>(freqs2[i])[j];
    //     }
    // }

}
