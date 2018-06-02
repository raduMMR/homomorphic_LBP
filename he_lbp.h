#pragma once
#include "basic_primitives.h"

typedef vector<Ctxt*> ENC_INT;

class HE_INT {

public:

    vector<Ctxt*> enc_bits;

    HE_INT();

    HE_INT(const HE_INT& he_int);

    HE_INT& operator=(const HE_INT& he_int);

    void clone(const HE_INT& he_int);

    // make HE_INT an encrypted representation of zero,
    // with enc_bits[0] = LSB and enc_bits[T_BITS-1] = MSB.
    void initialize_HE_INT();
    
    void add1Bit(Ctxt* ctxt);

    void setSkipping(const Ctxt &isSkipped);

    ENC_INT getNumber() const { return enc_bits; }

    ~HE_INT();

};

// typedef for encrypted histogram, the number of occurences for 
// each LBP code. 
// vector<Ctxt*> - the LBP code
// HE_INT - the number of occurences
typedef vector<tuple<vector<Ctxt*>, HE_INT>> ENC_HIST;

vector<Ctxt*> hom_LBP(vector<Ctxt*> enc_pixeli, vector<vector<Ctxt*>> vecini, int t_bits);

/**
 * @param enc_nums - vector of encrypted numbers on bits representation
 * @return - the number of occurences for each number in the vector
 */
vector<tuple<vector<Ctxt*>, HE_INT>> hom_counter(vector<ENC_INT> &enc_nums);


// returns the absolute value of the difference a-b, |a-b|
ENC_INT absDifference(const ENC_INT &a, const ENC_INT &b);

Ctxt* hom_xor(const Ctxt* a, const Ctxt* b);

// implements a part of the OR bit operation using XOR provided by ctxt addition.
// a || b = a x b x ((a x 0) * (b x 0)). This function implements only this part ((a x 0) * (b x 0)).
Ctxt* special_hom_or(const Ctxt* a, const Ctxt* b);

// the homomorphic OR implementation  using XOR provided by ctxt addition.
Ctxt* hom_OR(const Ctxt* a, const Ctxt *b);

// "sorts" the numbers a & b, putting in gt the greater than number
// and in lt the less than number.
void orderNumbers(const vector<Ctxt*> &a, const vector<Ctxt*> &b, 
                    vector<Ctxt*> &gt, vector<Ctxt*> &lt);

// we'll try to use the absolute value metric for comparing histograms.
// @param h1, h2 - the histograms to be compared
// @return the value of the absolute metric used for decision.
ENC_INT absoluteValueMetric(const ENC_HIST &h1, const ENC_HIST &h2);


// add to accumulator the term number
// numbers accumulator and term are in binary format.
void hom_binarySum(ENC_INT &accumulator, const ENC_INT &term);

/**
 * instead of bootstrapping, we will involve the user in a 
 * "refreshing" the ciphertext after many multiplications.
 * this method is very costly, need improvement.
*/
void pseudoRefreshCtxt(Ctxt *&ctxt);



// tests
void test_hom_counter();

void test_LBP();

void test_absDifference();

void test_absoluteValueMetric();

void test_homBinarySum();









