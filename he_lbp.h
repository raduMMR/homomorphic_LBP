#pragma once
#include "basic_primitives.h"

vector<Ctxt*> hom_LBP(vector<Ctxt*> enc_pixeli, vector<vector<Ctxt*>> vecini, int t_bits);

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

    ~HE_INT();

};

typedef vector<Ctxt*> ENC_INT;

/**
 * @param enc_nums - vector of encrypted numbers on bits representation
 * @return - the number of occurences for each number in the vector
 */
vector<tuple<vector<Ctxt*>, HE_INT>> hom_counter(vector<ENC_INT> &enc_nums);


// returns the absolute value of the difference a-b, |a-b|
ENC_INT absDifference(const ENC_INT &a, const ENC_INT &b);

Ctxt* hom_xor(Ctxt* a, Ctxt* b);

// implements OR bit operation using XOR provided by ctxt addition.
Ctxt* special_hom_or(Ctxt* a, Ctxt* b);

// "sorts" the numbers a & b, putting in gt the greater than number
// and in lt the less than number.
void orderNumbers(const vector<Ctxt*> &a, const vector<Ctxt*> &b, 
                    vector<Ctxt*> &gt, vector<Ctxt*> &lt);


/**
 * instead of bootstrapping, we will involve the user in a 
 * "refreshing" the ciphertext after many multiplications.
 * this method is very costly, need improvement.
*/
void pseudoRefreshCtxt(Ctxt *&ctxt);

void test_hom_counter();

void test_LBP();

void test_absDifference();









