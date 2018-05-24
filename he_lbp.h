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


/**
 * instead of bootstrapping, we will involve the user in a 
 * "refreshing" the ciphertext after many multiplications.
 * this method is very costly, need improvement.
*/
void pseudoRefreshCtxt(Ctxt *&ctxt);

void test_hom_counter();

void test_LBP();

