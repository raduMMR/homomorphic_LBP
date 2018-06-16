// HE_FR_LBP Homomorphic Encryption for Face Recognition using LBP operator
#pragma once
#include <iostream>
#include <vector>
#include "./../FHE.h"
using namespace std;


// we choose to make EncRegion and NumberOfOccurences classes to use 
// the destructor method of the objects. TO DO: protect more the fields of these 
// two classes
class EncRegion
{
    int nslots;

    void cloneRegion(const EncRegion& copy);

public:
    vector<Ctxt*> enc_pixels;
    vector<vector<Ctxt*>> enc_neighbours;

    EncRegion(int nslots){
        this->nslots = nslots;
    }

    EncRegion(const EncRegion &copy){
        cloneRegion(copy);
    }

    EncRegion& operator=(const EncRegion& copy){
        cloneRegion(copy);
    }

    ~EncRegion();
};

// encrypted number of occurences of the lbp_code.
class NumberOfOccurences{
public:

    // only reference. DO NOT FREE this vector.
    vector<Ctxt*> lbp_code;

    Ctxt* occurence_hits;

    ~NumberOfOccurences();
};

class EncHistogram{
public:
    EncHistogram();

    // vector<Ctxt*> contains a vector of tuples ( enc(lbp_code), enc(nb_of_occurences) )
    vector<NumberOfOccurences> enc_hist;

    // for each lbp code from enc_hist compute the number of occurences.
    // @param lbp_codes is matrix
    /*          bit_0 bit_1  ... bit_7
    lbp_codes = [                    ]  slot_0    => LBP code for pixel_0 
                        ...             slot_1
                        ...              ...
                [                    ]  slot_4096  => LBP code for pixel_4096
    */
    void computeNbOfOccurences(vector<Ctxt*> lbp_codes);

    ~EncHistogram();
};

class ImageProcessor{

    int nslots;

    vector<EncRegion*> encrypted_regions;

    void encodeImage(uint8_t **image);

    vector<Ctxt*> computeLBP4Region(int region_id);

    EncRegion *encryptRegion(vector<long> pixels, vector<vector<long>> neighbours);

public:

    void setNumberOfSlots(int nslots){
        this->nslots = nslots;
    }

    void encryptImage(uint8_t **image);

    EncHistogram computeHistogram4Region(int region_id);

    ~ImageProcessor(){
        for(int i=0; i<encrypted_regions.size(); i++){
            delete encrypted_regions[i];
        }
    }
};


void test_HE_FR_LBP();