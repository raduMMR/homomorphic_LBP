// HE_FR_LBP Homomorphic Encryption for Face Recognition using LBP operator
#pragma once
#include <iostream>
#include <vector>
#include "./../FHE.h"
using namespace std;

#define NUMBER_OF_REGIONS 256 // 128 // 64 // the number of regions the image is divided by
#define NUMBER_OF_PIXELS 256 //512 // 1024 // the number of pixels in each region
#define REGION_ROOT 16 // 

class ImageProcessor{

    int nslots;

    void encodeImage(uint8_t **image);

    vector<Ctxt*> computeLBP4Region(char *regionFile);

    void encryptRegion(vector<long> pixels, vector<vector<long>> neighbours, char *regionFile);

    void computeNbOfOccurences(vector<Ctxt*> LBP_codes, vector<long> histogram, char *filename);

public:

    void computeHistogram4Region(char *regionFile, char *regionHist);

    void setNumberOfSlots(int nslots){ this->nslots = nslots; }

    void encryptImage(uint8_t **image, char **histFiles);

    void computeHistogram(char **regionFiles, char **histFiles);
};

// compare two encrypted histograms of two images
// for each region of an image, the local histogram is saved in one file.
vector<Ctxt*> compareTwoImages(char **histogram1Files, char **histogram2Files, int nb_regions);

void test_HE_FR_LBP();

void myitoa(int number, char *buffer);

