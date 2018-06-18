// HE_FR_LBP Homomorphic Encryption for Face Recognition using LBP operator
#include "HE_FR_LBP.h"
#include "basic_primitives.h"
#include "he_lbp.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>


void add1Bit(vector<Ctxt*> &metric, Ctxt* &ctxt) {
    Ctxt *carry = encryptBitVal(vector<long>(NSLOTS, 0));
    *carry = (*metric[0]);
    carry->multiplyBy((*ctxt));

    metric[0]->addCtxt(*ctxt);

    Ctxt *new_carry = encryptBitVal(vector<long>(NSLOTS, 0));
    for(int i=1; i<metric.size(); i++) {
        *new_carry = (*metric[i]);
        new_carry->multiplyBy(*carry);
        metric[i]->addCtxt(*carry);

        if(needsBootstrapping(metric[i]))
            pseudoRefreshCtxt(metric[i]);

        *carry = *new_carry;
    }

    delete carry;
    delete new_carry;
}

vector<Ctxt*> computeRegionMetric(char *histFile1, char *histFile2){
    fstream hist1(histFile1, fstream::in);
    fstream hist2(histFile2, fstream::in);

    // 16 is chosen arbitrarly, hope it will suffice and at the same time is not too big.
    vector<Ctxt*> metric(16, encryptBitVal(vector<long>(NUMBER_OF_PIXELS, 0)));

    Ctxt* c1 = encryptBitVal(vector<long>(NUMBER_OF_PIXELS, 0));
    Ctxt* c2 = encryptBitVal(vector<long>(NUMBER_OF_PIXELS, 0));

    for(int i=0; i<256; i++){
        hist1 >> *c1;
        hist2 >> *c2;

        // compare the two amplitudes of the same LBP code.
        c2->addCtxt(*c1);
        // adding one amplitude to another, we obtain the number of differences


        // now add the differences to the global metric.
        add1Bit(metric, c2);
    }

    delete c1; delete c2;

    hist1.close();
    hist2.close();

    return metric;
}

vector<Ctxt*> getRow(const vector<Ctxt*> &global_metric, const int slot_index){
    vector<long> plain_mask(NUMBER_OF_PIXELS, 0);
    plain_mask[slot_index] = 1;
    // we will have the value 1 in slot 1 of each Ctxt so that I can obtain
    // each row of Ctxts and add to the final metric.
    vector<Ctxt*> row(global_metric.size(), encryptBitVal(plain_mask));

    // I will zeroes everywhere except in every slot_index of the Ctxt.
    for(int i=0; i<global_metric.size(); i++){
        row[i]->multiplyBy(*global_metric[i]);
    }
    
    // shift the row in the front position to add the returned value
    // to the final metric.
    for(int i=0; i<row.size(); i++){
        // shift Ctxt toward last slot.
        ea->shift(*row[i], NUMBER_OF_PIXELS-slot_index-1);
    }

    return row;
}

void addRow(vector<Ctxt*> &global_metric, vector<Ctxt*> &row){
    for(int i=0; i<row.size(); i++){
        add1Bit(global_metric, row[i]);
    }
}

vector<Ctxt*> compareTwoImages(char **histogramFiles1, char **histogramFiles2, int nb_regions){
    vector<Ctxt*> global_metric = encryptIntVal(vector<long>(NUMBER_OF_PIXELS, 0), 16);

    for(int i=0; i<nb_regions; i++){
        vector<Ctxt*> region_metric = computeRegionMetric(histogramFiles1[i], histogramFiles2[i]);
        hom_binarySum(global_metric, region_metric);
    }

    // accumulate in the last slot the global metric from partial metric of the region metrics
    
    // get all rows, except last one where the accumulation will take place.
    for(int i=0; i<NUMBER_OF_PIXELS-1; i++){
        vector<Ctxt*> row = getRow(global_metric, i);
        addRow(global_metric, row);
        for(int j=0; j<row.size(); j++){
            delete row[j];
        }
    }

    return global_metric;
}

void ImageProcessor::computeNbOfOccurences(vector<Ctxt*> LBP_codes, vector<long> histogram, char *filename){
    // save the histogram of this region in this file
    fstream localHistFile(filename, fstream::out|fstream::trunc);
    assert(localHistFile.is_open());

    for(int i=0; i<256; i++){
        vector<Ctxt*> lbp_code = encryptIntVal(vector<long>(NUMBER_OF_PIXELS,/*LBP_code=*/ i), /*T_BITS=*/8);

        Ctxt *eh = compute_z(0, 8, lbp_code, LBP_codes);

        for(int j=0; j<lbp_code.size(); j++){
            delete lbp_code[j];
        }

        localHistFile << *eh << endl;
        
        // testing with the plain histogram.
        // vector<long> dec_hist = decryptBitVal(eh);
        delete eh;

        // int hist_val = 0;
        // for(int j=0; j<NUMBER_OF_PIXELS; j++){
        //     hist_val += dec_hist[j];
        // }
        // if(hist_val != histogram[i]){
        //     cout << "Hist diff " << i << endl;
        // } else { cout << "Step " << i << " completed.\n"; }
    }

    localHistFile.close();
}

void ImageProcessor::encryptRegion(vector<long> pixels, vector<vector<long>> neighbours, char *regionFile){
    assert(neighbours.size() == 8);

    fstream regionStream(regionFile, fstream::out|fstream::trunc);

    // EncRegion* er = new EncRegion(this->nslots);
    // er->enc_pixels = encryptIntVal(pixels, 8);
    // er->enc_neighbours = vector<vector<Ctxt*>>(8);
    // for(int i=0; i<8; i++) {
    //     er->enc_neighbours[i] = encryptIntVal(neighbours[i], 8);
    // }

    vector<Ctxt*> enc_pixels = encryptIntVal(pixels, 8);
    for(int i=0; i<enc_pixels.size(); i++){
        regionStream << *enc_pixels[i] << endl;
        delete enc_pixels[i];
    }

    for(int i=0; i<8; i++){
        vector<Ctxt*> enc_neighbour = encryptIntVal(neighbours[i], 8);
        for(int j=0; j<enc_neighbour.size(); j++){
            regionStream << *enc_neighbour[j] << endl;
            delete enc_neighbour[j];
        }
    }

    regionStream.close();
    regionStream.flush();
}

// compute the LBP codes in the clear for comparison.
vector<long> clearLBPcodes(vector<long> pixels, vector<vector<long>> neighbours){
    vector<long> lbp_codes(NUMBER_OF_PIXELS);
    for(int i=0; i<NUMBER_OF_PIXELS; i++) {
        lbp_codes[i] = 0;
        for(int j=0; j<8; j++) {
            lbp_codes[i] |= (neighbours[j][i] >= pixels[i]) << j;
        }
    }
    return lbp_codes;
}

// compute the histogram of a region in the clear for testing purpose.
vector<long> computeClearHistogram(vector<long> lbp_codes){
    vector<long> histogram(256);

    for(int i=0; i<256; i++){
        int counter = 0;
        for(int j=0; j<lbp_codes.size(); j++){
            if(lbp_codes[j] == i){
                counter++;
            }
        }
        histogram[i] = counter;
    }

    return histogram;
}

void ImageProcessor::encryptImage(uint8_t **image, char **histFiles){
    assert(this->nslots == NUMBER_OF_PIXELS);
    // for an image of size 256x256 we split it in 64 regions, 4x4
    // such that a ctxt contain the pixels from an entire region, 1024

    // we added two extra columns and two extra lines to border 
    // the image with zeroes such that every pixel has 8 neighbours
    for(int i=0; i<8; i++){
        for(int j=0; j<8; j++){
            vector<long> pixels(NUMBER_OF_PIXELS);
            int k=0;
            vector<vector<long>> neighbours(8, vector<long>(NUMBER_OF_PIXELS));
            for(int line=1+256/8*i; line<256/8*(i+1)+1; line++){
                for(int col=1+256/8*j; col<256/8*(j+1)+1; col++){
                    pixels[k] = image[line][col];

                    // the neighbours of the pixel (i,j)
                    neighbours[0][k] = image[line-1][col-1];
                    neighbours[1][k] = image[line][col-1];
                    neighbours[2][k] = image[line+1][col-1];
                    neighbours[3][k] = image[line+1][col];
                    neighbours[4][k] = image[line+1][col+1];
                    neighbours[5][k] = image[line][col+1];
                    neighbours[6][k] = image[line-1][col+1];
                    neighbours[7][k] = image[line-1][col+1];

                    k++;
                }
            }  

            // encrypt the region
            cout << "se cripteaza regiunea (" << i << ", " << j << ")\n";
            encryptRegion(pixels, neighbours, "region.enc");      

            // in production, the next call could not be here, but for speeding
            // testing is good
            this->computeHistogram4Region("region.enc", histFiles[i*8+j]);
            cout << "Regiune criptata\n";
        }
    }
}

vector<Ctxt*> ImageProcessor::computeLBP4Region(char *regionFile){
    fstream regionStream(regionFile, fstream::in);
    vector<Ctxt*> enc_pixels = encryptIntVal(vector<long>(NUMBER_OF_PIXELS, 0), 8);
    vector<vector<Ctxt*>> enc_neighbours(8, enc_pixels);

    for(int i=0; i<8; i++){
        regionStream >> *enc_pixels[i];
    }

    for(int i=0; i<8; i++){
        for(int j=0; j<8; j++){
            regionStream >> *enc_neighbours[i][j];
        }
    }

    regionStream.close();

    vector<Ctxt*> enc_lbp = hom_LBP(enc_pixels, enc_neighbours, 8);

    // cleanup.
    for(int i=0; i<enc_pixels.size(); i++){
        delete enc_pixels[i];
    }
    for(int i=0; i<enc_neighbours.size(); i++){
        for(int j=0; j<enc_neighbours[i].size(); j++){
            delete enc_neighbours[i][j];
        }
    }

    return enc_lbp;
}

void ImageProcessor::computeHistogram4Region(char *regionFile, char *regionHist){
    vector<Ctxt*> LBP_codes = computeLBP4Region(regionFile);

    computeNbOfOccurences(LBP_codes, vector<long>(NUMBER_OF_PIXELS, 0), regionHist);
}

void ImageProcessor::computeHistogram(char **regionFiles,  char **histFiles){

    for(int i=0; i<NUMBER_OF_REGIONS; i++){
        this->computeHistogram4Region(regionFiles[i], histFiles[i]);
    }
    
}

void myitoa(int number, char *buffer){
    strcpy(buffer, "");

    vector<long> digits;
    do{
        digits.push_back(number%10);
        number = number/10;
    }while(number!=0);

    for(int i=digits.size()-1; i>=0; i--){
        switch(digits[i]){
            case 0: 
                strcat(buffer, "0"); break;
            case 1: 
                strcat(buffer, "1"); break;
            case 2: 
                strcat(buffer, "2"); break;
            case 3: 
                strcat(buffer, "3"); break;
            case 4: 
                strcat(buffer, "4"); break;
            case 5: 
                strcat(buffer, "5"); break;
            case 6: 
                strcat(buffer, "6"); break;
            case 7: 
                strcat(buffer, "7"); break;
            case 8: 
                strcat(buffer, "8"); break;
            case 9: 
                strcat(buffer, "9"); break;
            default: 
                cout << "Can't be! Digit not existing ??? \n"; 
                exit(-1);
        }
    }
}

void test_HE_FR_LBP(){
    uint8_t **image = new uint8_t*[258];
    for(int i=0; i<258; i++){
        image[i] = new uint8_t[258];
        for(int j=0; j<258; j++){
            image[i][j] = rand() % 256;
        }
    }

    // border the image with zeroes such that each pixel would have eight neighbours.
    for(int i=0; i<258; i++){
        image[i][0] = 0;
        image[257][i] = 0;
        image[i][257] = 0;
        image[0][i] = 0;
    }

    char **regionFiles = new char*[NUMBER_OF_REGIONS];
    char **histFiles = new char*[NUMBER_OF_REGIONS];
    for(int i=0; i<NUMBER_OF_REGIONS; i++){
        regionFiles[i] = new char[7];
        myitoa(i, regionFiles[i]);
        strcat(regionFiles[i], ".enc");

        histFiles[i] = new char[8];
        myitoa(i, histFiles[i]);
        strcat(histFiles[i], ".hist");
    }

    ImageProcessor ip;
    ip.setNumberOfSlots(NUMBER_OF_PIXELS);

    // cout << "Se cripteaza imaginea ...\n";
    ip.encryptImage(image, histFiles);
    // cout << "Imagine criptata.\n";

    // cout << "Se calculeaza histograma imaginii ...\n";
    // ip.computeHistogram(regionFiles, histFiles);
    // cout << "Histograma calculata.\n";


    cout << "Se compara doua histograme identice ...\n";
    vector<Ctxt*> result = compareTwoImages(histFiles, histFiles, NUMBER_OF_REGIONS);
    cout << "Terminat de comparat histogramele.\n";

    vector<long> decision = decryptIntVal(result);
    cout << "Metrica comparatiei histogramelor: " << decision[NUMBER_OF_PIXELS-1] << endl;

    // cleanup, delete image matrix.
    for(int i=0; i<result.size(); i++){
        delete result[i];
    }

    for(int i=0; i<258; i++){
        delete image[i];
    }
    delete image;

    for(int i=0; i<NUMBER_OF_REGIONS; i++){
        delete regionFiles[i];
        delete histFiles[i];
    }
    delete regionFiles;
    delete histFiles;
}