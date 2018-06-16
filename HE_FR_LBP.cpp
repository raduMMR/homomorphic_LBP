// HE_FR_LBP Homomorphic Encryption for Face Recognition using LBP operator
#include "HE_FR_LBP.h"
#include "basic_primitives.h"
#include "he_lbp.h"

#define NUMBER_OF_PIXELS 1024

void EncRegion::cloneRegion(const EncRegion& copy){
    for(int i=0; i<copy.enc_pixels.size(); i++){
        Ctxt *c = encryptBitVal(vector<long>(nslots, 0));
        *c = *copy.enc_pixels[i];
        this->enc_pixels.push_back(c);
    }

    cout << "Passed\n";

    this->enc_neighbours = vector<vector<Ctxt*>>(copy.enc_neighbours.size());
    for(int i=0; i<copy.enc_neighbours.size(); i++){
        for(int j=0; j<copy.enc_pixels.size(); j++){
            Ctxt *c = encryptBitVal(vector<long>(nslots, 0));
            *c = *copy.enc_neighbours[i][j];
            this->enc_neighbours[i].push_back(c);
        }
    }

}

EncRegion::~EncRegion(){
    for(int i=0; i<enc_pixels.size(); i++){
        delete enc_pixels[i];
    }

    for(int i=0; i<enc_neighbours.size(); i++){
        for(int j=0; j<enc_neighbours[i].size(); j++){
            delete enc_neighbours[i][j];
        }
    }

    cout << "Called EndRegion destructor\n";
}

EncHistogram::EncHistogram(){
    enc_hist = vector<Ctxt*>(256);
    // all lbp_codes from 0 to 255.
}

EncHistogram::~EncHistogram(){
    for(int i=0; i<enc_hist.size(); i++){
            delete enc_hist[i];
    }
}

void EncHistogram::computeNbOfOccurences(vector<Ctxt*> LBP_codes, vector<long> histogram){
    for(int i=0; i<256; i++){
	vector<Ctxt*> lbp_code = encryptIntVal(vector<long>(NUMBER_OF_PIXELS,/*LBP_code=*/ i), /*T_BITS=*/8);

        //enc_hist[i] = compute_z(0, 8, lbp_code, LBP_codes);

	Ctxt *eh = compute_z(0, 8, lbp_code, LBP_codes);

	for(int j=0; j<lbp_code.size(); j++){
		delete lbp_code[j];
	}
	
	vector<long> dec_hist = decryptBitVal(eh);
	delete eh;

	int hist_val = 0;
	for(int j=0; j<NUMBER_OF_PIXELS; j++){
		hist_val += dec_hist[j];
	}
	if(hist_val != histogram[i]){
		cout << "Hist diff " << i << endl;
	} else { cout << "Step " << i << " completed.\n"; }
    }
}

EncRegion* ImageProcessor::encryptRegion(vector<long> pixels, vector<vector<long>> neighbours){
    assert(neighbours.size() == 8);

    EncRegion* er = new EncRegion(this->nslots);
    er->enc_pixels = encryptIntVal(pixels, 8);
    er->enc_neighbours = vector<vector<Ctxt*>>(8);
    for(int i=0; i<8; i++) {
        er->enc_neighbours[i] = encryptIntVal(neighbours[i], 8);
    }

    return er;
}

// compute the LBP codes in the clear for comparison.
vector<long> clearLBPcodes(vector<long> pixels, vector<vector<long>> neighbours){
<<<<<<< HEAD
    vector<long> lbp_codes(1024);
    for(int i=0; i<1024; i++) {
=======
    vector<long> lbp_codes(NUMBER_OF_PIXELS);
    for(int i=0; i<NUMBER_OF_PIXELS; i++) {
>>>>>>> be047bb8c648a38be90f741d91f4759e27dfc54f
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

void ImageProcessor::encryptImage(uint8_t **image){
    assert(this->nslots == NUMBER_OF_PIXELS);
    // for an image of size 256x256 we split it in 16 regions, 4x4
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
            EncRegion *er = encryptRegion(pixels, neighbours);
            this->encrypted_regions.push_back(er);
            cout << "Regiune criptata\n";

            // compare the LBP codes.
            vector<long> plain_lbp = clearLBPcodes(pixels, neighbours);

            cout << "Se calculeaza codurile LBP pentru regiune ...\n";
            vector<Ctxt*> LBP_codes = computeLBP4Region(0);
            vector<long> dec_lbp = decryptIntVal(LBP_codes);
            assert( dec_lbp.size() >= plain_lbp.size() );
            for(int z=0; z<plain_lbp.size(); z++){
                if(plain_lbp[z] != dec_lbp[z]){
                    cout << "Coduri LBP diferite.\n";
			cout << "Index = " << z << ", plain = " << plain_lbp[z] << ", enc = " << dec_lbp[z] << endl;
                    break;
                }
            }
            cout << "Coduri LBP calculate.\n";
		delete er;
            // compare the two histograms.
            cout << "Se compara histogramele ...\n";
            vector<long> plain_hist = computeClearHistogram(plain_lbp);
	cout << "clear histogram done\n";

            EncHistogram eh;
            eh.computeNbOfOccurences(LBP_codes, plain_hist);
            assert(plain_hist.size() == 256);
           /* for(int z=0; z<plain_hist.size(); z++){
                if(plain_hist[z] != decryptBitVal(eh.enc_hist[z])[0] ){
                    cout << "Histograme diferite\n";
                }
            }*/
            cout << "Histograme comparate.\n";

            // cleanup.
            for(int z=0; z<LBP_codes.size(); z++){
                delete LBP_codes[z];
            }

            break;
        }
    }
}

vector<Ctxt*> ImageProcessor::computeLBP4Region(int region_id){
    return hom_LBP(encrypted_regions[region_id]->enc_pixels, encrypted_regions[region_id]->enc_neighbours, 8);
}

EncHistogram ImageProcessor::computeHistogram4Region(int region_id){
    vector<Ctxt*> LBP_codes = computeLBP4Region(region_id);

    EncHistogram eh;

    eh.computeNbOfOccurences(LBP_codes, vector<long>(NUMBER_OF_PIXELS, 0));

    return eh;
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

    ImageProcessor ip;
    if(NSLOTS >= 1024)
        ip.setNumberOfSlots(NUMBER_OF_PIXELS);

    // cout << "Se cripteaza imaginea ...\n";
    ip.encryptImage(image);
    // cout << "Imagine criptata.\n";

    // cout << "Se calculeaza histogramele ...\n";
    // for(int i=0; i<64; i++){
        // EncHistogram enc_hist = ip.computeHistogram4Region(0);
    // }
    // cout << "Histograme calculate.\n";

    // cleanup, delete image matrix.
    for(int i=0; i<258; i++){
        delete image[i];
    }
    delete image;
}
