#include "basic_primitives.h"

// Global variables.
FHEcontext *context;
FHESecKey *secretKey;
EncryptedArray *ea;
bool noPrint = false;
int NSLOTS = 0;

/*************************************************************************************/
void setGlobalVariables(long p, long r, long d, long c, long k, long w, 
               long L, long m, const Vec<long>& gens, const Vec<long>& ords) {
    
    vector<long> gens1, ords1;
    convert(gens1, gens);
    convert(ords1, ords);

    context = new FHEcontext(m, p, r, gens1, ords1);
    buildModChain(*context, L, c);

    ZZX G;
    if (d == 0)
        G = context->alMod.getFactorsOverZZ()[0];
    else
        G = makeIrredPoly(p, d); 

    secretKey = new FHESecKey(*context);
    // const FHEPubKey& publicKey = secretKey;
    secretKey->GenSecKey(w); // A Hamming-weight-w secret key
    addSome1DMatrices(*secretKey); // compute key-switching matrices that we need

    ea = new EncryptedArray(*context, G);

    NSLOTS = ea->size();
}

/*************************************************************************************/
void cleanGlobalVariables() {
    delete context;
    delete secretKey;
    delete ea;
}

/*************************************************************************************/
Ctxt* encryptBitVal (const vector<long> bits) {
    const FHEPubKey& publicKey = *secretKey;
    Ctxt* ctxt = new Ctxt(publicKey);
    ea->encrypt(*ctxt, publicKey, bits);
    return ctxt;
}

/*************************************************************************************/
vector<Ctxt*> encryptIntVal (const vector<long> val, int t_bits) {
    vector<Ctxt*> vec_ctxt(t_bits);
    for(int i=0; i<t_bits; i++){
        vector<long> bits(NSLOTS);
        for(int j=0; j<NSLOTS; j++) {
            bits[j] = (val[j] >> i) & 1;
        }
        vec_ctxt[i] = encryptBitVal(bits);
    }
    return vec_ctxt;
}

/*************************************************************************************/
vector<long> decryptBitVal (const Ctxt *ct) {
    vector<long> decs;
    ea->decrypt(*ct, *secretKey, decs);
    return decs;
}

/*************************************************************************************/
vector<long> decryptIntVal(const vector<Ctxt*> enc_bits) {
    vector<long> vals(NSLOTS, 0);
    for(int i=0; i<enc_bits.size(); i++) {
        vector<long> decrypted_bits = decryptBitVal(enc_bits[i]);
        for(int j=0; j<NSLOTS; j++) {
            vals[j] |= decrypted_bits[j] << i;
        }
    }
    return vals;
}

/*************************************************************************************/
Ctxt* compute_z (int i, int j, vector<Ctxt*>& ct_x, vector<Ctxt*>& ct_y)
{
	Ctxt *ret = NULL;
	if (j == 1)
	{
		ret = encryptBitVal(vector<long>(NSLOTS, 1));
		*ret += *ct_x[i]; 
		*ret += *ct_y[i];
		return ret;
	}
	
	int l;
	l = (j%2 == 0) ? j/2: j/2 + 1; 
	//cout << endl << "compute_z...." << "j="<<j<< "; l=" << l; 
		
	ret = compute_z(i+l, j-l, ct_x, ct_y);
	Ctxt *ct = compute_z (i, l, ct_x, ct_y);	
	*ret *= *ct;
	delete ct;
	
	return ret;
}

/*************************************************************************************/
Ctxt* compute_t (int i, int j, vector<Ctxt*>& ct_x, vector<Ctxt*>& ct_y)
{
	Ctxt *ret = NULL;
	if (j == 1)
	{
		ret  = new Ctxt (*ct_x[i]);
		*ret *= *ct_y[i]; 
		*ret += *ct_x[i];
		return ret;
	}
			
	int l;
	l = (j%2 == 0) ? j/2: j/2 + 1; 

	ret = compute_t(i+l, j-l, ct_x, ct_y);
	Ctxt *ct_z = compute_z (i+l, j-l, ct_x, ct_y);
	Ctxt *ct_t = compute_t (i, l, ct_x, ct_y);
	
	*ct_z *= *ct_t;		
	*ret += *ct_z;	
	
	delete ct_z;
	delete ct_t;	
	return ret;	
}

/*************************************************************************************/
Ctxt* compute_s (int i, int j, vector<Ctxt*>& ct_x, vector<Ctxt*>& ct_y)
{
	Ctxt *ret = NULL;
	if (j == 1)
	{
		Ctxt *ct_1 = encryptBitVal(vector<long>(NSLOTS, 1));
		ret  = new Ctxt (*ct_x[i]);
		*ret *= *ct_y[i]; 
		*ret += *ct_y[i];
		*ret += *ct_1;
		
		delete ct_1;		
		return ret;
	}
			
	int l;
	l = (j%2 == 0) ? j/2: j/2 + 1; 

	ret = compute_t(i+l, j-l, ct_x, ct_y);
	Ctxt *ct_z = compute_z (i+l, j-l, ct_x, ct_y);
	Ctxt *ct_s = compute_s (i, l, ct_x, ct_y);
	
	*ct_z *= *ct_s;	
	*ret += *ct_z;	
	
	delete ct_z;
	delete ct_s;	
	return ret;	
}