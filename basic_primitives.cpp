#include "basic_primitives.h"

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <NTL/ZZ.h>
#include <NTL/fileio.h>
#include <NTL/BasicThreadPool.h>
NTL_CLIENT
#include "./../EvalMap.h"
#include "./../powerful.h"

// Global variables.
FHEcontext *context;
FHESecKey *secretKey;
EncryptedArray *ea;
bool noPrint = false;
int NSLOTS = 0;
int T_BITS = 8;
ZZX G;

/*************************************************************************************/
double clock_diff(const clock_t &t1, const clock_t &t2){
    return double(t2 - t1) / CLOCKS_PER_SEC;
}

/*************************************************************************************/
void setGlobalVariables(long p, long r, long d, long c, long w, 
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
	cout << "NSLOTS = " << NSLOTS << endl;
}

/*************************************************************************************/
void writeContextToFile(const char *filename, int d, int p){
	fstream keyFile(filename, fstream::out|fstream::trunc);

	writeContextBase(keyFile, *context);
    keyFile << *context << endl;
	keyFile << *secretKey << endl;
	keyFile << d << endl;
	keyFile << p << endl;

	keyFile.close();
	keyFile.flush();
}

/*************************************************************************************/
void readContextFromFile(const char *filename){
	if(context != NULL) delete context;
	if(secretKey != NULL) delete secretKey;
	if( ea != NULL) delete ea;

	fstream keyFile(filename, fstream::in);

    unsigned long m1, p1, r1;
    vector<long> gens, ords;
    readContextBase(keyFile, m1, p1, r1, gens, ords);

	context = new FHEcontext(m1, p1, r1, gens, ords);
	keyFile >> *context;
	secretKey = new FHESecKey(*context);
	keyFile >> *secretKey;

	int d, p;
	keyFile >> d;
	keyFile >> p;
	ZZX G;
    if (d == 0)
        G = context->alMod.getFactorsOverZZ()[0];
    else
        G = makeIrredPoly(p, d); 
	ea = new EncryptedArray(*context, G);

	cout << "EA SLOTS = " << ea->size() << endl << endl;

	keyFile.close();
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
bool needsBootstrapping(const Ctxt* ctxt){
    double curNoise = log(ctxt->getNoiseVar())/2;
    double noiseThreshold = log(ctxt->modSwitchAddedNoiseVar())*0.55;

    if (curNoise>noiseThreshold && ctxt->log_of_ratio()>-0.5)
    {
        // cerr << "Ctxt::findBaseSet warning: already at lowest level\n";
        return true;
    }
    return false;
}

/*************************************************************************************/
Ctxt* compute_z (int i, int j, const vector<Ctxt*>& ct_x, const vector<Ctxt*>& ct_y)
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
Ctxt* compute_t (int i, int j, const vector<Ctxt*>& ct_x, const vector<Ctxt*>& ct_y)
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
Ctxt* compute_s (int i, int j, const vector<Ctxt*>& ct_x, const vector<Ctxt*>& ct_y)
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

/*************************************************************************************/
static long mValues[][14] = { 
//{ p, phi(m),  m,    d, m1,  m2, m3,   g1,    g2,    g3,ord1,ord2,ord3, c_m}
  {  2,    48,   105, 12,  3,  35,  0,    71,    76,    0,  2,  2,   0, 100},
  {  2,   600,  1023, 10, 11,  93,  0,   838,   584,    0, 10,  6,   0, 100}, // m=(3)*11*{31} m/phim(m)=1.7    C=24  D=2 E=1
  {  2,  1200,  1705, 20, 11, 155,  0,   156,   936,    0, 10,  6,   0, 100}, // m=(5)*11*{31} m/phim(m)=1.42   C=34  D=2 E=2
  {  2,  1728,  4095, 12,  7,  5, 117,  2341,  3277, 3641,  6,  4,   6, 100}, // m=(3^2)*5*7*{13} m/phim(m)=2.36 C=26 D=3 E=2
  {  2,  2304,  4641, 24,  7,  3, 221,  3979,  3095, 3760,  6,  2,  -8, 100}, // m=3*7*(13)*{17} :-( m/phim(m)=2.01 C=45 D=4 E=3
  {  2,  4096,  4369, 16, 17, 257,  0,   258,  4115,    0, 16,-16,   0, 100}, // m=17*(257) :-( m/phim(m)=1.06 C=61 D=3 E=4
  {  2, 12800, 17425, 40, 41, 425,  0,  5951,  8078,    0, 40, -8,   0, 100}, // m=(5^2)*{17}*41 m/phim(m)=1.36 C=93  D=3 E=3
  {  2, 15004, 15709, 22, 23, 683,  0,  4099, 13663,    0, 22, 31,   0, 100}, // m=23*(683) m/phim(m)=1.04      C=73  D=2 E=1
  {  2, 16384, 21845, 16, 17,   5,257,  8996, 17477, 21591, 16, 4, -16,1600}, // m=5*17*(257) :-( m/phim(m)=1.33 C=65 D=4 E=4
  {  2, 18000, 18631, 25, 31, 601,  0, 15627,  1334,    0, 30, 24,   0,  50}, // m=31*(601) m/phim(m)=1.03      C=77  D=2 E=0
  {  2, 18816, 24295, 28, 43, 565,  0, 16386, 16427,    0, 42, 16,   0, 100}, // m=(5)*43*{113} m/phim(m)=1.29  C=84  D=2 E=2
  {  2, 21168, 27305, 28, 43, 635,  0, 10796, 26059,    0, 42, 18,   0, 100}, // m=(5)*43*{127} m/phim(m)=1.28  C=86  D=2 E=2
  {  2, 23040, 28679, 24, 17,  7, 241, 15184,  4098,28204, 16,  6, -10,1000}, // m=7*17*(241) m/phim(m)=1.24    C=63  D=4 E=3
  {  2, 24000, 31775, 20, 41, 775,  0,  6976, 24806,    0, 40, 30,   0, 100}, // m=(5^2)*{31}*41 m/phim(m)=1.32 C=88  D=2 E=2
  {  2, 26400, 27311, 55, 31, 881,  0, 21145,  1830,    0, 30, 16,   0, 100}, // m=31*(881) m/phim(m)=1.03      C=99  D=2 E=0
  {  2, 27000, 32767, 15, 31,  7, 151, 11628, 28087,25824, 30,  6, -10, 150},
  {  2, 31104, 35113, 36, 37, 949,  0, 16134,  8548,    0, 36, 24,   0, 400}, // m=(13)*37*{73} m/phim(m)=1.12  C=94  D=2 E=2
  {  2, 34848, 45655, 44, 23, 1985, 0, 33746, 27831,    0, 22, 36,   0, 100}, // m=(5)*23*{397} m/phim(m)=1.31  C=100 D=2 E=2
  {  2, 42336, 42799, 21, 127, 337, 0, 25276, 40133,    0,126, 16,   0,  20}, // m=127*(337) m/phim(m)=1.01     C=161 D=2 E=0
  {  2, 45360, 46063, 45, 73, 631,  0, 35337, 20222,    0, 72, 14,   0, 100}, // m=73*(631) m/phim(m)=1.01      C=129 D=2 E=0
  {  2, 46080, 53261, 24, 17, 13, 241, 43863, 28680,15913, 16, 12, -10, 100}, // m=13*17*(241) m/phim(m)=1.15   C=69  D=4 E=3
  {  2, 49500, 49981, 30, 151, 331, 0,  6952, 28540,    0,150, 11,   0, 100}, // m=151*(331) m/phim(m)=1        C=189 D=2 E=1
  {  2, 54000, 55831, 25, 31, 1801, 0, 19812, 50593,    0, 30, 72,   0, 100}, // m=31*(1801) m/phim(m)=1.03     C=125 D=2 E=0
  {  2, 60016, 60787, 22, 89, 683,  0,  2050, 58741,    0, 88, 31,   0, 200}, // m=89*(683) m/phim(m)=1.01      C=139 D=2 E=1

  {  7,    36,    57,  3,  3,  19,  0,    20,    40,    0,  2, -6,   0, 100}, // m=3*(19) :-( m/phim(m)=1.58 C=14 D=3 E=0

  { 17,    48,   105, 12,  3,  35,  0,    71,    76,    0,  2,  2,   0, 100}, // m=3*(5)*{7} m/phim(m)=2.18 C=14 D=2 E=2
  { 17,   576,  1365, 12,  7,   3, 65,   976,   911,  463,  6,  2,   4, 100}, // m=3*(5)*7*{13} m/phim(m)=2.36  C=22  D=3
  { 17, 18000, 21917, 30, 101, 217, 0,  5860,  5455,    0, 100, 6,   0, 100}, // m=(7)*{31}*101 m/phim(m)=1.21  C=134 D=2 
  { 17, 30000, 34441, 30, 101, 341, 0,  2729, 31715,    0, 100, 10,  0, 100}, // m=(11)*{31}*101 m/phim(m)=1.14 C=138 D=2
  { 17, 40000, 45551, 40, 101, 451, 0, 19394,  7677,    0, 100, 10,  0,2000}, // m=(11)*{41}*101 m/phim(m)=1.13 C=148 D=2
  { 17, 46656, 52429, 36, 109, 481, 0, 46658,  5778,    0, 108, 12,  0, 100}, // m=(13)*{37}*109 m/phim(m)=1.12 C=154 D=2
  { 17, 54208, 59363, 44, 23, 2581, 0, 25811,  5199,    0, 22, 56,   0, 100}, // m=23*(29)*{89} m/phim(m)=1.09  C=120 D=2
  { 17, 70000, 78881, 10, 101, 781, 0, 67167, 58581,    0, 100, 70,  0, 100}, // m=(11)*{71}*101 m/phim(m)=1.12 C=178 D=2

  {127,   576,  1365, 12,  7,   3, 65,   976,   911,  463,  6,  2,   4, 100}, // m=3*(5)*7*{13} m/phim(m)=2.36   C=22  D=3
  {127,  1200,  1925, 20,  11, 175, 0,  1751,   199,    0, 10, 6,    0, 100}, //  m=(5^2)*{7}*11 m/phim(m)=1.6   C=34 D=2
  {127,  2160,  2821, 30,  13, 217, 0,   652,   222,    0, 12, 6,    0, 100}, // m=(7)*13*{31} m/phim(m)=1.3     C=46 D=2
  {127, 18816, 24295, 28, 43, 565,  0, 16386, 16427,    0, 42, 16,   0, 100}, // m=(5)*43*{113} m/phim(m)=1.29   C=84  D=2
  {127, 26112, 30277, 24, 17, 1781, 0, 14249, 10694,    0, 16, 68,   0, 100}, // m=(13)*17*{137} m/phim(m)=1.15  C=106 D=2
  {127, 31752, 32551, 14, 43,  757, 0,  7571, 28768,    0, 42, 54,   0, 100}, // m=43*(757) :-( m/phim(m)=1.02   C=161 D=3
  {127, 46656, 51319, 36, 37, 1387, 0, 48546, 24976,    0, 36, -36,  0, 200}, //m=(19)*37*{73}:-( m/phim(m)=1.09 C=141 D=3
  {127, 49392, 61103, 28, 43, 1421, 0,  1422, 14234,    0, 42, 42,   0,4000}, // m=(7^2)*{29}*43 m/phim(m)=1.23  C=110 D=2
  {127, 54400, 61787, 40, 41, 1507, 0, 30141, 46782,    0, 40, 34,   0, 100}, // m=(11)*41*{137} m/phim(m)=1.13  C=112 D=2
  {127, 72000, 77531, 30, 61, 1271, 0,  7627, 34344,    0, 60, 40,   0, 100}  // m=(31)*{41}*61 m/phim(m)=1.07   C=128 D=2
};
#define num_mValues (sizeof(mValues)/(14*sizeof(long)))

#define OUTER_REP (3)
#define INNER_REP (3)

static bool dry = false;

/*************************************************************************************/
void setContextWithBootstrapping(long idx, int p, int r, int L, int c, long B, long skHwt, bool cons, int build_cache)
{
	Vec<long> mvec;
	vector<long> gens;
	vector<long> ords;

	long phim = mValues[idx][1];
	long m = mValues[idx][2];
	assert(GCD(p, m) == 1);

	append(mvec, mValues[idx][4]);
	if (mValues[idx][5]>1) append(mvec, mValues[idx][5]);
	if (mValues[idx][6]>1) append(mvec, mValues[idx][6]);
	gens.push_back(mValues[idx][7]);
	if (mValues[idx][8]>1) gens.push_back(mValues[idx][8]);
	if (mValues[idx][9]>1) gens.push_back(mValues[idx][9]);
	ords.push_back(mValues[idx][10]);
	if (abs(mValues[idx][11])>1) ords.push_back(mValues[idx][11]);
	if (abs(mValues[idx][12])>1) ords.push_back(mValues[idx][12]);

	if (!noPrint) {
		cout << "*** TestIt";
		if (isDryRun()) cout << " (dry run)";
		cout << ": p=" << p
		<< ", r=" << r
		<< ", L=" << L
		<< ", B=" << B
		<< ", c=" << c
		<< ", m=" << m
		<< " (=" << mvec << "), gens="<<gens<<", ords="<<ords
		<< endl;
		cout << "Computing key-independent tables..." << std::flush;
	}
	setTimersOn();
	setDryRun(false); // Need to get a "real context" to test bootstrapping

	double t = -GetTime();
	context = new FHEcontext(m, p, r, gens, ords);
	context->bitsPerLevel = B;
	buildModChain(*context, L, c,/*extraBits=*/7);

	// FIXME: The extraBits is an exceedingly ugly patch, used to bypass the
	//   issue that buildModChain must be called BEFORE the context is made
	//   bootstrappable (else the "powerful" basis is not initialized correctly.)
	//   This is a bug, the value 7 is sometimes the right one, but seriously??

	context->makeBootstrappable(mvec, /*t=*/0, cons, build_cache);
	t += GetTime();

	if (skHwt>0) context->rcData.skHwt = skHwt;
	if (!noPrint) {
		cout << " done in "<<t<<" seconds\n";
		cout << "  e="    << context->rcData.e
		<< ", e'="   << context->rcData.ePrime
		<< ", alpha="<< context->rcData.alpha
		<< ", t="    << context->rcData.skHwt
		<< "\n  ";
		context->zMStar.printout();
	}
	setDryRun(dry); // Now we can set the dry-run flag if desired

	long nPrimes = context->numPrimes();
	IndexSet allPrimes(0,nPrimes-1);
	double bitsize = context->logOfProduct(allPrimes)/log(2.0);
	if (!noPrint)
		cout << "  "<<nPrimes<<" primes in chain, total bitsize="
		<< ceil(bitsize) << ", secparam="
		<< (7.2*phim/bitsize -110) << endl;

	long p2r = context->alMod.getPPowR();
	context->zMStar.set_cM(mValues[idx][13]/100.0);

	for (long numkey=0; numkey<OUTER_REP; numkey++) { // test with 3 keys

	t = -GetTime();
	if (!noPrint) cout << "Generating keys, " << std::flush;
	secretKey = new FHESecKey(*context);
	FHEPubKey& publicKey = *secretKey;
	secretKey->GenSecKey(64);      // A Hamming-weight-64 secret key
	addSome1DMatrices(*secretKey); // compute key-switching matrices that we need
	addFrbMatrices(*secretKey);
	if (!noPrint) cout << "computing key-dependent tables..." << std::flush;
	secretKey->genRecryptData();
	t += GetTime();
	if (!noPrint) cout << " done in "<<t<<" seconds\n";

	zz_p::init(p2r);
	zz_pX poly_p = random_zz_pX(context->zMStar.getPhiM());
	PowerfulConversion pConv(context->rcData.p2dConv->getIndexTranslation());
	HyperCube<zz_p> powerful(pConv.getShortSig());
	pConv.polyToPowerful(powerful, poly_p);
	ZZX ptxt_poly = conv<ZZX>(poly_p);
	PolyRed(ptxt_poly, p2r, true); // reduce to the symmetric interval

	ZZX poly2;
	Ctxt c1(publicKey);

	secretKey->Encrypt(c1,ptxt_poly,p2r);
	for (long num=0; num<INNER_REP; num++) { 
		publicKey.reCrypt(c1);
		secretKey->Decrypt(poly2,c1);

		if (ptxt_poly == poly2) cout << "  *** reCryption succeeds!!\n";
		else cout << "Eroare la recriptare\n";
	}
	}
}
