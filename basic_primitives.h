#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>
#include "./../FHE.h"
#include "./../timing.h"
#include "./../EncryptedArray.h"
#include <NTL/lzz_pXFactoring.h>
#include <cassert>
#include <assert.h>
#include <cstdio>

// Global variables.
extern FHEcontext *context;
extern FHESecKey *secretKey;
extern EncryptedArray *ea;
extern bool noPrint;
extern int NSLOTS;
extern int T_BITS;

double clock_diff(const clock_t &t1, const clock_t &t2);

void setGlobalVariables(long p, long r, long d, long c, long w, 
               long L, long m, const Vec<long>& gens, const Vec<long>& ords);

void setContextWithBootstrapping(long idx, int p, int r, int L, int c, long B, long skHwt, bool cons, int build_cache);

void cleanGlobalVariables();

Ctxt* encryptBitVal (const vector<long> bit);
vector<Ctxt*> encryptIntVal (const vector<long> val, int t_bits);

vector<long> decryptBitVal (const Ctxt *ct);
vector<long> decryptIntVal(const vector<Ctxt*> enc_bits);


Ctxt* compute_z (int i, int j, vector<Ctxt*>& ct_x, vector<Ctxt*>& ct_y);
Ctxt* compute_t (int i, int j, vector<Ctxt*>& ct_x, vector<Ctxt*>& ct_y);
Ctxt* compute_s (int i, int j, vector<Ctxt*>& ct_x, vector<Ctxt*>& ct_y);