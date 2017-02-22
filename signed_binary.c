#include "entete.h"

void signedBinary(mpz_t K,int w,mpz_t RES) {
	mpz_t k,res;
	mpz_inits(k,res,NULL);
	mpz_cmp_ui(K,0) < 0 ? mpz_neg(k,K) : mpz_set(k,K);
	mpz_set_ui(res,0);
	//i browses the bits of K
	int i,T;
	i = 0;
	//n is the length of K
	int n = (int)mpz_sizeinbase(K,2);
	//kk is the array of the signed decomposition of K
	int kk[n+1];
	//the last value is a particular case...
	kk[n] = 1;
	while(i < n) {
		//if it is the end, copy k
		begining :
		if(i+w-1>=n-w) {
			while(i<n) {
				kk[i] = mpz_tstbit(k,i);
				i++;
			}
			//when we doesn't go there, it is a particular case
			kk[n] = 0;
			break;
		}
		//if the bit is a 0, copy it.
		if(mpz_tstbit(k,i) == 0) {
			kk[i] = mpz_tstbit(k,i);
			i++;
			goto begining;
		}
		//if the bit is a 1, compute the window
		if(mpz_tstbit(k,i) == 1) {
			//compute fen
			int fen = 0;
			int puiss = 1;
			for(int l = 0; l < w; l++) {
				fen += mpz_tstbit(k,i+l)*puiss;
				puiss *= 2;
			}
			//if the first bit after the window is 0
			if(mpz_tstbit(k,i+w) == 0) {
				for(int l = i ; l <= i+w; l++) {
					kk[l] = mpz_tstbit(k,l);
				}
				i += w+1;
				goto begining;
			}
			//if the first bit after the window is 1
			if(mpz_tstbit(k,i+w) == 1) {
				//T is the rank of the first 0 following
				T = mpz_scan0(k,i+w+1);
				//we put 0 from i+w to T
				for(int l = i+w; l <T; l++) {
					kk[l] = 0;
				}
				//we set on T-th bit the 1 in 0 !!IN k!! (not in kk!)
				mpz_setbit(k,T);
				//Compensate with fen :
				fen = (1<<w) - fen;
				for(int m = 0; m <w; m++) {
					kk[i+m] = -((fen&(1<<m))>>m);//the opposite of the m-th bit of fen
				}
				i = T;
			}
		}
	}
	
	for(int j = 0 ; j < n+1 ; j++) {
		if(kk[j] == 0) {
			mpz_clrbit(res,2*j);
			mpz_clrbit(res,2*j+1);
		}
		if(kk[j] == 1) {
			mpz_clrbit(res,2*j);
			mpz_setbit(res,2*j+1);
		}
		if(kk[j] == -1) {
			mpz_setbit(res,2*j);
			mpz_setbit(res,2*j+1);
		}
	}
	
	if(mpz_cmp_ui(K,0) < 0) {
		mpz_neg(res,res);
	}
	mpz_set(RES,res);
	mpz_clears(k,res,NULL);
}
