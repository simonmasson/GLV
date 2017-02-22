#include "entete.h"

void double_and_add(mpz_t p,mpz_t a,mpz_t b,point_t P,mpz_t K,point_t R) {
	mpz_t k;
	mpz_init(k);
	//if K<0, we compute KP = |K|(-P)
	if(mpz_cmp_ui(K,0)<0) {
		mpz_neg(k,K);
		mpz_neg(P->Y,P->Y);
	}
	else {
		mpz_set(k,K);
	}
	//algortihm square and multiply from left to right
	mpz_set(R->X,P->X);
	mpz_set(R->Y,P->Y);
	mpz_set(R->Z,P->Z);
	int n = (int)mpz_sizeinbase(k,2);
	for(int i = n-2 ; i >=0; i--) {
		add_points(p,a,b,R,R,R);
		if(mpz_tstbit(k,i) == 1) {
			add_points(p,a,b,R,P,R);
		}
	}
	//reset P because it was maybe modified
	if(mpz_cmp_ui(K,0)<0) {
		mpz_neg(P->Y,P->Y);
	}
	mpz_clear(k);
}
