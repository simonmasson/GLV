#include "entete.h"

int multiple_point(mpz_t p,mpz_t a,mpz_t b,point_t P,mpz_t k,point_t R) {
	//Naive algorithm
	
	mpz_set(R->X,P->X);
	mpz_set(R->Y,P->Y);
	mpz_set(R->Z,P->Z);
	
	mpz_t i,kk,L;
	mpz_inits(i,kk,L,NULL);
	//if k<0, compute |k|P and return -R at the end
	if(mpz_cmp_ui(k,0)<0) {
		mpz_neg(kk,k);
	}
	else {
		mpz_set(kk,k);
	}
	mpz_sub_ui(L,k,1);
	for(mpz_set_ui(i,0);mpz_cmp(i,L)<0;mpz_add_ui(i,i,1)) {
		add_points(p,a,b,P,R,R);
	}
	if(mpz_cmp_ui(k,0)<0) {
		mpz_neg(R->Y,R->Y);
	}
	mpz_clears(i,kk,L,NULL);
}
