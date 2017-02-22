#include "entete.h"

void point_init(point_t P) {
	//init a point P
	mpz_inits(P->X,P->Y,P->Z,NULL);
}

void point_clear(point_t P) {
	//clear a point P
	mpz_clear(P->X);
	mpz_clear(P->Y);
	mpz_clear(P->Z);
}

void point_printf(point_t R) {
	//print a point R
	gmp_printf("[%Zd:%Zd:%Zd]\n",R->X,R->Y,R->Z);
}

void point_copy(point_t P, point_t R) {
	//copy a point P in a point R
	mpz_set(R->X,P->X);
	mpz_set(R->Y,P->Y);
	mpz_set(R->Z,P->Z);
}

void point_neg(mpz_t p, point_t P, point_t R) {
	//copy (-P) in R (and reduce the yR mod p)
	mpz_set(R->X,P->X);
	mpz_neg(R->Y,P->Y);
	mpz_mod(R->Y,R->Y,p);
	mpz_set(R->Z,P->Z);
}


