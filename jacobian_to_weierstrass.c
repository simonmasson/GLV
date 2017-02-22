#include "entete.h"

void convert_jac(mpz_t p,point_t P) {
	mpz_t tmp;
	mpz_init(tmp);
	
	if(mpz_cmp_ui(P->Z,0) != 0) {
		//compute 1/Z
		mpz_invert(tmp,P->Z,p);
	
		//x = X/Z^2
		mpz_mul(P->X,P->X,tmp);
		mpz_mul(P->X,P->X,tmp);
		mpz_mod(P->X,P->X,p);
	
		//y = Y/Z^3
		mpz_mul(P->Y,P->Y,tmp);
		mpz_mul(P->Y,P->Y,tmp);
		mpz_mul(P->Y,P->Y,tmp);
		mpz_mod(P->Y,P->Y,p);
	
		//z = 1
		mpz_set_ui(P->Z,1);
	}
	mpz_clear(tmp);
}
