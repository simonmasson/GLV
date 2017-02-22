#include "entete.h"

void add_points(mpz_t p,mpz_t a,mpz_t b,point_t P, point_t Q, point_t R) {
	//a and b such that y^2 = x^3 + ax + b, P and Q are points of the curve
	//fill R with coordinates of P+Q, with jacobian formulas.
	
	//We use an other point RR for problems like if we compute add(P,P,P) for example.
	point_t RR;
	point_init(RR);

	//P=O so P+Q = O+Q = Q
	if(mpz_cmp_ui(P->X,0) == 0 && mpz_cmp_ui(P->Y,1) == 0 && mpz_cmp_ui(P->Z,0) == 0) {
		point_copy(Q,RR);
	}
	//P=O so P+Q = P+O = P
	else if(mpz_cmp_ui(Q->X,0) == 0 && mpz_cmp_ui(Q->Y,1) == 0 && mpz_cmp_ui(Q->Z,0) == 0) {
		point_copy(P,RR);
	}

	//if P and Q have the same x
	else if(mpz_cmp(P->X,Q->X) == 0) {
		if(mpz_cmp(P->Y,Q->Y)==0) {
			//we double a point
			if(mpz_cmp_ui(P->Y,0) == 0) {
				//2-torsion point
				mpz_set_ui(RR->X,0);
				mpz_set_ui(RR->Y,1);
				mpz_set_ui(RR->Z,0);
			}
			else {
				//double formula
				mpz_t S,M,T,tmp;
				mpz_inits(S,M,T,tmp,NULL);
				
				//S = 4X1Y1^2
				mpz_set_ui(S,4);
				mpz_mul(S,S,P->X);
				mpz_mul(S,S,P->Y);
				mpz_mul(S,S,P->Y);
				mpz_mod(S,S,p);
				
				//M = 3X1^2 + aZ1^4
				mpz_set_ui(M,3);
				mpz_mul(M,M,P->X);
				mpz_mul(M,M,P->X);
				mpz_powm_ui(tmp,P->Z,4,p);
				mpz_addmul(M,a,tmp);
				mpz_mod(M,M,p);
				
				//RR->X = -2S+M^2
				mpz_mul_ui(RR->X,S,2);
				mpz_neg(RR->X,RR->X);
				mpz_mul(tmp,M,M);
				mpz_add(RR->X,RR->X,tmp);
				mpz_mod(RR->X,RR->X,p);
				
				
				//RR->Y = −8Y1^4 + M(S−(RR->X))
				mpz_powm_ui(tmp,P->Y,4,p);
				mpz_mul_ui(RR->Y,tmp,8);
				mpz_neg(RR->Y,RR->Y);
				mpz_sub(tmp,S,RR->X);
				mpz_mul(tmp,tmp,M);
				mpz_add(RR->Y,tmp,RR->Y);
				mpz_mod(RR->Y,RR->Y,p);
				
				//RR->Z = 2Y1Z1
				mpz_set_ui(RR->Z,2);
				mpz_mul(RR->Z,RR->Z,P->Y);
				mpz_mul(RR->Z,RR->Z,P->Z);
				mpz_mod(RR->Z,RR->Z,p);
				
				mpz_clears(S,M,T,tmp,NULL);
				
			}
		}
		else {
			//P = -Q so P+Q = O
			mpz_set_ui(RR->X,0);
			mpz_set_ui(RR->Y,1);
			mpz_set_ui(RR->Z,0);
		}
	}
	else {
		//general case
		mpz_t U1,U2,S1,S2,H,r,tmp;
		mpz_inits(U1,U2,S1,S2,H,r,tmp,NULL);
		
		//U1 = X1Z2^2
		mpz_set(U1,P->X);
		mpz_mul(U1,U1,Q->Z);
		mpz_mul(U1,U1,Q->Z);
		mpz_mod(U1,U1,p);
		
		//U2 = X2Z1^2
		mpz_set(U2,Q->X);
		mpz_mul(U2,U2,P->Z);
		mpz_mul(U2,U2,P->Z);
		mpz_mod(U2,U2,p);
		
		//S1 = Y1Z2^3
		mpz_set(S1,P->Y);
		mpz_mul(S1,S1,Q->Z);
		mpz_mul(S1,S1,Q->Z);
		mpz_mul(S1,S1,Q->Z);
		mpz_mod(S1,S1,p);
		
		//S2 = Y2Z1^3
		mpz_set(S2,Q->Y);
		mpz_mul(S2,S2,P->Z);
		mpz_mul(S2,S2,P->Z);
		mpz_mul(S2,S2,P->Z);
		mpz_mod(S2,S2,p);
		
		//H = U2-U1
		mpz_sub(H,U2,U1);
		mpz_mod(H,H,p);
		
		//r = S2-S1
		mpz_sub(r,S2,S1);
		mpz_mod(r,r,p);
		
		//X3 = -H^3-2U1H^2+r^2
		mpz_powm_ui(RR->X,H,3,p);
		mpz_neg(RR->X,RR->X);
		mpz_mul_ui(tmp,U1,2);
		mpz_mul(tmp,tmp,H);
		mpz_mul(tmp,tmp,H);
		mpz_sub(RR->X,RR->X,tmp);
		mpz_mul(tmp,r,r);
		mpz_add(RR->X,RR->X,tmp);
		mpz_mod(RR->X,RR->X,p);
				
		//Y3 = -S1H^3+r(U1H^2 - X3)
		mpz_mul(tmp,U1,H);
		mpz_mul(tmp,tmp,H);
		mpz_sub(tmp,tmp,RR->X);
		mpz_mul(RR->Y,tmp,r);
		mpz_mul(tmp,S1,H);
		mpz_mul(tmp,tmp,H);
		mpz_mul(tmp,tmp,H);
		mpz_sub(RR->Y,RR->Y,tmp);
		mpz_mod(RR->Y,RR->Y,p);
		
		//Z3 = Z1Z2H
		mpz_mul(RR->Z,P->Z,Q->Z);
		mpz_mul(RR->Z,RR->Z,H);
		mpz_mod(RR->Z,RR->Z,p);
		
		
		mpz_clears(U1,U2,S1,S2,H,r,tmp,NULL);
	}
	point_copy(RR,R);
	point_clear(RR);
}

