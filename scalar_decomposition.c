#include "entete.h"

void find_v1_v2(mpz_t n,mpz_t LAMBDA,mpz_t v1x,mpz_t v1y,mpz_t v2x,mpz_t v2y) {
	//compute extended euclidin algorithm in n and lambda. Stop at the last remainder grower than sqrt(n)
	//fill the values of v1 and v2
	mpz_t lambda,rem,rrem,u,v,uu,vv,tmp1,tmp2,tmp3,q;
	mpz_inits(lambda,rem,rrem,u,v,uu,vv,tmp1,tmp2,tmp3,q,NULL);
	//rem = 1*rem + 0*rrem
	mpz_set(rem,n);
	mpz_set_ui(u,1);
	mpz_set_ui(v,0);
	
	mpz_mod(lambda,LAMBDA,n);
	
	//rrem = 0*rem + 1*rrem
	mpz_set(rrem,lambda);
	mpz_set_ui(uu,0);
	mpz_set_ui(vv,1);
	
	//tmp2 = sqrt(n)
	mpz_sqrt(tmp2,n);
	
	while(mpz_cmp(rrem,tmp2) > 0) {
		//q = rem/rrem
		mpz_tdiv_q(q,rem,rrem);
		
		//[rrem,rem] = [rem - q*rem,rrem]
		mpz_set(tmp1,rrem);
		mpz_set(rrem,rem);		
		mpz_submul(rrem,q,tmp1);
		mpz_set(rem,tmp1);
		
		//[uu,u] = [u-q*uu,uu]
		mpz_set(tmp1,uu);
		mpz_set(uu,u);
		mpz_submul(uu,q,tmp1);
		mpz_set(u,tmp1);
		
		//[vv,v] = [v-q*vv,vv]
		mpz_set(tmp1,vv);
		mpz_set(vv,v);
		mpz_submul(vv,q,tmp1);
		mpz_set(v,tmp1);
	}
	//v1
	mpz_set(v1x,rrem);
	mpz_neg(v1y,vv);	
	
	//v2 = min ( [rem,-v] , [rem - q*rrem,v-q*vv] )    (with q = rem/rrem)
	
	//Choose the minimal length (their square)
	//tmp1 = rem^2 + u^2
	//tmp2 = (rem -qrrem)^2 +(v-qvv)^2
	
	//tmp1 = rem^2 + u^2
	mpz_pow_ui(tmp1,rem,2);
	mpz_pow_ui(tmp3,u,2);
	mpz_add(tmp1,tmp1,tmp3);
	
	//q = rem/rrem
	mpz_tdiv_q(q,rem,rrem);

	//tmp3 = rem - q * rrem
	mpz_set(tmp3,rem);
	mpz_submul(tmp3,q,rrem);
	
	//(trick) v2y = -(v - q*vv)
	mpz_set(tmp2,v);
	mpz_submul(tmp2,q,vv);
	//trick for less calculus
	mpz_neg(v2y,tmp2);
	
	//tmp2 = (v-q*vv)^2 + (rem - q*rrem)^2
	mpz_pow_ui(tmp2,tmp2,2);
	mpz_addmul(tmp2,tmp3,tmp3);
	
	//v2 is the shortest
	if(mpz_cmp(tmp1,tmp2) < 0) {
		mpz_set(v2x, rem);
		mpz_neg(v2y,v);
	}
	else {
		mpz_set(v2x,tmp3);
		//v2y is already defined		
	}
	mpz_clears(rem,rrem,u,v,uu,vv,tmp1,tmp2,tmp3,q,NULL);
}

void find_v(mpz_t k,mpz_t v1x,mpz_t v1y,mpz_t v2x,mpz_t v2y,mpz_t vx,mpz_t vy) {
	//fill b1 and b2 by the nearest integers of β1 and β2 such that (k,0) = β1*v1+β2*v2	
	//Suppose that (v1,v2) is a free family.
	mpq_t beta1, beta2,ddet;
	mpq_inits(beta1,beta2,ddet,NULL);
	//for fractions calculus and rounding
	mpz_t n,d,b1,b2,tmp,det;
	mpz_inits(n,d,b1,b2,tmp,det,NULL);
	
	//det = v1x*v2y - v1y*v2x
	mpz_mul(det,v1x,v2y);
	mpz_submul(det,v1y,v2x);
	//ddet = det
	mpq_set_z(ddet,det);
	
	//Mat inverse = 1/det [[v2y,-v2x],[-v1y,v1x]]
	//β1 = v2y*k/det
	mpz_mul(tmp,v2y,k);
	mpq_set_z(beta1,tmp);
	mpq_div(beta1,beta1,ddet);
	//β2 = -v1y*k/det
	mpz_mul(tmp,v1y,k);
	mpq_set_z(beta2,tmp);
	mpq_div(beta2,beta2,ddet);
	mpq_neg(beta2,beta2);
		
	//β1 := n/d, b1 the nearest integer from β1
	mpz_set(n,mpq_numref(beta1));
	mpz_set(d,mpq_denref(beta1));
	mpz_mul_ui(n,n,2);
	mpz_add(n,n,d);
	mpz_mul_ui(d,d,2);
	mpz_fdiv_q(b1,n,d);
	//beta2 := n/d, b2 the nearest integer from β2
	mpz_set(n,mpq_numref(beta2));
	mpz_set(d,mpq_denref(beta2));
	mpz_mul_ui(n,n,2);
	mpz_add(n,n,d);
	mpz_mul_ui(d,d,2);
	mpz_fdiv_q(b2,n,d);
	
	//v = b1*v1+b2*v2
	mpz_mul(vx,b1,v1x);
	mpz_addmul(vx,b2,v2x);
	mpz_mul(vy,b1,v1y);
	mpz_addmul(vy,b2,v2y);
	
	mpq_clears(beta1,beta2,ddet,NULL);
	mpz_clears(n,d,b1,b2,det,NULL);
}

void find_u(mpz_t k,mpz_t vx,mpz_t vy,mpz_t ux,mpz_t uy) {
	//u = (k,0) - v
	mpz_sub(ux,k,vx);
	mpz_neg(uy,vy);
}

void decompose(mpz_t K,mpz_t n,mpz_t lambda,mpz_t k1,mpz_t k2) {
	mpz_t k,v1x,v1y,v2x,v2y,vx,vy,tmp;
	mpz_inits(k,v1x,v1y,v2x,v2y,vx,vy,tmp,NULL);
	find_v1_v2(n,lambda,v1x,v1y,v2x,v2y);
	mpz_mod(k,K,n);
	find_v(k,v1x,v1y,v2x,v2y,vx,vy);
	find_u(k,vx,vy,k1,k2);
	
	mpz_set(tmp,k1);
	mpz_addmul(tmp,lambda,k2);
	mpz_mod(tmp,tmp,n);
	if (mpz_cmp(tmp,k) != 0) {
		printf("Error during the decomposition of k !\n");
	}
	mpz_clears(k,v1x,v1y,v2x,v2y,vx,vy,tmp,NULL);
}

int test_decompose(mpz_t K,mpz_t n,mpz_t lambda,mpz_t k1,mpz_t k2) {
	mpz_t tmp,k;
	mpz_inits(tmp,k,NULL);
	mpz_set(tmp,k1);
	mpz_mod(k,K,n);
	mpz_addmul(tmp,lambda,k2);
	mpz_mod(tmp,tmp,n);
	if (mpz_cmp(tmp,k) == 0) {
		mpz_clears(tmp,k,NULL);
		return 1;
	}
	mpz_clears(tmp,k,NULL);
	return 0;
}
