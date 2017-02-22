#include "entete.h"

void main(int argc, char** argv) {


	mpz_t p,a,b,n,k;
	point_t P;
	mpz_inits(p,a,b,n,k,NULL);
	point_init(P);
	//p from example 7.
	mpz_set_str(p,"1461501637330902918203684832716283019655932313743",10);
	mpz_set_ui(a,0);
	mpz_set_ui(b,3);
	//n = #(E) (prime)
	mpz_set_str(n,"1461501637330902918203687013445034429194588307251",10);
	//k is a big integer
	mpz_set_str(P->X,"1",10);
	mpz_set_str(P->Y,"2",10);
	mpz_set_str(P->Z,"1",10);
	
////////////////////////////////////////////////////////////////////////////////////
	
	mpz_t pp,lambda,beta,k1,k2;
	mpz_inits(pp,lambda,beta,k1,k2,NULL);
	
	//On vérifie qu'on est dans le cas p = 1 mod 3
	mpz_mod_ui(pp,p,3);
	if(mpz_cmp_ui(pp,1) != 0) {printf("NON!\n");}

	//Variables temporaires utiles
	mpz_t tmp1,tmp2,tmp3;
	mpz_inits(tmp1,tmp2,tmp3,NULL);
	
	//tmp1 = (p+1)/4
	mpz_add_ui(tmp1,p,1);
	mpz_tdiv_q_ui(tmp1,tmp1,4);
	//β est une racine cubique de 1 :  ( (p-3)^((p+1)/4) - 1 ) /2 (voir rapport)
	mpz_sub_ui(tmp3,p,3);
	mpz_powm(tmp2,tmp3,tmp1,p);
	mpz_sub_ui(tmp2,tmp2,1);
	mpz_set_ui(tmp3,2);
	mpz_invert(tmp3,tmp3,p);
	mpz_mul(tmp2,tmp2,tmp3);
	mpz_mod(tmp2,tmp2,p);
	mpz_set(beta,tmp2);
	
	//tmp1 = (n+1)/4
	mpz_add_ui(tmp1,n,1);
	mpz_tdiv_q_ui(tmp1,tmp1,4);
	//λ vérifie λ^2+λ+1 = 0 mod N : ( -(n-3)^((n+1)/4) - 1 ) /2 (voir rapport)
	mpz_sub_ui(tmp3,n,3);
	mpz_powm(tmp2,tmp3,tmp1,n);
	//Il faut choisir la bonne racine !
	mpz_neg(tmp2,tmp2);
	mpz_sub_ui(tmp2,tmp2,1);
	mpz_set_ui(tmp3,2);
	mpz_invert(tmp3,tmp3,n);
	mpz_mul(tmp2,tmp2,tmp3);
	mpz_mod(tmp2,tmp2,n);
	mpz_set(lambda,tmp2);

	//Q = φ(P) avec φ endomorphisme de l'exemple 4.	
	point_t Q;
	point_init(Q);
	mpz_mul(Q->X,beta,P->X);
	mpz_set(Q->Y,P->Y);
	mpz_set_ui(Q->Z,1);
	
	point_t R,R2;
	point_init(R);
	point_init(R2);

	point_t RR,RR2,PP;
	point_init(RR);
	point_init(RR2);
	point_init(PP);

	for(int w = 2; w<8 ; w++) {printf("w=%d\n",w);

		//PRECOMPUTE
		//four tables to store the {iP+jQ} where i and j can be negative. 
		point_t precomp[1<<w][1<<w];
		point_t neg_precomp[1<<w][1<<w];
		point_t negI_precomp[1<<w][1<<w];
		point_t negJ_precomp[1<<w][1<<w];
		mpz_set_ui(RR->X,0);
		mpz_set_ui(RR->Y,1);
		mpz_set_ui(RR->Z,0);

		mpz_set_ui(RR2->X,0);
		mpz_set_ui(RR2->Y,1);
		mpz_set_ui(RR2->Z,0);

		mpz_set_ui(PP->X,0);
		mpz_set_ui(PP->Y,1);
		mpz_set_ui(PP->Z,0);

		//precompute iP+jQ, -iP-jQ,-iP+jQ, and iP-jQ/////////////////////////
		for(int i = 0 ; i < 1<<w ; i ++) {
			for(int j = 0 ; j < 1<<w ; j++) {
				point_init(precomp[i][j]);
				point_init(neg_precomp[i][j]);
				mpz_set(precomp[i][j]->X,RR->X);
				mpz_set(precomp[i][j]->Y,RR->Y);
				mpz_set(precomp[i][j]->Z,RR->Z);
				point_neg(p,precomp[i][j],neg_precomp[i][j]);
				add_points(p,a,b,RR,Q,RR);
		
				point_init(negI_precomp[i][j]);
				point_init(negJ_precomp[i][j]);
				mpz_set(negI_precomp[i][j]->X,RR2->X);
				mpz_set(negI_precomp[i][j]->Y,RR2->Y);
				mpz_set(negI_precomp[i][j]->Z,RR2->Z);
				point_neg(p,negI_precomp[i][j],negJ_precomp[i][j]);
				add_points(p,a,b,RR2,Q,RR2);
			}
		add_points(p,a,b,PP,P,PP);
		point_copy(PP,RR);
		point_neg(p,RR,RR2);
		}
		//end of PRECOMPUTE////////////////////////////////////////////////


		gmp_randstate_t state;
		gmp_randinit_default(state);
	
		clock_t start, end;
		double cpu_time_used;
		double aver_daa=0;
		double aver_glv=0;
	
		FILE *file;
		file = fopen("bench_aver.txt", "a");
		fprintf(file, "%d,		",w);
	
		//gmp_printf("λ = %Zd\nn = %Zd\nβ = %Zd\n\n",lambda,n,beta);
	
		for(int j = 0; j < 5 ; j++){
			mpz_urandomm(k,state,n);
		
			start = clock();
			for(int i = 0 ; i < 1000; i++){
				double_and_add(p,a,b,P,k,R);
			}
			end = clock();
			cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
			aver_daa += cpu_time_used;

			convert_jac(p,R);
			//gmp_printf("R = [%Zd,%Zd,%Zd] (daa)\n",R->X,R->Y,R->Z);
		
			start = clock();
			for(int i = 0 ; i < 1000; i++){
				decompose(k,n,lambda,k1,k2);
				multiple_multiplication_without_precompute(p,a,b,w,k1,k2,P,Q,R,precomp,neg_precomp,negI_precomp,negJ_precomp);
			}
			end = clock();
			cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
			aver_glv += cpu_time_used;
			
			convert_jac(p,R);
			//gmp_printf("R = [%Zd,%Zd,%Zd]\n\n",R->X,R->Y,R->Z);
		}
		fprintf(file,"%f,	%f\n",aver_daa/5.0,aver_glv/5.0);
		
		fclose(file);
	
		//CLEAR PRECOMPUTE/////////////////////////////////
		for(int i = 0 ; i < 1<<w ; i ++) {
			for(int j = 0 ; j < 1<<w ; j++) {
				point_clear(precomp[i][j]);
				point_clear(neg_precomp[i][j]);
				point_clear(negI_precomp[i][j]);
				point_clear(negJ_precomp[i][j]);
			}
		}
		//////////////////////////////////////////////////	
	}

	point_clear(Q); point_clear(R);	point_clear(R2);
	point_clear(RR); point_clear(RR2);point_clear(PP);
	
	
	mpz_clears(tmp1,tmp2,tmp3,NULL);	
	
	mpz_clears(pp,lambda,beta,k1,k2,NULL);
	
////////////////////////////////////////////////////////////////////////////////////
	
	point_clear(P);
	mpz_clears(p,a,b,n,k,NULL);


}
