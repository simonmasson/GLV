#include "entete.h"

void multiple_multiplication(mpz_t p,mpz_t a,mpz_t b,int w, mpz_t U, mpz_t V,point_t P, point_t Q, point_t R) {
	//log(w) < 16 for the precomputing
	
	point_t RR,RR2,PP,Pcopy,QQ,Qcopy;
	mpz_t u,v,ubloc, vbloc,ufin,vfin, carre,tmp;
	int iu, jv;
	
	point_init(PP); point_init(QQ); point_init(RR);point_init(RR2);point_init(Pcopy);point_init(Qcopy);
	
	point_copy(P,Pcopy);
	point_copy(Q,Qcopy);
	
	mpz_inits(u,v,ubloc,vbloc,ufin,vfin,carre,tmp,NULL);
	
	//We can suppose that u and v are positive (changing P, Q !)
	if(mpz_cmp_ui(U,0) < 0) {
		mpz_neg(u,U);
		mpz_neg(P->Y,P->Y);
	}
	else {
		mpz_set(u,U);
	}
	
	if(mpz_cmp_ui(V,0) < 0) {
		mpz_neg(v,V);
		mpz_neg(Q->Y,Q->Y);
	}
	else {
		mpz_set(v,V);
	}
	
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
	//precompute iP+jQ, -iP-jQ,-iP+jQ, and iP-jQ
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
	//end of PRECOMPUTE
	
	//u and v are transformed with their signed binary representation and then in 2n-bit form corresponding.
	signedBinary(u,w,u);
	signedBinary(v,w,v);
	
	//init R = O
	mpz_set_ui(RR->X,0);
	mpz_set_ui(RR->Y,1);
	mpz_set_ui(RR->Z,0);
	
	//t is the length of max(u,v)
	mpz_cmp(u,v)>0 ? mpz_set(tmp,u) : mpz_set(tmp,v);
	int t = (int)(mpz_sizeinbase(tmp,2));
	//d is the number of step to do
	int d = t/(2*w);
	if(t%(2*w) != 0) {
		d ++;
	}
	mpz_ui_pow_ui(carre,2,w);
	for(int i=d-1 ; i>=0 ; i--) {
		//R = 2^w * R
		double_and_add(p,a,b,RR,carre,RR);
		mpz_ui_pow_ui(tmp,4,w*i);//4 because the blocs have twice the number of bits
		mpz_tdiv_q(ubloc,u,tmp);
		mpz_tdiv_q(vbloc,v,tmp);
		for(int j = 0 ; j < 2*w ; j++) {
			mpz_clrbit(u,2*w*i + j);
			mpz_clrbit(v,2*w*i + j);
		}
		mpz_set_ui(ufin,0);
		mpz_set_ui(vfin,0);
		//ufin and vfin take the real values of ubloc and vbloc : 00, 10 and 11 are converted in 0,1,-1.
		for(int s = 0 ; s < w ; s++ ) {
			if(mpz_tstbit(ubloc,2*s) == 0 && mpz_tstbit(ubloc,2*s+1) == 1) {
				mpz_add_ui(ufin,ufin,1<<s);
			}
			if(mpz_tstbit(ubloc,2*s) == 1 && mpz_tstbit(ubloc,2*s+1) == 1) {
				mpz_sub_ui(ufin,ufin,1<<s);
			}
			if(mpz_tstbit(vbloc,2*s) == 0 && mpz_tstbit(vbloc,2*s+1) == 1) {
				mpz_add_ui(vfin,vfin,1<<s);
			}
			if(mpz_tstbit(vbloc,2*s) == 1 && mpz_tstbit(vbloc,2*s+1) == 1) {
				mpz_sub_ui(vfin,vfin,1<<s);
			}
		}
		iu = mpz_get_ui(ufin);
		jv = mpz_get_ui(vfin);
		
		//use the precompute table in function of the sign of ufin and vfin.
		if(mpz_cmp_ui(ufin,0) >= 0 && mpz_cmp_ui(vfin,0) >= 0) {
			add_points(p,a,b,RR,precomp[iu][jv],RR);
		}
		if(mpz_cmp_ui(ufin,0) < 0 && mpz_cmp_ui(vfin,0) >= 0) {
			add_points(p,a,b,RR,negI_precomp[iu][jv],RR);
		}
		if(mpz_cmp_ui(ufin,0) >= 0 && mpz_cmp_ui(vfin,0) < 0) {
			add_points(p,a,b,RR,negJ_precomp[iu][jv],RR);
		}
		if(mpz_cmp_ui(ufin,0) < 0 && mpz_cmp_ui(vfin,0) < 0) {
			add_points(p,a,b,RR,neg_precomp[iu][jv],RR);
		}
	}
	point_copy(RR,R);
	
	point_clear(PP);
	point_clear(QQ);
	point_clear(RR);
	point_clear(RR2);
	
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
	
	point_copy(Pcopy,P);
	point_copy(Qcopy,Q);
	
	point_clear(Pcopy);
	point_clear(Qcopy);

	mpz_clears(u,v,ubloc,vbloc,ufin,vfin,carre,tmp,NULL);
}


//////////////////////third with precompute before !


void multiple_multiplication_without_precompute(mpz_t p,mpz_t a,mpz_t b,int w, mpz_t U, mpz_t V,point_t P, point_t Q, point_t R,
point_t p1[1<<w][1<<w],
point_t p2[1<<w][1<<w],
point_t p3[1<<w][1<<w],
point_t p4[1<<w][1<<w]) {
	point_t RR;
	mpz_t u,v,ubloc, vbloc,ufin,vfin, carre,tmp;
	int iu, jv;
	
	point_init(RR);
	
	mpz_inits(u,v,ubloc,vbloc,ufin,vfin,carre,tmp,NULL);
	
	//We can suppose that u and v are positive (changing P, Q when needed !)
	if(mpz_cmp_ui(U,0) < 0) {
		mpz_neg(u,U);
	}
	else {
		mpz_set(u,U);
	}
	
	if(mpz_cmp_ui(V,0) < 0) {
		mpz_neg(v,V);
	}
	else {
		mpz_set(v,V);
	}
	
	//u and v are transformed with their signed binary representation and then in 2n-bit form corresponding.
	signedBinary(u,w,u);
	signedBinary(v,w,v);
	
	//init R = O
	mpz_set_ui(RR->X,0);
	mpz_set_ui(RR->Y,1);
	mpz_set_ui(RR->Z,0);
	
	//t is the length of max(u,v)
	mpz_cmp(u,v)>0 ? mpz_set(tmp,u) : mpz_set(tmp,v);
	int t = (int)(mpz_sizeinbase(tmp,2));
	//d is the number of step to do
	int d = t/(2*w);
	if(t%(2*w) != 0) {
		d ++;
	}
	mpz_ui_pow_ui(carre,2,w);
	for(int i=d-1 ; i>=0 ; i--) {
		//R = 2^w * R
		double_and_add(p,a,b,RR,carre,RR);
		mpz_ui_pow_ui(tmp,4,w*i);//4 because the blocs have twice the number of bits
		mpz_tdiv_q(ubloc,u,tmp);
		mpz_tdiv_q(vbloc,v,tmp);
		for(int j = 0 ; j < 2*w ; j++) {
			mpz_clrbit(u,2*w*i + j);
			mpz_clrbit(v,2*w*i + j);
		}
		mpz_set_ui(ufin,0);
		mpz_set_ui(vfin,0);
		//ufin and vfin take the real values of ubloc and vbloc : 00, 10 and 11 are converted in 0,1,-1.
		for(int s = 0 ; s < w ; s++ ) {
			if(mpz_tstbit(ubloc,2*s) == 0 && mpz_tstbit(ubloc,2*s+1) == 1) {
				mpz_add_ui(ufin,ufin,1<<s);
			}
			if(mpz_tstbit(ubloc,2*s) == 1 && mpz_tstbit(ubloc,2*s+1) == 1) {
				mpz_sub_ui(ufin,ufin,1<<s);
			}
			if(mpz_tstbit(vbloc,2*s) == 0 && mpz_tstbit(vbloc,2*s+1) == 1) {
				mpz_add_ui(vfin,vfin,1<<s);
			}
			if(mpz_tstbit(vbloc,2*s) == 1 && mpz_tstbit(vbloc,2*s+1) == 1) {
				mpz_sub_ui(vfin,vfin,1<<s);
			}
		}
		iu = mpz_get_ui(ufin);
		jv = mpz_get_ui(vfin);
		
		//use the precompute table in function of the sign of ufin and vfin.
		if(mpz_cmp_ui(U,0) < 0) {
			mpz_neg(ufin,ufin);
		}
		if(mpz_cmp_ui(V,0) < 0) {
			mpz_neg(vfin,vfin);
		}
		if(mpz_cmp_ui(ufin,0) >= 0 && mpz_cmp_ui(vfin,0) >= 0) {add_points(p,a,b,RR,p1[iu][jv],RR);}
		if(mpz_cmp_ui(ufin,0) < 0 && mpz_cmp_ui(vfin,0) >= 0) {add_points(p,a,b,RR,p3[iu][jv],RR);}
		if(mpz_cmp_ui(ufin,0) >= 0 && mpz_cmp_ui(vfin,0) < 0) {add_points(p,a,b,RR,p4[iu][jv],RR);}
		if(mpz_cmp_ui(ufin,0) < 0 && mpz_cmp_ui(vfin,0) < 0) {add_points(p,a,b,RR,p2[iu][jv],RR);}
	}
	point_copy(RR,R);
	
	point_clear(RR);

	mpz_clears(u,v,ubloc,vbloc,ufin,vfin,carre,tmp,NULL);
}
