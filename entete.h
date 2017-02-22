#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include<gmp.h>
#include <time.h>

/* Un point représenté par ses coordonnées projectives */
typedef struct {
  mpz_t X;
  mpz_t Y;
  mpz_t Z;
} point_struct;

/* Astuce à la gmp pour passer des pointeurs qui allouent la mémoire sur le stack. */
typedef point_struct point_t[1];

void point_init(point_t P);
void point_clear(point_t P);
void point_printf(point_t R);
void point_copy(point_t P, point_t R);
void point_neg(mpz_t p,point_t P, point_t R);

void add_points(mpz_t p,mpz_t a,mpz_t b,point_t P, point_t Q, point_t R);
int double_point(mpz_t p,mpz_t a,mpz_t b,point_t P,point_t R);
int multiple_point(mpz_t p,mpz_t a,mpz_t b,point_t P,mpz_t k,point_t R);
void double_and_add(mpz_t p,mpz_t a,mpz_t b,point_t P,mpz_t K,point_t R);

void find_v1_v2(mpz_t n,mpz_t LAMBDA,mpz_t v1x,mpz_t v1y,mpz_t v2x,mpz_t v2y);
void find_v(mpz_t k,mpz_t v1x,mpz_t v1y,mpz_t v2x,mpz_t v2y,mpz_t vx,mpz_t vy);
void find_u(mpz_t k,mpz_t vx,mpz_t vy,mpz_t ux,mpz_t uy);
void decompose(mpz_t K,mpz_t n,mpz_t lambda,mpz_t k1,mpz_t k2);
int test_decompose(mpz_t K,mpz_t n,mpz_t lambda,mpz_t k1,mpz_t k2);

void signedBinary(mpz_t K,int w,mpz_t res);
void multiple_multiplication(mpz_t p,mpz_t a,mpz_t b,int w, mpz_t U, mpz_t V,point_t P, point_t Q, point_t R);

void convert_jac(mpz_t p,point_t P);

void multiple_multiplication_without_precompute(mpz_t p,mpz_t a,mpz_t b,int w, mpz_t U, mpz_t V,point_t P, point_t Q, point_t R,
point_t p1[1<<w][1<<w],
point_t p2[1<<w][1<<w],
point_t p3[1<<w][1<<w],
point_t p4[1<<w][1<<w]);
