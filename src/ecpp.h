#ifndef ECPP_H
#define ECPP_H

#include <stdlib.h>

#include "gmp.h"

/* primes.c */
mpz_t* init_primes(int B);
void compute_primes(mpz_t *primes, int B);
void clear_primes(mpz_t *primes, int B);

int* init_signs(int B);
void compute_signs(int *signs, mpz_t *primes, int B);
void clear_signs(int *signs);

int* init_symbols(int B);
void compute_symbols(int *symbols, int *signs, mpz_t *primes, int B, mpz_t p);
void clear_symbols(int *symbols);

int nb_sqprimes(int *specials, int *odds, int *symbols, int B);
int* init_sqprimes(int A);
void compute_sqprimes(int *sqprimes, int A, int *symbols, int B);
void clear_sqprimes(int *sqprimes);

mpz_t* init_sqrtmods(int A);
void compute_sqrtmods(mpz_t *sqrtmods, int* sqprimes, int A, int *signs, mpz_t *primes, mpz_t p);
void clear_sqrtmods(mpz_t *sqrtmods, int A);

/* disc.c */
void fund_disc_from_indices(mpz_t d, int even, int *indices, int max_factors, mpz_t *primes);
int fund_disc_and_sqrt_from_indices(mpz_t d, mpz_t sqrtmd, int* indices, int max_factors, mpz_t *sqrtmods, int *sqprimes, int *signs, mpz_t *primes, mpz_t p);

/* indices.c */
int* init_indices(int max_factors);
void first_indices(int *indices, int specials, int min_factors, int max_factors);
int next_indices(int *indices, int specials, int max_factors, int A);
void clear_indices(int *indices);

/* sqrt_mod.c */
int sqrt_mod_4p(mpz_t u, mpz_t d, mpz_t p);
void sqrt_mod_p(mpz_t b, mpz_t a, mpz_t p);

/* cornacchia2.c */
int cornacchia2(mpz_t u, mpz_t v, mpz_t d, mpz_t sqrtmd, mpz_t p);

/* factor.c */
void factor_trial(mpz_t c, mpz_t q, mpz_t n, mpz_t* primes, int B);

/* step.c */
void init_step(int *indices, int *specials, int *odds, int *A, mpz_t *sqrtmods, int *sqprimes, int *symbols, int *signs, mpz_t *primes, int B, int min_factors, int max_factors, mpz_t p);

int make_step(mpz_t c, mpz_t q, mpz_t u, mpz_t v, mpz_t d, mpz_t sqrtmd, unsigned long gap, int sbound, int *indices, int specials, int max_factors, int A, mpz_t *sqrtmods, int *sqprimes, int *signs, mpz_t *primes, mpz_t p);

#endif /* ECPP_H */
