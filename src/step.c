#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "ecpp.h"

void init_step(int *indices, int *specials, int* odds, int *A, mpz_t *sqrtmods, int *sqprimes, int *symbols, int *signs, mpz_t *primes, int B, int min_factors, int max_factors, mpz_t p)
{
    clock_t start;

    start = clock();
    compute_symbols(symbols, signs, primes, B, p);
    printf("Symbols: %lf\n", ((double) clock() - start) / CLOCKS_PER_SEC);

    start = clock();
    *A = nb_sqprimes(specials, odds, symbols, B);
    printf("Number of special/odd squares: %d/%d (%d)\n", *specials, *odds, *A);
    if (max_factors > 1 + *odds)
        abort();
    compute_sqprimes(sqprimes, *A, symbols, B);
    printf("Square primes: %lf\n", ((double) clock() - start) / CLOCKS_PER_SEC);

    start = clock();
    compute_sqrtmods(sqrtmods, sqprimes, *A, signs, primes, p);
    printf("Sqrts: %lf\n", ((double) clock() - start) / CLOCKS_PER_SEC);

    first_indices(indices, *specials, min_factors, max_factors);
}

int make_step(mpz_t c, mpz_t q, mpz_t u, mpz_t v, mpz_t d, mpz_t sqrtmd, unsigned long gap, int sbound, int *indices, int specials, int max_factors, int A, mpz_t *sqrtmods, int *sqprimes, int *signs, mpz_t *primes, mpz_t p)
{
    int found;
    found = 0;

    clock_t start, tcornacchia, tcheck, tsqrt, tprimality;
    tcornacchia = tcheck = tsqrt = tprimality = 0;

    long i, j, k, l, m;
    i = 0L;
    j = (1L << 10) - 1L;
    k = 0L;
    l = 0L;
    m = 0L;

    do
    {
        i += 1L;
        if (((i & j) ^ j) == 0L)
        {
            k += 1L;
            printf("\r%ld kloops %ld discs %ld norms", k, l, m);
            fflush(stdout);
        }

        /* Check and compute disc */
        start = clock();
        found = fund_disc_and_sqrt_from_indices(d, sqrtmd, indices, max_factors, sqrtmods, sqprimes, signs, primes, p);
        tcheck += clock() - start;
        if (!found)
            continue;
        l += 1L;

        start = clock();
        sqrt_mod_4p(sqrtmd, d, p);
        tsqrt += clock() - start;

        start = clock();
        found = cornacchia2(u, v, d, sqrtmd, p);
        tcornacchia += clock() - start;

        // Check for next step
        if (found)
        {
            m += 1L;
            start = clock();
            found = 0;
            mpz_add_ui(q, p, 1UL);
            mpz_sub(q, q, u);
            factor_trial(c, q, q, primes, sbound);
            if (mpz_cmp_ui(c, gap) > 0)
            {
                found = mpz_probab_prime_p(q, 1);
            }
            /*
            if (mpz_even_p(q))
            {
                mpz_fdiv_q_2exp(q, q, 1UL);
                found = mpz_probab_prime_p(q, 1);
            }
            */
            tprimality += clock() - start;

            if (found)
                break;

            start = clock();
            mpz_add_ui(q, p, 1UL);
            mpz_add(q, q, u);
            factor_trial(c, q, q, primes, sbound);
            if (mpz_cmp_ui(c, gap) > 0)
            {
                found = mpz_probab_prime_p(q, 1);
            }
            /*
            if mpz_even_p(q))
            {
                mpz_fdiv_q_2exp(q, q, 1UL);
                found = mpz_probab_prime_p(q, 1);
            }
            */
            tprimality += clock() - start;
        }
    } while(!found && next_indices(indices, specials, max_factors, A));
    if (k > 0)
        printf("\n");

    printf("Number of loops: %ld\n", i);
    printf("Number of discriminants: %ld\n", l);
    printf("Number of norms: %ld\n", m);

    printf("Check: %lf\n", ((double) tcheck) / CLOCKS_PER_SEC);
    printf("Sqrt mod p to mod 4 p: %lf\n", ((double) tsqrt) / CLOCKS_PER_SEC);
    printf("Cornacchia: %lf\n", ((double) tcornacchia) / CLOCKS_PER_SEC);
    printf("Primality: %lf\n", ((double) tprimality) / CLOCKS_PER_SEC);

    return found;
}
