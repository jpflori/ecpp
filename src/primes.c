#include <omp.h>

#include "ecpp.h"

/*
  Absolute values of special and odd primes
*/
mpz_t* init_primes(int B)
{
    int i;
    mpz_t *primes;

    primes = (mpz_t *) malloc((B+3)*sizeof(mpz_t));

    for (i = 0; i < B+3; i++)
    {
        mpz_init(primes[i]);
    }

    return primes;
}

void compute_primes(mpz_t *primes, int B)
{
    int i;

    /* special primes */
    /* -4, 8, -8 */
    {
        mpz_set_ui(primes[0], 4UL);
        mpz_set_ui(primes[1], 8UL);
        mpz_set_ui(primes[2], 8UL);
    }
    /* odd primes */
    {
        mpz_set_ui(primes[3], 3UL);
    }
    for (i = 4; i < B+3; i++)
    {
        mpz_nextprime(primes[i], primes[i-1]);
    }
}

void clear_primes(mpz_t *primes, int B)
{
    int i;

    for (i = 0; i < B+3; i++)
    {
        mpz_clear(primes[i]);
    }

    free(primes);
}

/*
  Legendre symbol (-1 / q )
*/
int* init_signs(int B)
{
    int *signs;

    signs = (int *) malloc((B+3)*sizeof(int));

    return signs;
}

void compute_signs(int *signs, mpz_t *primes, int B)
{
    int i;

    /* special primes */
    /* -4, 8, -8 */
    {
        signs[0] = -1;
        signs[1] = 1;
        signs[2] = -1;
    }
    /* odd primes */
    for (i = 3; i < B+3; i++)
    {
        /* signs[i] = ((-mpz_tstbit(primes[i], 1)) << 1) + 1 */
        /* if (mpz_tstbit(primes[i], 1) == 1) */
        if (mpz_tstbit(primes[i], 1))
            signs[i] = -1;
        else
            signs[i] = 1;
    }
}

void clear_signs(int *signs)
{
    free(signs);
}

/*
  Legendre symbol (q* / p) = (p / q)
  Special primes {-4, 8, -8} and odd primes {q*}
 */
int* init_symbols(int B)
{
    int *symbols;

    symbols = (int*) malloc((B+3)*sizeof(int));

    return symbols;
}

void compute_symbols(int *symbols, int *signs, mpz_t *primes, int B, mpz_t p)
{
    int i;
    int sm1p, s2p;

    /* (-1/p) */
    if (mpz_tstbit(p, 1))
        sm1p = -1;
    else
        sm1p = 1;

    /* (2/p) */
    if (mpz_tstbit(p, 1) == mpz_tstbit(p, 2))
        s2p = 1;
    else
        s2p = -1;

    /* special primes */
    /* -4, 8, -8 */
    {
        symbols[0] = sm1p;
        symbols[1] = s2p;
        symbols[2] = sm1p*s2p;
    }

    /* odd primes */
    for (i = 3; i < B+3; i++)
    {
        /* symbols[i] = mpz_legendre(primes[i], p)*(sm1p==-1?signs[i]:1) */
        symbols[i] = mpz_legendre(p, primes[i]);
    }
}

void clear_symbols(int *symbols)
{
    free(symbols);
}

/*
  List of square primes
*/
int nb_sqprimes(int *specials, int *odds, int *symbols, int B)
{
    int i;

    /* special primes */
    (*specials) = 0;
    for (i = 0; i < 3; i++)
        if (symbols[i] == 1)
            (*specials)++;

    /* odd primes */
    (*odds) = 0;
    for (; i < B+3; i++)
        if (symbols[i] == 1)
            (*odds)++;

    return (*odds) + (*specials);
}

int* init_sqprimes(int A)
{
    return (int*) malloc(A*sizeof(int));
}

void compute_sqprimes(int* sqprimes, int A, int *symbols, int B)
{
    int i, j;

    j = 0;
    for (i = 0; i < A; i++)
    {
        while (symbols[j] != 1)
            j++;
        sqprimes[i] = j;
        j++;
    }
}

void clear_sqprimes(int *sqprimes)
{
    free(sqprimes);
}

/*
  Modular sqrt precomputation
  sqrt(q*) mod p
 */
mpz_t* init_sqrtmods(int A)
{
    int i;
    mpz_t *sqrtmods;

    sqrtmods = (mpz_t *) malloc(A*sizeof(mpz_t));

    for (i = 0; i < A; i++)
    {
        mpz_init(sqrtmods[i]);
    }

    return sqrtmods;
}
#include <stdio.h>
void compute_sqrtmods(mpz_t *sqrtmods, int* sqprimes, int A, int *signs, mpz_t *primes, mpz_t p)
{
    int i;
    int j;

    mpz_t tmp;

#pragma omp parallel default(shared) private(tmp, i, j)
    {
    mpz_init(tmp);

    /* odd primes */
#pragma omp for schedule(dynamic,20) nowait
    for (i = 0; i < A; i++)
    {
        j = sqprimes[i];
        if (signs[j] == 1)
            mpz_set(tmp, primes[j]);
        else
            mpz_sub(tmp, p, primes[j]);

        sqrt_mod_p(sqrtmods[i], tmp, p);
    }

    mpz_clear(tmp);
    } /* end of parallel section */
}

void clear_sqrtmods(mpz_t *sqrtmods, int A)
{
    int i;

    for (i = 0; i < A; i++)
    {
        mpz_clear(sqrtmods[i]);
    }

    free(sqrtmods);
}
