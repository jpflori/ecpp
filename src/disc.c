#include "ecpp.h"

/*
  Computes the (opposite of the) fundamental discriminant corresponding to
  - a set of odd primes,
  - a parity.
  So this is positive and already multiplied by 4 if needed.
*/
void fund_disc_from_indices(mpz_t d, int even, int *indices, int max_factors, mpz_t *primes)
{
    int i;

    /* parity */
    if (even)
        mpz_set_ui(d, 2UL);
    else
        mpz_set_ui(d, 1UL);

    /* odd primes */
    for (i = 0; i < max_factors && indices[i] != -1; i++)
    {
        mpz_mul(d, d, primes[indices[i]]);
    }

    /*  multiply by 4 if d % 4 = 1 or 2 */
    /* if (mpz_tstbit(d, 0) != 1 || mpz_tstbit(d, 1) != 1) */
    if (mpz_fdiv_ui(d, 4UL) != 3UL)
        mpz_mul_2exp(d, d, 2);
}

/*
  Check that a set of special/odd primes yields a negative discriminant.
  If so computes the (opposite of the) fundamental discriminant and
  its sqrt mod p.
*/
int fund_disc_and_sqrt_from_indices(mpz_t d, mpz_t sqrtmd, int* indices, int max_factors, mpz_t *sqrtmods, int *sqprimes, int *signs, mpz_t *primes, mpz_t p)
{
    int i;
    int j;
    int sign;

    /* compute sign */
    sign = 1;
    for (i = 0; i < max_factors && indices[i] != -1; i++)
    {
        j = sqprimes[indices[i]];
        if (signs[j] == -1)
            sign = -sign;
    }
    /* check sign */
    if (sign == 1)
        return 0;

    /* compute discriminant and its square root */
    {
        mpz_set_ui(d, 1UL);
        mpz_set_ui(sqrtmd, 1UL);
    }
    for (i = 0; i < max_factors && indices[i] != -1; i++)
    {
        j = sqprimes[indices[i]];
        mpz_mul(d, d, primes[j]);
        mpz_mul(sqrtmd, sqrtmd, sqrtmods[indices[i]]);
        mpz_mod(sqrtmd, sqrtmd, p);
    }

    return 1;
}
