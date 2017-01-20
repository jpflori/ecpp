#include "ecpp.h"

void factor_trial(mpz_t c, mpz_t q, mpz_t n, mpz_t* primes, int B)
{
    int i;
    unsigned long bit;

    mpz_set_ui(c, 1UL);
    mpz_set(q, n);

    /* even part */
    bit = mpz_scan1(q, 0);
    mpz_mul_2exp(c, c, bit);
    mpz_fdiv_q_2exp(q, q, bit);

    /* odd part */
    for (i = 0; i < B; i++)
    {
        while(mpz_divisible_p(q, primes[i]))
        {
            mpz_mul(c, c, primes[i]);
            mpz_divexact(q, q, primes[i]);
        }
    }
}
