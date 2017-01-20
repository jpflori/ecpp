#include "ecpp.h"

/*
Tries to solve 4 p = u^2 + d v^2
Assumes 0  <= d < p and the sqrt of d mod p is precomputed and is of the same parity as d.
*/
int cornacchia2(mpz_t u, mpz_t v, mpz_t d, mpz_t sqrtmd, mpz_t p)
//int cornacchia2_cohen(mpz_t u, mpz_t v, mpz_t d, mpz_t sqrtmd, mpz_t p)
{
    mpz_t sqrtp;

    // We compute p/x as a continued fraction
    // This is the euclidean algorithm
    // Should use a half-gcd computation
    // And low precision estimate before going to multiprecision

    // 2 sqrt(p)
    mpz_init(sqrtp);
    mpz_sqrt(sqrtp, p);
    mpz_mul_2exp(sqrtp, sqrtp, 1);

    // sqrt(-d)
    mpz_set(u, sqrtmd);

    // 2 p
    mpz_mul_2exp(v, p, 1);

    // Euclidean algorithm
    while (mpz_cmp(u, sqrtp) > 0)
    {
        mpz_mod(v, v, u);
        mpz_swap(u, v);
    }

    // Divisibility condition
    mpz_mul(sqrtp, u, u);
    mpz_mul_2exp(v, p, 2);
    mpz_sub(v, v, sqrtp);
    mpz_clear(sqrtp);
    if (mpz_divisible_p(v, d))
    {
        // Square condition
        mpz_divexact(v, v, d);
        if (mpz_perfect_square_p(v))
        {
            mpz_sqrt(v, v);
            return 1;
        }
    }

    // No solution
    return 0;
}

int cornacchia2_morain(mpz_t u, mpz_t v, mpz_t d, mpz_t sqrtmd, mpz_t p)
//int cornacchia2(mpz_t u, mpz_t v, mpz_t d, mpz_t sqrtmd, mpz_t p)
{
    mpz_t sqrtp, uu, vv, q, m;
    int solution;

    // We compute p/x as a continued fraction
    // This is the euclidean algorithm
    // Should use a half-gcd computation
    // And low precision estimate before going to multiprecision

    // 2 sqrt(p)
    mpz_init(sqrtp);
    mpz_sqrt(sqrtp, p);
    mpz_mul_2exp(sqrtp, sqrtp, 1);

    // sqrt(-d)
    mpz_set(u, sqrtmd);

    // 2 p
    mpz_init(uu);
    mpz_mul_2exp(uu, p, 1);

    // reduced
    mpz_set_ui(v, 1UL);
    mpz_init(vv);
    mpz_set_ui(vv, 0UL);

    // quotient
    mpz_init(q);

    // m
    mpz_init(m);

    // Euclidean algorithm
    while (mpz_cmp(u, sqrtp) > 0)
    {
        mpz_fdiv_qr(q, uu, uu, u);
        mpz_swap(u, uu);
        mpz_mul(m, q, v);
        mpz_add(vv, vv, m);
        mpz_swap(v, vv);
    }

    // check equality
    mpz_mul(uu, u, u);
    mpz_mul(vv, v, v);
    mpz_mul(vv, vv, d);
    mpz_add(m, uu, vv);
    mpz_mul_2exp(q, p, 2);

    solution = mpz_cmp(m, q);

    mpz_clear(sqrtp);
    mpz_clear(uu);
    mpz_clear(vv);
    mpz_clear(q);
    mpz_clear(m);

    return solution==0?1:0;
}
