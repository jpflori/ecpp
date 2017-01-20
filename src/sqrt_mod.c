#include "ecpp.h"

/*
Takes an sqrt mod p and make it a sqrt mod 4p
*/
int sqrt_mod_4p(mpz_t u, mpz_t d, mpz_t p)
{
    // Check we got the right root mod 4 p
    if (mpz_tstbit(u, 0) != mpz_tstbit(d, 0))
    {
        mpz_sub(u, p, u);
        return 1;
    }

    return 0;
}

/*
Assumes a is reduced mod p already.
*/
void sqrt_mod_p(mpz_t b, mpz_t a, mpz_t p)
{
    unsigned long p_mod_16 = mpz_fdiv_ui(p, 16UL);

    if (!(p_mod_16 & 1UL))  // p == ...xxx0 => p == 2
    {
        mpz_set(b, a);
        return;
    }

    mpz_t exp, x, y, z;
    mpz_init(exp);
    mpz_init(x);
    mpz_init(y);
    mpz_init(z);

    long bits = mpz_sizeinbase(p, 2);

    if (p_mod_16 & 2UL) // p == ...xx11
    {
        mpz_add_ui(exp, p, 1UL); // p+1
        mpz_fdiv_q_2exp(exp, exp, 2UL); // (p+1) / 4
        mpz_powm(b, a, exp, p); // a**((p+1)/4)
    }
    else if (p_mod_16 & 4UL) // p == ...x101
    {
        mpz_fdiv_q_2exp(exp, p, 3UL); // (p-5)/8
        mpz_mul_2exp(x, a, 1UL); // 2 a
        mpz_powm(z, x, exp, p); // (2 a)**((p-5)/8)
        mpz_mul(y, z, z);
        mpz_mod(y, y, p); // (2 a)**((p-5)/4)
        mpz_mul(y, y, x);
        mpz_mod(y, y, p); // (2 a)**((p-1)/4)
        mpz_sub_ui(y, y, 1UL); // (2 a)**((p-1)/4) - 1
        mpz_mul(b, a, y);
        mpz_mod(b, b, p); // a * ((2 a)**((p-1)/4) - 1)
        mpz_mul(b , b, z);
        mpz_mod(b, b, p); // (2 a)**((p-5)/8) * a * ((2 a)**((p-1)/4) - 1)
    }
    else if (p_mod_16 == 9UL && bits < 500) // p = ...1001
    {
        mpz_mul_2exp(x, a, 1UL); // 2 a
        mpz_fdiv_q_2exp(exp, p, 2UL); // (p-1)/4
        mpz_powm(z, x, exp, p); // (2 a)**((p-1)/4)
        mpz_fdiv_q_2exp(exp, p, 4UL); // (p-9)/16

        if (mpz_cmp_ui(z, 1UL) == 0)
        {
            mpz_set_ui(y, 2UL);
            while (mpz_legendre(y, p) != -1) // QNR
                mpz_add_ui(y, y, 1UL);
            mpz_mul(z, y, y);
            mpz_mod(z, z, p); // QNR**2
            mpz_mul(z, z, x);
            mpz_mod(z, z, p); // 2 a * QNR**2
            mpz_powm(x, z, exp, p); // (2 a * QNR**2)**((p-9)/16)
            mpz_mul(exp, x, x);
            mpz_mod(exp, exp, p); // (2 a * QNR**2)**((p-9)/8)
            mpz_mul(z, z, exp);
            mpz_mod(z, z, p); // (2 a * QNR**2)**((p-1)/8)
            mpz_sub_ui(z, z, 1UL);  // (2 a * QNR**2)**((p-1)/8) - 1
            mpz_mul(b, a, y);
            mpz_mod(b, b, p); // a * QNR
            mpz_mul(b, b, x);
            mpz_mod(b, b, p); // a * QNR * (2 a * QNR**2)**((p-9)/16)
            mpz_mul(b, b, z);
            mpz_mod(b, b, p); // a * QNR * (2 a * QNR**2)**((p-9)/16) * ((2 a * QNR**2)**((p-1)/8) - 1)
        }
        else
        {
            mpz_powm(z, x, exp, p); // (2 a)**((p-9)/16)
            mpz_mul(y, z, z);
            mpz_mod(y, y, p); // (2 a)**((p-9)/8)
            mpz_mul(y, y, x);
            mpz_mod(y, y, p); // (2 a)**((p-1)/8)
            mpz_sub_ui(y, y, 1UL); // (2 a)**((p-1)/8) - 1
            mpz_mul(b, a, z);
            mpz_mod(b, b, p); // a * (2 a)**((p-9)/16)
            mpz_mul(b, b, y);
            mpz_mod(b, b, p); // a * (2 a)**((p-9)/16) * ((2 a)**((p-1)/8) - 1)
        }
    }
    else
    {
        long m, r;

        mpz_sub_ui(exp, p, 1UL);
        r = mpz_scan1(exp, 0); // p - 1 = ...10...0
        mpz_fdiv_q_2exp(exp, exp, r); // q = ...1

        mpz_set_ui(y, 2UL);
        while (mpz_legendre(y, p) != -1) // QNR
            mpz_add_ui(y, y, 1UL);
        mpz_powm(y, y, exp, p); // QNR**q

        mpz_fdiv_q_2exp(exp, exp, 1); // (q-1)/2
        mpz_powm(x, a, exp, p); // a**((q-1)/2)
        mpz_mul(b, a, x);
        mpz_mod(b, b, p); // a**((q+1)/2)
        mpz_mul(x, x, b);
        mpz_mod(x, x, p); // a**q

        while (mpz_cmp_ui(x, 1UL) != 0)
        {
            m = 1;
            mpz_mul(z, x, x);
            mpz_mod(z, z, p);
            while (mpz_cmp_ui(z, 1UL) != 0)
            {
                mpz_mul(z, z, z);
                mpz_mod(z, z, p);
                m += 1;
            }
            mpz_set_ui(exp, 0UL);
            mpz_setbit(exp, r-m-1);
            mpz_powm(z, y, exp, p);
            mpz_mul(b, b, z);
            mpz_mod(b, b, p);
            mpz_mul(x, x, z);
            mpz_mod(x, x, p);
            mpz_mul(x, x, z);
            mpz_mod(x, x, p);
        }
    }

    mpz_clear(exp);
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(z);
    return;
}
