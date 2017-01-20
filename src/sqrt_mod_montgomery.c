#include "ecpp.h"
/*
Assumes p odd and a is reduced mod p already.
*/
void sqrt_mod_p_montgomery(mpz_t b, mpz_t a, mpz_t p)
{
    unsigned long p_mod_16 = mpz_fdiv_ui(p, 16);

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
        /* REDCIFY */
        mpm_double(x, a); // 2 a
        /* check 2 a < p */
        mpm_pow(z, x, exp, p); // (2 a)**((p-5)/8)
        mpm_mul(y, z, z); // (2 a)**((p-5)/4)
        mpm_mul(y, y, x); // (2 a)**((p-1)/4)
        mpm_sub_ui(y, y, 1UL); // (2 a)**((p-1)/4) - 1
        mpm_mul(b, a, y); // a * ((2 a)**((p-1)/4) - 1)
        mpm_mul(b , b, z); // (2 a)**((p-5)/8) * a * ((2 a)**((p-1)/4) - 1)
    }
    else if (p_mod_16 == 9UL && bits < 500) // p = ...1001
    {
        /* REDCIFY */
        mpm_double(x, a); // 2 a
        mpz_fdiv_q_2exp(exp, p, 2UL); // (p-1)/4
        mpm_pow(z, x, exp, p); // (2 a)**((p-1)/4)

        mpz_fdiv_q_2exp(exp, p, 4UL); // (p-9)/16
        /* RECIFY 1UL or UNDRECIFY z */
        if (mpz_cmp_ui(z, 1UL) == 0)
        {
            mpz_set_ui(y, 2UL);
            while (mpz_legendre(y, p) != -1) // QNR
                mpz_add_ui(y, y, 1UL);
            /* REDCIFY y */
            mpm_mul(z, y, y); // QNR**2
            mpm_mul(z, z, x); // 2 a * QNR**2
            mpm_pow(x, z, exp, p); // (2 a * QNR**2)**((p-9)/16)
            mpm_mul(b, x, x); // (2 a * QNR**2)**((p-9)/8)
            mpm_mul(z, z, b); // (2 a * QNR**2)**((p-1)/8)
            mpm_sub_ui(z, z, 1UL);  // (2 a * QNR**2)**((p-1)/8) - 1
            mpm_mul(b, x, z); // (2 a * QNR**2)**((p-9)/16) * ((2 a * QNR**2)**((p-1)/8) - 1)
            mpm_mul(b, b, y); // QNR * (2 a * QNR**2)**((p-9)/16) * ((2 a * QNR**2)**((p-1)/8) - 1)
            mpm_mul(b, b, a); // a * QNR * (2 a * QNR**2)**((p-9)/16) * ((2 a * QNR**2)**((p-1)/8) - 1)
        }
        else
        {
            mpm_pow(z, x, exp, p); // (2 a)**((p-9)/16)
            mpm_mul(y, z, z); // (2 a)**((p-9)/8)
            mpm_mul(y, y, x); // (2 a)**((p-1)/8)
            mpm_sub_ui(y, y, 1UL); // (2 a)**((p-1)/8) - 1
            mpm_mul(b, a, z); // a * (2 a)**((p-9)/16)
            mpm_mul(b, b, y); // a * (2 a)**((p-9)/16) * ((2 a)**((p-1)/8) - 1)
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
        /* REDCIFY y */
        mpm_pow(y, y, exp, p); // QNR**q

        mpz_fdiv_q_2exp(exp, exp, 1); // (q-1)/2
        /* REDCIFY a */
        mpm_pow(x, a, exp, p); // a**((q-1)/2)
        mpm_mul(b, a, x); // a**((q+1)/2)
        mpm_mul(x, x, b); // a**q

        /* REDCIFY 1UL */
        while (mpm_cmp_ui(x, 1UL) != 0)
        {
            m = 1;
            mpm_mul(z, x, x);
            while (mpm_cmp_ui(z, 1UL) != 0)
            {
                mpm_mul(z, z, z);
                m += 1;
            }
            mpz_set_ui(exp, 0UL);
            mpz_setbit(exp, r-m-1);
            mpm_pow(z, y, exp, p);
            mpm_mul(b, b, z);
            mpm_mul(x, x, z);
            mpm_mul(x, x, z);
        }
    }

    mpz_clear(exp);
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(z);
    return;
}
