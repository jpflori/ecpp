#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "ecpp.h"

int main(int argc, char* argv[])
{
    if (argc != 7)
    {
        printf("Usage: ecpp <bits> <offset> <nb of odd primes> <min factors> <max factors> <gap>\n");
        return 1;
    }

    int B;
    B = atoi(argv[3]);

    int min_factors, max_factors;
    min_factors = atoi(argv[4]);
    max_factors = atoi(argv[5]);

    unsigned long gap = 1UL << atoi(argv[6]);

    clock_t start, ttotal;
    start = clock();
    ttotal = start;

    mpz_t b;
    mpz_init(b);
    mpz_set_ui(b, 1UL);
    mpz_mul_2exp(b, b, atoi(argv[1]));
    mpz_add_ui(b, b, atoi(argv[2]));

    mpz_t p;
    mpz_init(p);
    mpz_nextprime(p, b);
    printf("p: ");
    mpz_out_str(stdout, 10, p);
    printf("\n");
    printf("p: %lf\n", ((double) clock() - start) / CLOCKS_PER_SEC);

    mpz_t c, q;
    mpz_init(c);
    mpz_init(q);

    start = clock();
    mpz_t *primes;
    primes = init_primes(B);
    compute_primes(primes, B);
    printf("Primes: %lf\n", ((double) clock() - start) / CLOCKS_PER_SEC);

    start = clock();
    int *signs;
    signs = init_signs(B);
    compute_signs(signs, primes, B);
    printf("Signs: %lf\n", ((double) clock() - start) / CLOCKS_PER_SEC);

    start = clock();
    int *symbols;
    symbols = init_symbols(B);
    compute_symbols(symbols, signs, primes, B, p);
    printf("Symbols: %lf\n", ((double) clock() - start) / CLOCKS_PER_SEC);

    start = clock();
    int specials, odds, A;
    int *sqprimes;
    A = nb_sqprimes(&specials, &odds, symbols, B);
    printf("Number of special/odd squares: %d/%d (%d)\n", specials, odds, A);
    if (max_factors > 1 + odds)
        max_factors = 1 + odds;
    sqprimes = init_sqprimes(A);
    compute_sqprimes(sqprimes, A, symbols, B);
    printf("Square primes: %lf\n", ((double) clock() - start) / CLOCKS_PER_SEC);

    start = clock();
    mpz_t *sqrtmods;
    sqrtmods = init_sqrtmods(A);
    compute_sqrtmods(sqrtmods, sqprimes, A, signs, primes, p);
    printf("Sqrts: %lf\n", ((double) clock() - start) / CLOCKS_PER_SEC);

    int *indices;
    indices = init_indices(max_factors);
    first_indices(indices, specials, min_factors, max_factors);

    mpz_t d, sqrtmd;
    mpz_init(d);
    mpz_init(sqrtmd);

    mpz_t u, v;
    mpz_init(u);
    mpz_init(v);

    int found;

    found = make_step(c, q, u, v, d, sqrtmd, gap, B, indices, specials, max_factors, A, sqrtmods, sqprimes, signs, primes, p);

    if (found)
    {
        printf("Found solution (B = %d) (proof = %d)\n", B, found);
        printf("p: ");
        mpz_out_str(stdout, 10, p);
        printf("\n");
        printf("c: ");
        mpz_out_str(stdout, 10, c);
        printf("\n");
        printf("q: ");
        mpz_out_str(stdout, 10, q);
        printf("\n");
        printf("d: ");
        mpz_out_str(stdout, 10, d);
        printf("\n");
        printf("sqrt(-d): ");
        mpz_out_str(stdout, 10, sqrtmd);
        printf("\n");
        printf("u: ");
        mpz_out_str(stdout, 10, u);
        printf("\n");
        printf("v: ");
        mpz_out_str(stdout, 10, v);
        printf("\n");
    }
    else
    {
        printf("Not enough discriminants (B = %d)\n", B);
        printf("p: ");
        mpz_out_str(stdout, 10, p);
        printf("\n");
    }

    printf("Total: %lf\n", ((double) clock() - ttotal) / CLOCKS_PER_SEC);

    clear_sqprimes(sqprimes);
    clear_indices(indices);
    clear_sqrtmods(sqrtmods, A);
    clear_symbols(symbols);
    clear_signs(signs);
    clear_primes(primes, B);

    mpz_clear(sqrtmd);
    mpz_clear(d);
    mpz_clear(u);
    mpz_clear(v);
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(c);
    mpz_clear(b);

    return 0;
}
