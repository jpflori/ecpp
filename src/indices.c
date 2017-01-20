#include "ecpp.h"

int* init_indices(int max_factors)
{
    return (int*) malloc(max_factors*sizeof(int));
}

void first_indices(int *indices, int specials, int min_factors, int max_factors)
{
    int i;

    /* non-last slot cannot include special primes */
    for (i = 0; i < min_factors - 1; i++)
        indices[i] = min_factors - i + specials - 1;
    /* last slot includes special primes */
    if (min_factors != 0)
        indices[i] = 0;
    /* unused slots */
    for (i = min_factors; i < max_factors; i++)
        indices[i] = -1;
}

int next_indices(int *indices, int specials, int max_factors, int A)
{
    int i, j;

    /* does not make sense with special primes */
    /*
    if (max_factors == 0)
        return 0;
    */

    /* find last used index and increment it */
    i = 0;
    while ((i < max_factors-1) && (indices[i+1] != -1))
        i++;
    indices[i] += 1;

    /* time for a new factor? */
    if (indices[i] == A - i)
    {
        i++;

        /* reached max number of factors? */
        if (i == max_factors)
            return 0;

        /* slot including special primes */
        {
            indices[i] = 0;
            i--;
        }
        /* slots without special primes */
        j = specials;
        for (; i >= 0; i--)
            indices[i] = j++;
    }
    else
    {
        j = i;

        /* increment as needed */
        while (i > 0 && indices[i] == indices[i-1])
        {
            i--;
            indices[i] += 1;
        }

        /* slots without special primes */
        while (i < j - 1)
        {
            indices[i+1] = j - i - 1 + specials - 1;
            i++;
        }
        /* slot with special primes */
        if (i == j-1)
            indices[j] = 0;
    }

    return 1;
}

void clear_indices(int *indices)
{
    free(indices);
}
