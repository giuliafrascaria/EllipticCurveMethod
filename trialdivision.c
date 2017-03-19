//
// Created by root on 19/03/17.
//
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

int main(int argc, char ** argv)
{
    if(argc != 2)
    {
        perror("missing number");
        return EXIT_FAILURE;
    }
    //convert input into internal representation
    //declarations and init
    mpz_t n;
    mpz_init(n);
    //conversion
    printf("%s\n", argv[1]);
    if(mpz_set_str(n, (const char *) argv[1], 10) != 0)
    {
        perror("invalid input type");
        return EXIT_FAILURE;
    }

    mpz_t divisor;
    mpz_init(divisor);
    mpz_set_ui(divisor, 2);

    mpz_t limit;
    mpz_init(limit);
    mpz_sqrt(limit, n);

    mpz_t mod;
    mpz_init(mod);

    while(mpz_cmp(divisor, limit) <= 0)
    {
        mpz_mod(mod, n, divisor);
        if(mpz_cmp_ui(mod, 0) == 0)
        {
            gmp_printf("found factor %Zd\n", divisor);
        }

        mpz_nextprime(divisor, divisor);
    }
}
