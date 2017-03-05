//
// Created by root on 05/03/17.
//


#include <stdlib.h>
#include <stdio.h>
#include "dataStructures.h"


void randomEC(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd);


int main(int argc, char * argv)
{
    if(argc != 2)
    {
        return EXIT_FAILURE;
    }
    //convert input into internal representation
    struct problemData pd;
    if(mpz_set_str(pd.n, &argv[1], 10) != 1)
    {
        perror("invalid input type");
        return EXIT_FAILURE;
    }
    //choose criteria

    //choose random EC
    struct ECpoint startP;
    struct ellipticCurve EC;
    randomEC(EC, startP, pd);
    //stage one

    //stage two

    return EXIT_SUCCESS;
}

void randomEC(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd)
{
    //random sigma, derive u v anc C from this

    //initialize random state
    gmp_randstate_t state;
    gmp_randinit_mt(state);

    //generate random sigma in [6, n-1]
    unsigned long int six = 6;
    mpz_t maxsigma;
    mpz_sub_ui(maxsigma, pd.n, six);
    mpz_urandomm(EC.sigma, state, maxsigma);
    mpz_add_ui(EC.sigma, EC.sigma, six);

    //the curve is determined by these values

    //save coordinates for initial point Q (MOntgomery representation)
}
