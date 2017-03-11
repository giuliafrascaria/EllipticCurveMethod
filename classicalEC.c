//
// Created by root on 10/03/17.
//


#include "dataStructures.h"

void randomECtraditional(struct weirstrassEC EC, struct ECpoint Q, struct problemData pd)
{
    //random sigma, derive u v anc C from this

    //initialize random state
    gmp_randstate_t state;
    gmp_randinit_mt(state);

    //generate random x, y, a in [0, n-1]

    mpz_urandomm(Q.X, state, pd.n);
    mpz_urandomm(Q.Y, state, pd.n);
    mpz_urandomm(EC.a, state, pd.n);
    //the curve is determined by these values


    mpz_t squareY, cubeX, aX;
    mpz_init(squareY);
    mpz_init(cubeX);
    mpz_init(aX);
    mpz_pow_ui(squareY, Q.Y, 2);
    mpz_pow_ui(cubeX, Q.X, 3);
    mpz_mul(aX, Q.X, EC.a);


    mpz_t temp;
    mpz_init(temp);
    mpz_sub(temp, squareY, cubeX);
    mpz_sub(EC.b, temp, aX);
    mpz_mod(EC.b, EC.b, pd.n);

}

