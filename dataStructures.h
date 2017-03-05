//
// Created by root on 04/03/17.
//

#ifndef ECM_DATASTRUCTURES_H
#define ECM_DATASTRUCTURES_H

//library for large integer operations
#include <gmp.h>

struct ECpoint
{
    //montgomery coordinates [X:Z]
    mpz_t X;
    mpz_t Z;
};

struct problemData
{
    mpz_t n;                //large composite number to be factored

    mpz_t stageOneB;        //must be even
    mpz_t stageTwoB;        //100B1
    mpz_t D;                //total memory
};

struct ellipticCurve
{
    //not weirstrass form, this is the invertionless algorithm
    mpz_t sigma;
    mpz_t u;
    mpz_t v;
    mpz_t C;

};

#endif //ECM_DATASTRUCTURES_H
