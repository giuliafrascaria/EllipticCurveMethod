//
// Created by root on 04/03/17.
//

#ifndef ECM_DATASTRUCTURES_H
#define ECM_DATASTRUCTURES_H

//libreria per interi di grandi dimensioni
#include <gmp.h>

struct ECpoint
{
    mpz_t x;
    mpz_t y;
};

#endif //ECM_DATASTRUCTURES_H
