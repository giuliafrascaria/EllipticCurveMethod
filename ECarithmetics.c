//
// Created by root on 04/03/17.
//

#include "dataStructures.h"

/*
 * definire qui le operazioni di somma sulla curva ellittica, calcolo del coefficiente angolare m,
 *
 * doubleh
 * addh
 * elliptic multiplication (Montgomery method)
 * */

struct ECpoint ECmultiplyMontgomery(struct ECpoint Q, mpz_t p)
{
    //initialization
    if(mpz_cmp_ui(p, 0) == 0)
    {
        struct ECpoint infinity;    //the point to infinity is recognized with the pair 0,0
        mpz_init(infinity.X);
        mpz_set_ui(infinity.X, 0);

        mpz_init(infinity.Z);
        mpz_set_ui(infinity.Z, 0);
        return infinity;
    }
    else if(mpz_cmp_ui(p, 1) == 0)
    {
        return Q;
    }
    else if(mpz_cmp_ui(p, 2) == 0)
    {
        //doubleh
        return Q;       //delete
    }
    //begin Montgomery adding ladder

    //loop over bits

    //final calculation

}

struct ECpoint ECmultiplyTraditional(struct ECpoint Q, mpz_t p, struct problemData pd)
{
    //initialize
    if(mpz_cmp_ui(p, 0))
    {
        //point to infinity
        struct ECpoint infinity;
        mpz_init(infinity.X);
        mpz_init(infinity.Y);
        mpz_init(infinity.Z);

        mpz_set_ui(infinity.X, 0);
        mpz_set_ui(infinity.Y, 1);
        mpz_set_ui(infinity.Z, 0);

        return infinity;
    }
    else
    {
        //compare bits of [3n, n]
        mpz_t loopIndex;
        mpz_sub_ui(loopIndex, pd.stageOneB, 2);

        struct ECpoint P;

        while(mpz_cmp_ui(loopIndex, 1) >= 0)
        {
            //P = double(P)
            //bitwise operation
                //int mpz_tstbit (const mpz t op, mp bitcnt t bit_index) [Function]
                //Test bit bit index in op and return 0 or 1 accordingly.
            //if(mj, nj) == (1, 0) P = add(P, Q)
            //if(mj, nj) == (0, 1) p = sub(P, Q)
            mpz_sub_ui(loopIndex, loopIndex, 1);
        }
        return P;
    }
}

struct ECpoint add(struct ECpoint P, struct ECpoint Q)
{
    if(mpz_cmp_ui(P.Z, 0) == 0)
    {
        return Q;
    }
    if(mpz_cmp_ui(Q.Z, 0) == 0)
    {
        return P;
    }
    if(mpz_cmp(P.X, Q.X) == 0)
    {
        mpz_t sumOfY;
        mpz_init(sumOfY);
        mpz_add(sumOfY, P.Y, Q.Y);
        if(mpz_cmp_ui(sumOfY, 0))
        {
            //infinity
            struct ECpoint infinity;
            mpz_init(infinity.X);
            mpz_init(infinity.Y);
            mpz_init(infinity.Z);

            mpz_set_ui(infinity.X, 0);
            mpz_set_ui(infinity.Y, 1);
            mpz_set_ui(infinity.Z, 0);

            return infinity;
        }
        //calculate m
    }
    else
    {
        //calculate m
    }
    //calculate x3
    //return (x3, y3)
}

void sub(struct ECpoint P, struct ECpoint Q)
{

}

struct ECpoint doubleec(struct ECpoint P)
{
    return add(P, P);
}

void addh()
{

}

void doubleh()
{

}

void doubleOP()
{

}


