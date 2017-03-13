//
// Created by root on 04/03/17.
//

#include <stdio.h>
#include "dataStructures.h"

struct ECpoint negate(struct ECpoint P, struct nonInvertibleD d);
struct ECpoint doubleec(struct ECpoint P, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD d);
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

struct ECpoint ECmultiplyTraditional(struct ECpoint Q, mpz_t p, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD d)
{
    //initialize
    printf("initializing\n");
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
        printf("infinity\n");

        return infinity;
    }
    else
    {
        //compare bits of [3n, n]
        printf("starting ladder\n");
        mpz_t loopIndex;
        mpz_sub_ui(loopIndex, pd.stageOneB, 2);

        struct ECpoint P;

        while(mpz_cmp_ui(loopIndex, 1) >= 0)
        {

            P = doubleec(Q, EC, pd, d);
            //bitwise operation
                //int mpz_tstbit (const mpz t op, mp bitcnt t bit_index) [Function]
                //Test bit bit index in op and return 0 or 1 accordingly.
            //if(mj, nj) == (1, 0) P = add(P, Q)
            //if(mj, nj) == (0, 1) p = sub(P, Q)
            //if(mpz_tstbit())
            mpz_sub_ui(loopIndex, loopIndex, 1);
        }
        return P;
    }
}

struct ECpoint add(struct ECpoint P, struct ECpoint Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD d)
{
    if(mpz_cmp_ui(P.Z, 0) == 0)
    {
        return Q;
    }
    if(mpz_cmp_ui(Q.Z, 0) == 0)
    {
        return P;
    }
    mpz_t m;
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
        mpz_t squareX, threeSquareX, firstterm, doubley, invertY;
        mpz_init(squareX);
        mpz_init(threeSquareX);
        mpz_init(firstterm);
        mpz_init(doubley);
        mpz_init(invertY);
        mpz_init(m);

        mpz_pow_ui(squareX, P.X, 2);
        mpz_mul_ui(threeSquareX, squareX, 3);
        mpz_add(firstterm, threeSquareX, EC.a);

        mpz_mul_ui(doubley, P.Y, 2);
        int result;
        result = mpz_invert(invertY, doubley, pd.n);     //if the return value is zero the invert is not defined
        if(result == 0)
        {
            //non invertible d
            d.flag = 1;
            mpz_set(d.d, doubley);
            //I should break the cycle
            return P;
        }

        mpz_mul(m, firstterm, invertY);
    }
    else
    {
        //calculate m
        mpz_t firstterm, secondterm, invertsecond;
        mpz_init(firstterm);
        mpz_init(secondterm);
        mpz_init(invertsecond);
        mpz_init(m);

        mpz_sub(firstterm, Q.Y, P.Y);
        mpz_sub(secondterm, Q.X, P.X);
        int result;
        result = mpz_invert(invertsecond, secondterm, pd.n);
        if(result == 0)
        {
            //non invertible d
            d.flag = 1;
            mpz_set(d.d, secondterm);
            //I should break the cycle
            return P;
        }
        mpz_mul(m, firstterm, invertsecond);
    }
    //calculate x3
    struct ECpoint PandQ;
    mpz_init(PandQ.X);
    mpz_init(PandQ.Y);
    mpz_init(PandQ.Z);

    mpz_t squarem, partial;
    mpz_init(squarem);
    mpz_init(partial);

    mpz_pow_ui(squarem, m, 2);
    mpz_sub(partial, squarem, P.X);
    mpz_sub(PandQ.X, partial, Q.X);

    //calculate y3
    //m(x3 - x1)
    mpz_t partial2;
    mpz_init(partial2);

    mpz_sub(partial2, PandQ.X, P.X);
    mpz_mul(PandQ.Y, m, partial2);

    mpz_set_ui(PandQ.Z, 1);

    return PandQ;
    //return (x3, y3)
}

struct ECpoint sub(struct ECpoint P, struct ECpoint Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD d)
{
    struct ECpoint Qnegate = negate(Q, d);
    return add(P, Qnegate, EC, pd, d);
}

struct ECpoint doubleec(struct ECpoint P, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD d)
{
    return add(P, P, EC, pd, d);
}

struct ECpoint negate(struct ECpoint P, struct nonInvertibleD d)
{
    //return (x -y z)
    //void mpq_neg (mpq t negated_operand, const mpq t operand)
    struct ECpoint Pnegate;
    mpz_init(Pnegate.X);
    mpz_init(Pnegate.Y);
    mpz_init(Pnegate.Z);

    mpz_set(Pnegate.X, P.X);
    mpz_set(Pnegate.Z, P.Z);
    mpz_neg(Pnegate.Y, P.Y);

    return Pnegate;
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


