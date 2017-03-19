//
// Created by root on 04/03/17.
//

#include <stdio.h>
#include "dataStructures.h"
#include <unistd.h>
#include <stdlib.h>
//#include <gc/gc.h>


struct ECpoint * add(struct ECpoint *P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res);
struct ECpoint * sub(struct ECpoint *P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res);
/*
 * definire qui le operazioni di somma sulla curva ellittica, calcolo del coefficiente angolare m,
 *
 * doubleh
 * addh
 * elliptic multiplication (Montgomery method)
 * */

/*
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
*/

struct ECpoint ECmultiplyTraditional(struct ECpoint * Q, mpz_t p, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res)
{
    //initialize
    //printf("initializing ladder\n");
    //gmp_printf("x = %Zd , y= %Zd , z= %Zd \n", Q->X, Q->Y, Q->Z);
    //sleep(3);

    //gmp_printf("considero il primo %Zd\n", p);
    if(mpz_cmp_ui(p, 0) == 0)
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
        //printf("starting ladder\n");
        mpz_t n;
        mpz_init(n);
        mpz_set(n, pd.n);

        mpz_t m;
        mpz_init(m);
        mpz_mul_ui(m, n, 3);


        //size_t mpz_sizeinbase (const mpz t op, int base)
        ssize_t size = mpz_sizeinbase(m, 2);
        //printf("size of m = %d\n", size);

        mpz_realloc2(n, (mp_bitcnt_t) size);

        mp_bitcnt_t j = (mp_bitcnt_t) size-2;

        /*struct ECpoint * P = malloc(sizeof(struct ECpoint));
        mpz_init(P->X);
        mpz_init(P->Y);
        mpz_init(P->Z);

        mpz_set(P->X, Q->X);
        mpz_set(P->Y, Q->Y);
        mpz_set(P->Z, Q->Z);*/

        mpz_set(res->X, Q->X);
        mpz_set(res->Y, Q->Y);
        mpz_set(res->Z, Q->Z);

        /*void mpz_realloc2 (mpz t x, mp bitcnt t n) [Function]
Change the space allocated for x to n bits. The value in x is preserved if it fits, or is set to
0 if not*/

        while(j >= 1 && (d->flag == 0))
        {
            //printf("doubling P\n");
            res = doubleec(res, EC, pd, d, res);
            int mj, nj;
            mj = mpz_tstbit(m, j);
            nj = mpz_tstbit(n, j);

            //printf("%d, %d\n", mj, nj);

            if((mj == 1) && (nj == 0))
            {
                //printf("Adding P and Q\n");
                res = add(res, Q, EC, pd, d, res);
            }
            if((mj == 0) && (nj == 1))
            {
                //printf("subtracting p and q\n");
                res = sub(res, Q, EC, pd, d, res);
            }
            //bitwise operation
                //int mpz_tstbit (const mpz t op, mp bitcnt t bit_index) [Function]
                //Test bit bit index in op and return 0 or 1 accordingly.
            //if(mj, nj) == (1, 0) P = add(P, Q)
            //if(mj, nj) == (0, 1) p = sub(P, Q)
            //if(mpz_tstbit())
            j = j-1;
        }
        return *res;
    }
}

struct ECpoint * add(struct ECpoint * P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint *res)
{

    if(mpz_cmp_ui(P->Z, 0) == 0)
    {
        printf("null P\n");
        return Q;
    }

    if(mpz_cmp_ui(Q->Z, 0) == 0)
    {
        printf("null Q\n");
        return P;
    }
    mpz_t m;

    if(mpz_cmp(P->X, Q->X) == 0)
    {

   /*     gmp_printf("Qx %Zd\n", Q->X);
        gmp_printf("Qy %Zd\n", Q->Y);
        gmp_printf("Px %Zd\n", P->X);
        gmp_printf("Py %Zd\n", P->Y);
*/

        mpz_t sumOfY;
        mpz_init(sumOfY);

        mpz_add(sumOfY, P->Y, Q->Y);
        if(mpz_cmp_ui(sumOfY, 0) == 0)
        {
            //infinity
            printf("point to infinity\n");
            /*struct ECpoint * infinity = malloc(sizeof(struct ECpoint));
            mpz_init(infinity->X);
            mpz_init(infinity->Y);
            mpz_init(infinity->Z);

            mpz_set_ui(infinity->X, 0);
            mpz_set_ui(infinity->Y, 1);
            mpz_set_ui(infinity->Z, 0);

            return infinity;*/
            mpz_set_ui(res->X, 0);
            mpz_set_ui(res->Y, 1);
            mpz_set_ui(res->Z, 0);

            return res;
        }

        mpz_clear(sumOfY);
        //calculate m
        //printf("calculating m\n");
        mpz_t squareX, threeSquareX, firstterm, doubley, invertY;
        mpz_init(squareX);
        mpz_init(threeSquareX);
        mpz_init(firstterm);
        mpz_init(doubley);
        mpz_init(invertY);
        mpz_init(m);

        mpz_pow_ui(squareX, P->X, 2);
        mpz_mul_ui(threeSquareX, squareX, 3);
        mpz_add(firstterm, threeSquareX, EC.a);

        mpz_mul_ui(doubley, P->Y, 2);
        int result;
        result = mpz_invert(invertY, doubley, pd.n);     //if the return value is zero the invert is not defined
        if(result == 0)
        {
            //non invertible d
            printf("found non invertible den\n");
            d->flag = 1;
            mpz_set(d->d, doubley);
            //I should break the cycle
            return P;
        }

        mpz_mul(m, firstterm, invertY);

        mpz_clear(firstterm);
        mpz_clear(invertY);
        mpz_clear(doubley);
        mpz_clear(threeSquareX);
        mpz_clear(squareX);
    }
    else
    {
        //calculate m
        //printf("calculating m 1\n");
        mpz_t firstterm, secondterm, invertsecond;
        mpz_init(firstterm);
        mpz_init(secondterm);
        mpz_init(invertsecond);
        mpz_init(m);

        mpz_sub(firstterm, Q->Y, P->Y);
        mpz_sub(secondterm, Q->X, P->X);
        int result;
        result = mpz_invert(invertsecond, secondterm, pd.n);
        if(result == 0)
        {
            //non invertible d
            printf("noninvertible den 1\n");
            d->flag = 1;
            mpz_set(d->d, secondterm);

            //I should break the cycle
            return P;
        }
        mpz_mul(m, firstterm, invertsecond);
        mpz_clear(firstterm);
        mpz_clear(invertsecond);
        mpz_clear(secondterm);
    }
    //printf("calculating x3\n");
    //calculate x3
    mp_bitcnt_t size = mpz_size(pd.n);
    //mpz_realloc2(PandQ->Y, size);
    struct ECpoint * PandQ = malloc(sizeof(struct ECpoint));        //unica malloc che rimane, unica malloc che riempie
    mpz_init2(PandQ->X, size);
    mpz_init2(PandQ->Y, size);
    mpz_init2(PandQ->Z, size);

    mpz_t squarem, partial;
    mpz_init(squarem);
    mpz_init(partial);

    mpz_pow_ui(squarem, m, 2);
    mpz_sub(partial, squarem, P->X);
    mpz_sub(PandQ->X, partial, Q->X);

    mpz_clear(partial);
    mpz_clear(squarem);

    mpz_mod(PandQ->X, PandQ->X, pd.n);



    //calculate y3
    //m(x3 - x1)
    //printf("calculating y3\n");
    mpz_t partial2;
    mpz_init(partial2);

    mpz_sub(partial2, PandQ->X, P->X);
    mpz_mul(PandQ->Y, m, partial2);

    mpz_clear(partial2);
    mpz_clear(m);

    mpz_mod(PandQ->Y, PandQ->Y, pd.n);

    /*mp_bitcnt_t size = mpz_size(pd.n);
    mpz_realloc2(PandQ->Y, size);*/

    mpz_set_ui(PandQ->Z, 1);

    //free(P);

    return PandQ;
    //return (x3, y3)
}

struct ECpoint * sub(struct ECpoint *P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res)
{
    struct ECpoint Qnegate = negate(Q);
    return add(P, &Qnegate, EC, pd, d, res);

}

struct ECpoint * doubleec(struct ECpoint * P, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD *d, struct ECpoint *res)
{
    return add(P, P, EC, pd, d, res);
}

struct ECpoint negate(struct ECpoint *P)
{
    //return (x -y z)
    //void mpq_neg (mpq t negated_operand, const mpq t operand)
    /*struct ECpoint Pnegate;
    mpz_init(Pnegate.X);
    mpz_init(Pnegate.Y);
    mpz_init(Pnegate.Z);*/

    /*mpz_set(Pnegate.X, P->X);
    mpz_set(Pnegate.Z, P->Z);
    mpz_neg(Pnegate.Y, P->Y);*/

    
    mpz_neg(P->Y, P->Y);

    //return Pnegate;
    return *P;
}




