//
// Created by root on 04/03/17.
//

#include <stdio.h>
#include "dataStructures.h"
#include <unistd.h>
#include <stdlib.h>
#include "gc/include/gc/gc.h"


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
    //gmp_printf("\n\n\nx = %Zd , y= %Zd , z= %Zd \n\n\n", Q->X, Q->Y, Q->Z);
    //sleep(1);

    //gmp_printf("considero il primo %Zd\n", p);
    if(mpz_cmp_ui(p, 0) == 0)
    {
        //printf("sono qui\n");
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
        //printf("non sono andato all'infinito subito\n");
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

        /*mpz_set(res->X, Q->X);
        mpz_set(res->Y, Q->Y);
        mpz_set(res->Z, Q->Z);*/

        /*void mpz_realloc2 (mpz t x, mp bitcnt t n) [Function]
Change the space allocated for x to n bits. The value in x is preserved if it fits, or is set to
0 if not*/
        struct ECpoint term1, result;
        mpz_init(term1.X);
        mpz_init(term1.Y);
        mpz_init(term1.Z);


        mpz_init(result.X);
        mpz_init(result.Y);
        mpz_init(result.Z);

        //sleep(2);

        while(j >= 1 && (d->flag == 0))
        {

            //printf("doubling P\n");
            //res = doubleec(res, EC, pd, d, res);

            doubleec2(Q, EC, pd, d, &result);

            /*gmp_printf("resresresresQx %Zd\n", result.X);
            gmp_printf("resresresresy %Zd\n\n", result.Y);*/


            if(mpz_cmp_ui(result.X, 0) == 0)
            {
                if(mpz_cmp_ui(result.Y, 1) == 0)
                {
                    if(mpz_cmp_ui(result.Z, 0) == 0)
                    {
                        //printf("point to infinity for EZn, bad luck check 1\n");
                        //printf("ritorno 1\n");
                        return result;
                    }
                }
            }

            int mj, nj;
            mj = mpz_tstbit(m, j);
            nj = mpz_tstbit(n, j);

            //printf("%d, %d\n", mj, nj);

            if((mj == 1) && (nj == 0))
            {
                mpz_set(term1.X, result.X);
                mpz_set(term1.Y, result.Y);
                mpz_set(term1.Z, result.Z);


                //printf("Adding P and Q\n");
                //res = add(res, Q, EC, pd, d, res);
                add2(&term1, Q, EC, pd, d, &result);

                //free(res);
            }
            if((mj == 0) && (nj == 1))
            {
                mpz_set(term1.X, result.X);
                mpz_set(term1.Y, result.Y);
                mpz_set(term1.Z, result.Z);


                //printf("subtracting p and q\n");
                //res = sub(res, Q, EC, pd, d, res);
                sub2(&term1, Q, EC, pd, d, &result);

                //free(res);
            }


            mpz_set(Q->X, result.X);
            mpz_set(Q->Y, result.Y);
            mpz_set(Q->Z, result.Z);


            /*gmp_printf("bbbbbbbQx %Zd\n", Q->X);
            gmp_printf("bbbbbbbQy %Zd\n\n", Q->Y);*/


            if(mpz_cmp_ui(result.X, 0) == 0)
            {
                if(mpz_cmp_ui(result.Y, 1) == 0)
                {
                    if(mpz_cmp_ui(result.Z, 0) == 0)
                    {
                        printf("point to infinity for EZn, bad luck check 2\n");
                        //printf("ritorno 2\n");
                        return result;
                    }
                }
            }

            //bitwise operation
            //int mpz_tstbit (const mpz t op, mp bitcnt t bit_index) [Function]
            //Test bit bit index in op and return 0 or 1 accordingly.
            //if(mj, nj) == (1, 0) P = add(P, Q)
            //if(mj, nj) == (0, 1) p = sub(P, Q)
            //if(mpz_tstbit())
            j = j-1;
        }
        //printf("ritorno res\n");
        return result;
    }
}

struct ECpoint * add(struct ECpoint * P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint *res)
{

    if(mpz_cmp_ui(P->Z, 0) == 0)
    {
        //printf("null P\n");
        //sleep(1);
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
            //printf("point to infinity add1\n");
            //sleep();
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
    //struct ECpoint * PandQ = GC_MALLOC(sizeof(struct ECpoint));
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

void add2(struct ECpoint * P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint *res)
{

    /*gmp_printf("Px %Zd\n", P->X);
    gmp_printf("Py %Zd\n\n", P->Y);
    gmp_printf("Qx %Zd\n", Q->X);
    gmp_printf("Qy %Zd\n\n", Q->Y);*/


    //printf("heila\n");
    if(mpz_cmp_ui(P->Z, 0) == 0)
    {
        //printf("null P\n");
        mpz_set(res->X, Q->X);
        mpz_set(res->Y, Q->Y);
        mpz_set(res->Z, Q->Z);

        /*mpz_set(P->X, Q->X);
        mpz_set(P->Y, Q->Y);
        mpz_set(P->Z, Q->Z);*/

        //sleep(1);
        return;
    }

    if(mpz_cmp_ui(Q->Z, 0) == 0)
    {
        printf("null Q\n");
        //return P;

        mpz_set(res->X, P->X);
        mpz_set(res->Y, P->Y);
        mpz_set(res->Z, P->Z);
        return;
    }
    mpz_t m;

    if(mpz_cmp(P->X, Q->X) == 0)
    {

        /*gmp_printf("resx %Zd\n", res->X);
        gmp_printf("resy %Zd\n\n", res->Y);*/

        mpz_t sumOfY;
        mpz_init(sumOfY);

        mpz_add(sumOfY, P->Y, Q->Y);
        //gmp_printf("sumofY = %Zd\n\n", sumOfY);
        if(mpz_cmp_ui(sumOfY, 0) == 0)
        {
            //infinity
            printf("point to infinity add\n");
            //sleep(1);
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

            mpz_set_ui(P->X, 0);
            mpz_set_ui(P->Y, 1);
            mpz_set_ui(P->Z, 0);

            return;
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
            printf("found non invertible den add2\n");
            d->flag = 1;
            mpz_set(d->d, doubley);
            //I should break the cycle
            return;
        }

        mpz_mul(m, firstterm, invertY);

        //gmp_printf("valore di m %Zd\n", m);

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
            return;
        }
        mpz_mul(m, firstterm, invertsecond);
        mpz_clear(firstterm);
        mpz_clear(invertsecond);
        mpz_clear(secondterm);
    }
    //printf("calculating x3\n");
    //calculate x3
    //mp_bitcnt_t size = mpz_size(pd.n);
    //mpz_realloc2(PandQ->Y, size);
    //struct ECpoint * PandQ = malloc(sizeof(struct ECpoint));        //unica malloc che rimane, unica malloc che riempie
    //struct ECpoint * PandQ = GC_MALLOC(sizeof(struct ECpoint));


    mpz_t squarem, partial, x3;
    mpz_init(squarem);
    mpz_init(partial);
    mpz_init(x3);

    mpz_pow_ui(squarem, m, 2);
    mpz_sub(partial, squarem, P->X);
    //mpz_sub(P->X, partial, Q->X);
    mpz_sub(x3, partial, Q->X);

   // gmp_printf("X3 = %Zd\nPx = %Zd\n", x3, P->X);



    mpz_clear(partial);
    mpz_clear(squarem);

    mpz_mod(res->X, x3, pd.n);

    mpz_clear(x3);


    //gmp_printf("resX after mod = %Zd\nPx = %Zd\n", res->X, P->X);

    //calculate y3
    //m(x3 - x1)
    //printf("calculating y3\n");
    mpz_t partial2;
    mpz_init(partial2);

    //mpz_sub(partial2, P->X, P->X);
    //mpz_mul(P->Y, m, partial2);

    mpz_sub(partial2, res->X, P->X);
    //gmp_printf("partial2 = %Zd\n", partial2);
    mpz_mul(res->Y, m, partial2);


    mpz_clear(partial2);
    mpz_clear(m);

    //mpz_mod(P->Y, P->Y, pd.n);
    mpz_mod(res->Y, res->Y, pd.n);

    /*mp_bitcnt_t size = mpz_size(pd.n);
    mpz_realloc2(PandQ->Y, size);*/

    //mpz_set_ui(P->Z, 1);
    //mpz_set(P->X, res->X);
    //mpz_set(P->Y, res->Y);
    mpz_set_ui(res->Z, 1);

    //free(P);
    //sleep(1);
    return;
    //return (x3, y3)
}


struct ECpoint * sub(struct ECpoint *P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res)
{
    struct ECpoint Qnegate = negate(Q);
    return add(P, &Qnegate, EC, pd, d, res);

}

void sub2(struct ECpoint *P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res)
{
    struct ECpoint Qnegate = negate(Q);
    //printf("I negated QY!!\n\n\n");
    return add2(P, &Qnegate, EC, pd, d, res);

}


struct ECpoint * doubleec(struct ECpoint * P, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD *d, struct ECpoint *res)
{
    return add(P, P, EC, pd, d, res);
}

void doubleec2(struct ECpoint * P, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD *d, struct ECpoint *res)
{
    /*gmp_printf("Px %Zd\n", P->X);
    gmp_printf("Py %Zd\n\n", P->Y);*/

    return add2(P, P, EC, pd, d, res);
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

void montgomeryLadder()
{
    //efficient calculation of kP
}


