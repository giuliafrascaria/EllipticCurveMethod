//
// Created by root on 04/03/17.
//

#include <stdio.h>
#include "dataStructures.h"
#include <unistd.h>
#include <stdlib.h>


struct ECpoint * sub(struct ECpoint *P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res);
/*
 * definire qui le operazioni di somma sulla curva ellittica, calcolo del coefficiente angolare m,
 *
 * doubleh
 * addh
 * elliptic multiplication (Montgomery method)
 * */

struct ECpoint doubleAndAdd(struct ECpoint * P, mpz_t p,  struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res)
{
    //usare l'espansione binaria di p

    if(mpz_cmp_ui(P->X, 0) == 0)
    {
        printf("point to infinity doubleandAdd\n\n");
        return *P;
    }

    char * binaryp = malloc(sizeof(*p));
    //char * binaryP = mpz_get_str(binaryP, 2, p);
    binaryp = mpz_get_str(binaryp, 2, p);

    mpz_t binP;
    mpz_init(binP);

    //if(mpz_set_str(binP, binaryP, 2) != 0 )
    if(mpz_set_str(binP, binaryp, 2) != 0 )
    {
        printf("failed to convert\n");
        exit(EXIT_FAILURE);
    }

    //sleep(1);

    struct ECpoint Q, R, res1, res2;
    mpz_init(Q.X);
    mpz_init(Q.Y);
    mpz_init(Q.Z);

    mpz_init(R.X);
    mpz_init(R.Y);
    mpz_init(R.Z);

    mpz_init(res1.X);
    mpz_init(res1.Y);
    mpz_init(res1.Z);

    mpz_init(res2.X);
    mpz_init(res2.Y);
    mpz_init(res2.Z);

    mpz_set(Q.X, P->X); //Q = P0
    mpz_set(Q.Y, P->Y);
    mpz_set(Q.Z, P->Z);

    doubleec(P, EC, pd, d, &R); //R = 2P
    //checkIfCurve(R, EC, pd);


    ssize_t size = mpz_sizeinbase(binP, 2);

    mp_bitcnt_t j = (mp_bitcnt_t) size-2;
    long index = j;

    while(index >= 0 && (d->flag == 0))
    {

        int pj;
        pj = mpz_tstbit(binP, j);

        if(pj == 1)
        {
            add2(&R, &Q, EC, pd, d, &res1);    //Q = Q + R
            doubleec(&R, EC, pd, d, &res2) ;     // R = 2R
            printf("bla\n");
        }
        else
        {
            doubleec(&Q, EC, pd, d, &res1);      //Q = 2Q
            add2(&Q, &R, EC, pd, d, &res2);    //R = Q + R
            printf("bla2\n");
        }
        //sleep(1);

        mpz_set(Q.X, res1.X);
        mpz_set(Q.Y, res1.Y);
        mpz_set(Q.Z, res1.Z);

        mpz_set(R.X, res2.X);
        mpz_set(R.Y, res2.Y);
        mpz_set(R.Z, res2.Z);

        j = j-1;
        index = index-1;

        if(mpz_cmp_ui(res1.X, 0) == 0)
        {
            /*if(mpz_cmp_ui(res1.Y, 1) == 0)
            {
                if(mpz_cmp_ui(res1.Z, 0) == 0)
                {
                    //printf("point to infinity for EZn, bad luck check 1\n");
                    //printf("ritorno 1\n");
                    return;
                }
            }*/
            return res1;
        }

        if(mpz_cmp_ui(res2.X, 0) == 0)
        {
            /*if(mpz_cmp_ui(res2.Y, 1) == 0)
            {
                if(mpz_cmp_ui(res2.Z, 0) == 0)
                {
                    //printf("point to infinity for EZn, bad luck check 1\n");
                    //printf("ritorno 1\n");
                    return;
                }
            }*/
            return res2;
        }


    }

    mpz_set(res->X, Q.X);
    mpz_set(res->Y, Q.Y);
    mpz_set(res->Z, Q.Z);

    mpz_clear(res1.X);
    mpz_clear(res1.Y);
    mpz_clear(res1.Z);

    mpz_clear(res2.X);
    mpz_clear(res2.Y);
    mpz_clear(res2.Z);

    /*mpz_clear(Q.X);
    mpz_clear(Q.Y);
    mpz_clear(Q.Z);*/

    mpz_clear(R.X);
    mpz_clear(R.Y);
    mpz_clear(R.Z);

    return Q;

}

struct ECpoint ECmultiplyTraditional(struct ECpoint * Q, mpz_t p, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res)
{

    struct ECpoint P;
    mpz_init(P.X);
    mpz_init(P.Y);
    mpz_init(P.Z);

    mpz_set(P.X, Q->X);
    mpz_set(P.Y, Q->Y);
    mpz_set(P.Z, Q->Z);

    if(mpz_cmp_ui(p, 0) == 0)
    {

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

        mpz_t n;
        mpz_init(n);
        mpz_set(n, p); //!!!!!!!!

        mpz_t m;
        mpz_init(m);
        mpz_mul_ui(m, n, 3);
        size_t binarySizeM = mpz_sizeinbase(m, 2);

        char * binaryn = NULL;
        binaryn = mpz_get_str(binaryn, 2, n);
        char * binarym = NULL;
        binarym = mpz_get_str(binarym, 2, m);

//        printf("M : %d\n",(int) binarySizeM);
//        printf("N : %d\n",(int) mpz_sizeinbase(n, 2));
//        for(int i = 0; i < mpz_sizeinbase(p, 2); i++){
/*
            printf("binaryn = %s\n", binaryn);
        }
*/


        ssize_t size = mpz_sizeinbase(m, 2);

        mpz_realloc2(n, (mp_bitcnt_t) size);

//        mp_bitcnt_t j = (mp_bitcnt_t) size-2;


        struct ECpoint term1, result;
        mpz_init(term1.X);
        mpz_init(term1.Y);
        mpz_init(term1.Z);


        mpz_init(result.X);
        mpz_init(result.Y);
        mpz_init(result.Z);

        for(int i = (int) (binarySizeM - 2); i >= 1; i--){

            doubleec2(&P, EC, pd, d, &result);
            checkIfCurve(result, EC, pd);
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
            char ni = binaryn[i];
            char mi = binarym[i];
            if(mi == '1' && ni == '0')
            {
//                printf("mi = %c , ni = %c\n", mi, ni);
                mpz_set(term1.X, result.X);
                mpz_set(term1.Y, result.Y);
                mpz_set(term1.Z, result.Z);
                add2(&term1, Q, EC, pd, d, &result);
                //result = add3(&term1, Q, EC, pd, d);
                checkIfCurve(result, EC, pd);
            }
            else if(mi == '0' && ni == '1')
            {
//                printf("mi = %c , ni = %c\n", mi, ni);
                mpz_set(term1.X, result.X);
                mpz_set(term1.Y, result.Y);
                mpz_set(term1.Z, result.Z);
                sub2(&term1, Q, EC, pd, d, &result);
            }
            mpz_set(P.X, result.X);
            mpz_set(P.Y, result.Y);
            mpz_set(P.Z, result.Z);
        }

        return result;
    }
}

struct ECpoint ECmultiplyTraditional2(struct ECpoint * Q, mpz_t *p, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res)
{

    //gmp_printf("p = %Zd\n", *p);

    struct ECpoint P;
    mpz_init(P.X);
    mpz_init(P.Y);
    mpz_init(P.Z);

    mpz_set(P.X, Q->X);
    mpz_set(P.Y, Q->Y);
    mpz_set(P.Z, Q->Z);

    if(mpz_cmp_ui(*p, 0) == 0)
    {

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

        mpz_t n;
        mpz_init(n);
        mpz_set(n, *p); //!!!!!!!!

        mpz_t m;
        mpz_init(m);
        mpz_mul_ui(m, n, 3);
        size_t binarySizeM = mpz_sizeinbase(m, 2);

        char * binaryn = NULL;
        binaryn = mpz_get_str(binaryn, 2, n);
        char * binarym = NULL;
        binarym = mpz_get_str(binarym, 2, m);

//        printf("M : %d\n",(int) binarySizeM);
//        printf("N : %d\n",(int) mpz_sizeinbase(n, 2));
//        for(int i = 0; i < mpz_sizeinbase(p, 2); i++){
/*
            printf("binaryn = %s\n", binaryn);
        }
*/


        ssize_t size = mpz_sizeinbase(m, 2);

        mpz_realloc2(n, (mp_bitcnt_t) size);

//        mp_bitcnt_t j = (mp_bitcnt_t) size-2;


        struct ECpoint term1, result;
        mpz_init(term1.X);
        mpz_init(term1.Y);
        mpz_init(term1.Z);


        mpz_init(result.X);
        mpz_init(result.Y);
        mpz_init(result.Z);

        for(int i = (int) (binarySizeM - 2); i >= 1; i--){

            doubleec2(&P, EC, pd, d, &result);
            checkIfCurve(result, EC, pd);
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
            char ni = binaryn[i];
            char mi = binarym[i];
            if(mi == '1' && ni == '0')
            {
//                printf("mi = %c , ni = %c\n", mi, ni);
                mpz_set(term1.X, result.X);
                mpz_set(term1.Y, result.Y);
                mpz_set(term1.Z, result.Z);
                add2(&term1, Q, EC, pd, d, &result);
                //result = add3(&term1, Q, EC, pd, d);
                checkIfCurve(result, EC, pd);
            }
            else if(mi == '0' && ni == '1')
            {
//                printf("mi = %c , ni = %c\n", mi, ni);
                mpz_set(term1.X, result.X);
                mpz_set(term1.Y, result.Y);
                mpz_set(term1.Z, result.Z);
                sub2(&term1, Q, EC, pd, d, &result);
            }
            mpz_set(P.X, result.X);
            mpz_set(P.Y, result.Y);
            mpz_set(P.Z, result.Z);
        }

        return result;
    }
}

struct ECpoint doubleAndAdd2(struct ECpoint * P, mpz_t p, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d)
{
    struct ECpoint Q, N;
    mpz_init(Q.X);
    mpz_init(Q.Y);
    mpz_init(Q.Z);

    mpz_init(N.X);
    mpz_init(N.Y);
    mpz_init(N.Z);

    mpz_set(N.X, P->X);
    mpz_set(N.Y, P->Y);
    mpz_set(N.Z, P->Z);

    mpz_set_ui(Q.X, 0);
    mpz_set_ui(Q.Y, 1);
    mpz_set_ui(Q.Z, 0);


    size_t binarySizep = mpz_sizeinbase(p, 2);

    char * binaryp = NULL;
    binaryp = mpz_get_str(binaryp, 2, p);   //ha un bit in più per il segno lie

    struct ECpoint result1, result2;



    mpz_init(result1.X);
    mpz_init(result1.Y);
    mpz_init(result1.Z);

    mpz_init(result2.X);
    mpz_init(result2.Y);
    mpz_init(result2.Z);

    for(int i = 0; i < (int) binarySizep; i++)
    {
        if(binaryp[i] == '1')
        {
            add2(&Q, &N, EC, pd, d, &result1);
        }
        doubleec2(&N, EC, pd, d, &result2);

        mpz_set(Q.X, result1.X);
        mpz_set(Q.Y, result1.Y);
        mpz_set(Q.Z, result1.Z);

        mpz_set(N.X, result2.X);
        mpz_set(N.Y, result2.Y);
        mpz_set(N.Z, result2.Z);

    }

    return Q;

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

    if(mpz_cmp_ui(P->Z, 0) == 0)
    {
        printf("null P\n");
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
        //printf("tangente verticale\n");


        //m=(3x1^2 + a)/2y1

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
        //printf("tangente non verticale\n");

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

    //mpz_mod(res->X, x3, pd.n);
    mpz_set(res->X, x3);
    mpz_mod(res->X, res->X, pd.n);

    mpz_clear(x3);


    //gmp_printf("resX after mod = %Zd\nPx = %Zd\n", res->X, P->X);

    //calculate y3
    //m(x1 - x3) - y1
    //printf("calculating y3\n");
    mpz_t partial2, partial3;
    mpz_init(partial2);
    mpz_init(partial3);

    //mpz_sub(partial2, P->X, P->X);
    //mpz_mul(P->Y, m, partial2);

    //mpz_sub(partial2, res->X, P->X);
    mpz_sub(partial2,  P->X, res->X);

    mpz_mul(partial3, m, partial2);

    mpz_sub(res->Y, partial3, P->Y);


    //gmp_printf("partial2 = %Zd\n", partial2);
    //mpz_mul(res->Y, m, partial2);


    mpz_clear(partial2);
    mpz_clear(partial3);
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

struct ECpoint add3(struct ECpoint * P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d)
{


    /*gmp_printf("Px %Zd\n", P->X);
    gmp_printf("Py %Zd\n\n", P->Y);
    gmp_printf("Qx %Zd\n", Q->X);
    gmp_printf("Qy %Zd\n\n", Q->Y);*/

    if(mpz_cmp_ui(P->Z, 0) == 0)
    {
        struct ECpoint Qr;
        mpz_init(Qr.X);
        mpz_init(Qr.Y);
        mpz_init(Qr.Z);

        mpz_set(Qr.X, Q->X);
        mpz_set(Qr.Y, Q->Y);
        mpz_set(Qr.Z, Q->Z);
        printf("Q\n");

        return Qr;
    }

    if(mpz_cmp_ui(Q->Z, 0) == 0)
    {
        struct ECpoint Pr;
        mpz_init(Pr.X);
        mpz_init(Pr.Y);
        mpz_init(Pr.Z);

        mpz_set(Pr.X, P->X);
        mpz_set(Pr.Y, P->Y);
        mpz_set(Pr.Z, P->Z);
        printf("P\n");

        return Pr;
    }
    mpz_t m;
    struct ECpoint result;
    mpz_init(result.X);
    mpz_init(result.Y);
    mpz_init(result.Z);

    if(mpz_cmp(P->X, Q->X) == 0)
    {
        //printf("tangente verticale\n");


        //m=(3x1^2 + a)/2y1

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



            mpz_set_ui(result.X, 0);
            mpz_set_ui(result.Y, 1);
            mpz_set_ui(result.Z, 0);

            /*mpz_set_ui(P->X, 0);
            mpz_set_ui(P->Y, 1);
            mpz_set_ui(P->Z, 0);*/

            return result;
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
        int resultval;
        resultval = mpz_invert(invertY, doubley, pd.n);     //if the return value is zero the invert is not defined
        if(resultval == 0)
        {
            //non invertible d
            printf("found non invertible den add2\n");
            d->flag = 1;
            mpz_set(d->d, doubley);
            //I should break the cycle
            return result;
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
        //printf("tangente non verticale\n");

        mpz_t firstterm, secondterm, invertsecond;
        mpz_init(firstterm);
        mpz_init(secondterm);
        mpz_init(invertsecond);
        mpz_init(m);

        mpz_sub(firstterm, Q->Y, P->Y);
        mpz_sub(secondterm, Q->X, P->X);
        int resultval;
        resultval = mpz_invert(invertsecond, secondterm, pd.n);
        if(resultval == 0)
        {
            //non invertible d
            printf("noninvertible den 1\n");
            d->flag = 1;
            mpz_set(d->d, secondterm);

            //I should break the cycle
            return result;
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

    //mpz_mod(res->X, x3, pd.n);
    mpz_set(result.X, x3);
    mpz_mod(result.X, result.X, pd.n);

    mpz_clear(x3);


    //gmp_printf("resX after mod = %Zd\nPx = %Zd\n", res->X, P->X);

    //calculate y3
    //m(x1 - x3) - y1
    //printf("calculating y3\n");
    mpz_t partial2, partial3;
    mpz_init(partial2);
    mpz_init(partial3);

    //mpz_sub(partial2, P->X, P->X);
    //mpz_mul(P->Y, m, partial2);

    //mpz_sub(partial2, res->X, P->X);
    mpz_sub(partial2,  P->X, result.X);

    mpz_mul(partial3, m, partial2);

    mpz_sub(result.Y, partial3, P->Y);


    //gmp_printf("partial2 = %Zd\n", partial2);
    //mpz_mul(res->Y, m, partial2);


    mpz_clear(partial2);
    mpz_clear(partial3);
    mpz_clear(m);

    //mpz_mod(P->Y, P->Y, pd.n);
    mpz_mod(result.Y, result.Y, pd.n);

    /*mp_bitcnt_t size = mpz_size(pd.n);
    mpz_realloc2(PandQ->Y, size);*/

    //mpz_set_ui(P->Z, 1);
    //mpz_set(P->X, res->X);
    //mpz_set(P->Y, res->Y);
    mpz_set_ui(result.Z, 1);

    //free(P);
    //sleep(1);
    return result;
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


int checkIfCurve(struct ECpoint P, struct weirstrassEC EC, struct problemData pd)
{
    //y^2 = x^3 + ax + b
    //verifico l'uguaglianza

    mpz_t squarey, cubex, tot, sum;
    mpz_init(squarey);
    mpz_init(cubex);
    mpz_init(tot);
    mpz_init(sum);

    mpz_pow_ui(squarey, P.Y, 2);
    mpz_pow_ui(cubex, P.X, 3);
    //mpz_mul(ax, P.X, EC.a);

    mpz_add(sum, cubex, P.X);
    mpz_add(sum, sum, EC.b);

    mpz_sub(tot, squarey, sum);

    mpz_mod(tot, tot, pd.n);

    mpz_clear(squarey);
    mpz_clear(cubex);
    mpz_clear(sum);

    if(mpz_cmp_ui(tot, 0) == 0)
    {
        //printf("verified\n");
        mpz_clear(tot);
        return 1;
    }
    else
    {
        printf("failed\n");
        sleep(4);
        mpz_clear(tot);

        //gmp_printf("Y^2 = %Zd\n", squarey);
        //gmp_printf("sum = %Zd\n", sum);

        return 0;
    }


}