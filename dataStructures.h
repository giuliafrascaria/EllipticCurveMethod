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
    mpz_t Y;                //used to indicate the point to infinity
    mpz_t Z;
};

struct problemData
{
    mpz_t n;                //large composite number to be factored

    mpz_t stageOneB;        //must be even
    mpz_t stageTwoB;        //100B1
    mpz_t D;                //total memory
    unsigned long Dint;
    long iterations;
};

struct ellipticCurve
{
    //not weirstrass form, this is the invertionless algorithm
    mpz_t sigma;
    mpz_t u;
    mpz_t v;
    mpz_t C;

};

struct weirstrassEC
{
    //weirstrass form, this is the invertionless algorithm
    mpz_t a;
    mpz_t b;

};

struct nonInvertibleD
{
    int flag;
    mpz_t d;
};

struct JsElem
{
    unsigned long index;
    struct JsElem * next;
};

void addh();
void doubleh();
struct ECpoint * add(struct ECpoint *P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res);
struct ECpoint negate(struct ECpoint *P);
struct ECpoint *doubleec(struct ECpoint *P, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD *d, struct ECpoint * res);
struct ECpoint ECmultiplyTraditional(struct ECpoint * Q, mpz_t p, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res);
struct ECpoint ECmultiplyMontgomery(struct ECpoint Q, mpz_t p);
void add2(struct ECpoint * P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint *res);
void sub2(struct ECpoint *P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res);
void doubleec2(struct ECpoint * P, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD *d, struct ECpoint *res);
struct ECpoint doubleAndAdd(struct ECpoint * P, mpz_t p,  struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res);
int checkIfCurve(struct ECpoint P, struct weirstrassEC EC, struct problemData pd);
struct ECpoint doubleAndAdd2(struct ECpoint * P, mpz_t p, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d);
struct ECpoint add3(struct ECpoint * P, struct ECpoint *Q, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d);
struct ECpoint ECmultiplyTraditional2(struct ECpoint * Q, mpz_t *p, struct weirstrassEC EC, struct problemData pd, struct nonInvertibleD * d, struct ECpoint * res);


#endif //ECM_DATASTRUCTURES_H
