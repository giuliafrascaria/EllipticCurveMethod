//
// Created by root on 05/03/17.
//


#include <stdlib.h>
#include <stdio.h>
#include "dataStructures.h"
#include "gc/include/gc/gc.h"
#include <unistd.h>
#include <math.h>

int classicalECM(struct problemData pd, mpz_t * factor, gmp_randstate_t state);
int traditionalStageOne(struct weirstrassEC EC, struct ECpoint Q, struct problemData pd, mpz_t  *factor, struct ECpoint * returnQ);
int stageTwo(struct weirstrassEC EC, struct ECpoint Q, struct problemData pd);
void getBBest(int logp, mpz_t * B);


//void randomEC(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd);
//void stageOne(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd);
void randomECtraditional(struct weirstrassEC * EC, struct ECpoint  *Q, struct problemData pd, gmp_randstate_t state);


int main(int argc, char ** argv)
{
    //GC_INIT();
    if(argc != 2)
    {
        perror("missing number");
        return EXIT_FAILURE;
    }
    //convert input into internal representation
        //declarations and init
    struct problemData pd;
    mpz_init(pd.n);
        //conversion
    printf("%s\n", argv[1]);
    if(mpz_set_str(pd.n, (const char *) argv[1], 10) != 0)
    {
        perror("invalid input type");
        return EXIT_FAILURE;
    }
    //choose criteria
        //init
    mpz_init(pd.stageOneB);
    mpz_init(pd.stageTwoB);
    mpz_init(pd.D);
        //assign values
    mpz_set_ui(pd.stageOneB, 100);


    mpz_sqrt(pd.stageTwoB, pd.n);
    gmp_printf("square root of n = %Zd\n", pd.stageTwoB);


    //initialize random state
    gmp_randstate_t state;
    gmp_randinit_mt(state);

    //sleep(1);

    mpz_t factor;
    mpz_init(factor);
    int success = 0;
    int primeLog = 9;   //starting to find primes with 3 digits (this is natural log)
    double i, ww, w;
    long iterations;

    //get BBEst

   /* int maxlog = primeLog+10;

    pid_t pid;
    pid = fork();


    if(pid == 0)
    {
        primeLog = primeLog + 10;
        maxlog = maxlog + 10;
    }
*/


    //while((mpz_cmp(pd.stageOneB, pd.stageTwoB) < 0 && primeLog < maxlog) || pid == 0)
    while((mpz_cmp(pd.stageOneB, pd.stageTwoB) < 0))
    {
        getBBest(primeLog, &pd.stageOneB);

        w = primeLog/log(mpz_get_ui(pd.stageOneB));

        ww = pow(w, w);
        iterations = lround(ww);
        //printf("process %d looking for factors with approximately %d digits, iterating for %ld times\n", pid, primeLog/3, iterations);



        gmp_printf("BBest = %Zd\n\n", pd.stageOneB);
        for(i = 0; i < iterations; i++)       //repeat w^w times
        {
            success = classicalECM(pd, &factor, state);
            if(success)
            {
                exit(EXIT_SUCCESS);
            }
        }

        primeLog+=3;
        //get new BBest

        //printf("failure, try again and eventually increment B\n\n");
        //mpz_mul_ui(pd.stageOneB, pd.stageOneB, 100);
        //sleep(1);
    }
    //printf("the size of B1 exceeded max B, process %d quitting\n", pid);


    //loop through the next steps when there is a significant chance that there are no factors with logB digits

    exit(EXIT_FAILURE);
}

void getBBest(int logp, mpz_t * B)
{
    double w0; //sqrt of(2logp/loglogp)
    double den = log(logp);
    w0 = sqrt((2*logp)/den);

    double Bbest = pow(M_E, logp/w0);

    long integerBBest = lround(Bbest);

    mpz_set_ui(*B, (unsigned long) integerBBest);

}

int classicalECM(struct problemData pd, mpz_t *factor, gmp_randstate_t state) {
    struct ECpoint Q;
    struct weirstrassEC EC;
    int success = 0;

    mpz_init(Q.X);
    mpz_init(Q.Y);
    mpz_init(Q.Z);

    mpz_set_ui(Q.Z, 1);

    mpz_init(EC.a);
    mpz_init(EC.b);

    randomECtraditional(&EC, &Q, pd, state);
    /*printf("random init terminated\n");

    gmp_printf("%Zd\n", Q.X);

    gmp_printf("%Zd\n", Q.Y);

    gmp_printf("%Zd\n", Q.Z);

    gmp_printf("%Zd\n", EC.a);*/

    //sleep(1);


    //g calculus
    mpz_t g;
    mpz_init(g);

    mpz_t cubea, term1, squareb, term2, gcdterm;
    mpz_init(cubea);
    mpz_init(term1);
    mpz_init(squareb);
    mpz_init(term2);
    mpz_init(gcdterm);

    mpz_pow_ui(cubea, EC.a, 3);
    mpz_mul_ui(term1, cubea, 4);

    mpz_pow_ui(squareb, EC.b, 2);
    mpz_mul_ui(term2, squareb, 27);

    mpz_add(gcdterm, term1, term2);

    mpz_gcd(g, gcdterm, pd.n);


    if (mpz_cmp(g, pd.n) == 0)
    {
        //find new curve

        printf("bad EC, restarting\n");
        classicalECM(pd, factor, state);
    }
    else if (mpz_cmp_ui(g, 1) > 0)
    {
        //return gcdterm as factor
        printf("found factor with random choice\n");
        gmp_printf("%Zd\n\n", g);
        return 1;
    }
    else
    {
        //prime power multipliers
        //printf("starting step 1 in seconds\n");
        //sleep(2);
        struct ECpoint result;
        mpz_init(result.X);
        mpz_init(result.Y);
        mpz_init(result.Z);

        success = traditionalStageOne(EC, Q, pd, factor, &result);

        if(success == 1)
        {
            return success;
        }
        else
        {
            printf("trying stage 2\n");
            success = stageTwo(EC, result, pd);
            if(success)
            {
                printf("successful stage two\n");
            }
        }

        mpz_clear(result.X);
        mpz_clear(result.Y);
        mpz_clear(result.Z);

    }
    mpz_clear(Q.X);
    mpz_clear(Q.Y);
    mpz_clear(Q.Z);

    mpz_clear(EC.b);
    mpz_clear(EC.a);

    mpz_clear(cubea);
    mpz_clear(term1);
    mpz_clear(squareb);
    mpz_clear(term2);
    mpz_clear(gcdterm);
    mpz_clear(g);


    return success;
}

/*void invertionlessECM(struct problemData pd)
{
    //choose random EC
    //declaration and init
    struct ECpoint startP;
    struct ellipticCurve EC;
    mpz_init(startP.X);
    mpz_init(startP.Z);
    mpz_init(EC.sigma);
    mpz_init(EC.C);
    mpz_init(EC.u);
    mpz_init(EC.v);
    //generation
    randomEC(EC, startP, pd);
    //stage one
    stageOne(EC, startP, pd);
    //stage two
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

    //u
    mpz_t squaresigma, nomodval;
    mpz_init(squaresigma);
    mpz_init(nomodval);
    mpz_pow_ui(squaresigma, EC.sigma, 2);
    unsigned long int five = 5;
    mpz_sub_ui(nomodval, squaresigma, five);
    mpz_mod(EC.u, nomodval, pd.n);

    //v
    mpz_t foursigma;
    mpz_init(foursigma);
    unsigned long int four = 4;
    mpz_mul_ui(foursigma, EC.sigma, four);
    mpz_mod(EC.v, foursigma, pd.n);

    //C


    //coordinates for initial point Q (Montgomery representation)
    mpz_powm_ui(Q.X, EC.u, 3, pd.n);
    mpz_powm_ui(Q.Z, EC.u, 3, pd.n);
}

void stageOne(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd)
{
    //for cycle between all primes <= B1
    //      pi(B) is the number of primes less than B

    //funzioni utili
    //void mpz_nextprime (mpz t rop, const mpz t op)
    //void mpz_gcd (mpz t rop, const mpz t op1, const mpz t op2)

    mpz_t primen;
    mpz_init(primen);
    mpz_set_str(primen, "1", 10);

    *//*
     * int mpz_cmp_ui (const mpz t op1, unsigned long int op2) [Macro]
Compare op1 and op2. Return a positive value if op1 > op2, zero if op1 = op2, or a negative
value if op1 < op2.
     * *//*

    while(mpz_cmp(primen, pd.stageOneB) <= 0)
    {
        //perform operations
        //find largest integer a such that pi^a <= B1
        unsigned long exp = 1;
        int flag = 1;
        mpz_t power;
        mpz_init(power);
        while(flag)
        {
            mpz_pow_ui(power, primen, exp);
            if(mpz_cmp(power, pd.stageOneB) <= 0)
            {
                exp += 1;
                //looking for a bigger exponent
            }
            else
            {
                flag = 0;
                //Q = [pi^a]Q
                struct ECpoint nQ = ECmultiplyMontgomery(Q, power);
            }
        }

        //find next prime
        mpz_nextprime(primen, primen);

    }
    //gcd
}*/

int traditionalStageOne(struct weirstrassEC EC, struct ECpoint Q, struct problemData pd, mpz_t * factor, struct ECpoint * returnQ)
{
    //for cycle between all primes <= B1
    //      pi(B) is the number of primes less than B

    //funzioni utili
    //void mpz_nextprime (mpz t rop, const mpz t op)
    //void mpz_gcd (mpz t rop, const mpz t op1, const mpz t op2)

    mpz_t primen;
    mpz_init(primen);
    mpz_set_str(primen, "2", 10);

    /*
     * int mpz_cmp_ui (const mpz t op1, unsigned long int op2) [Macro]
Compare op1 and op2. Return a positive value if op1 > op2, zero if op1 = op2, or a negative
value if op1 < op2.
     * */


    struct ECpoint P;

    mpz_init(P.X);
    mpz_init(P.Y);
    mpz_init(P.Z);

    mpz_set(P.X, Q.X);
    mpz_set(P.Y, Q.Y);
    mpz_set(P.Z, Q.Z);
//
//
//    gmp_printf("x = %Zd , y= %Zd , z= %Zd \n", P.X, P.Y, P.Z);
    //sleep(3);
    while(mpz_cmp(primen, pd.stageOneB) <= 0)
    {
        //perform operations
        //find largest integer a such that pi^a <= B1
        unsigned long exp = 1;
        int flag = 1;
        mpz_t power;
        mpz_init(power);
        //printf("i'm inside the loop 1\n");
        while(flag)
        {
            mpz_pow_ui(power, primen, exp);
            //printf("i'm inside the loop 2\n");
            if(mpz_cmp(power, pd.stageOneB) <= 0)
            {
                exp += 1;
                //looking for a bigger exponent
            }
            else
            {
                flag = 0;
                //printf("the exponent is %lu\n", exp);
            }
        }
        mpz_clear(power);


        unsigned int i;

        struct nonInvertibleD d;
        mpz_init(d.d);
        d.flag = 0;

       // printf("starting prime power multipliers\n\n");
        for(i = 1; i <= exp; i++)
        {
            //Q = [pi]Q

            //gmp_printf("moltiplico per il primo %Zd\n", primen);
            P = ECmultiplyTraditional(&P, primen, EC, pd, &d, &P);

            //gmp_printf("px %Zd, py %Zd, pz %Zd\n", P.X, P.Y, P.Z);
            //sleep(1);

            if(d.flag)
            {
                mpz_t newFactor;
                mpz_init(newFactor);
                mpz_gcd(newFactor, pd.n, d.d);
                gmp_printf("found factor %Zd\n", newFactor);

                mpz_set(*factor, newFactor);

                return 1;
            }

            //the cycle stops if I find a non invertible denominator in the addition slope
        }

        mpz_clear(d.d);
        //printf("on with another prime\n");


        //find next prime
        mpz_nextprime(primen, primen);
        //gmp_printf("%Zd\n\n", primen);

    }
    //printf("cycle failed\n");
    mpz_set(returnQ->X, P.X);
    mpz_set(returnQ->Y, P.Y);
    mpz_set(returnQ->Z, P.Z);

    mpz_clear(primen);

    mpz_clear(P.X);
    mpz_clear(P.Y);
    mpz_clear(P.Z);

    return 0;

}

int stageTwo(struct weirstrassEC EC, struct ECpoint Q, struct problemData pd)
{
    //succede se per nessun fattore primo p di n la curva E( ZZp) ha ordine B-smooth.
    // si dovrebbero individuare fattori primi p di n per cui l’ordine della stessa
    //curva E su ZZp sia della forma
    //#E( ZZp) = q · m,

    //la fase 2 ha successo se esiste un p divisore di n per cui ordine della curva ellittica è Bsmooth = m*q
    //q è un primo in un intervallo B1 B2 con B2 uguale a circa 100 B1

    //precalcolo: lista dei primi tra B1 e B2, differenze tra adiacenti

    //dalla fase uno so che KP = Q
    //Ri = di*Q
    //Q1 = q1*Q
    //Q2 = Q1 + R1 = (q1 + d1)Q = q2*Q ecc

    //se tutti questi sono al finito ha fallito anche la fase 2
    //altrimenti esiste qi*kP che va a infinito -> ordine di E Zp divide qi*k e il denominatore non invertibile è un fattore non banale di n

    mpz_t q;
    mpz_t nextq;
    mpz_t delta;
    mpz_init(q);
    mpz_init(nextq);
    mpz_init(delta);

    mpz_t B2;
    mpz_init(B2);
    mpz_mul_ui(B2, pd.stageOneB, 10);

    struct ECpoint Qi, term2;

    struct ECpoint * previousCalc;

    previousCalc = &Qi;

    mpz_init(Qi.X);
    mpz_init(Qi.Y);
    mpz_init(Qi.Z);

    mpz_set(Qi.X, Q.X);
    mpz_set(Qi.Y, Q.Y);
    mpz_set(Qi.Z, Q.Z);

    mpz_init(term2.X);
    mpz_init(term2.Y);
    mpz_init(term2.Z);

    mpz_set(term2.X, Q.X);
    mpz_set(term2.Y, Q.Y);
    mpz_set(term2.Z, Q.Z);

    struct nonInvertibleD d;
    mpz_init(d.d);
    d.flag = 0;


    //gmp_printf("qx %Zd, qy %Zd, qz %Zd\n", Q.X, Q.Y, Q.Z);
    //sleep(1);

    mpz_nextprime(nextq, pd.stageOneB);
    mpz_nextprime(q, pd.stageOneB);

    *previousCalc = ECmultiplyTraditional(previousCalc, q, EC, pd, &d, previousCalc);
    //P = ECmultiplyTraditional(&P, primen, EC, pd, &d, &P);

    //gmp_printf("x %Zd, y %Zd, z %Zd\n", Qi.X, Qi.Y, Qi.Z);
    //sleep(1);



    while(mpz_cmp(q, B2) <= 0)
    {

        mpz_nextprime(nextq, q);
        mpz_sub(delta, nextq, q);

        term2 = ECmultiplyTraditional(&term2, delta, EC, pd, &d, &term2);

        previousCalc = add(previousCalc, &term2, EC, pd, &d, previousCalc);
        //gmp_printf("x %Zd, y %Zd, z %Zd\n", previousCalc->X, previousCalc->Y, previousCalc->Z);
        //sleep(1);
        //devo tenere in memoria la precedente moltiplicazione e fare quella nuova, moltiplicando Q per deltai = qi+1 - qi

        if(d.flag == 1)
        {
            //denominatore non invertibile, found factor
            mpz_t newFactor;
            mpz_init(newFactor);
            mpz_gcd(newFactor, pd.n, d.d);
            gmp_printf("found factor during stage two %Zd\n", newFactor);

            return 1;

        }
        else
        {
            mpz_swap(q, nextq);

            /*mpz_nextprime(q, q);

            mpz_set(Qi.X, Q.X);
            mpz_set(Qi.Y, Q.Y);
            mpz_set(Qi.Z, Q.Z);*/

            //gmp_printf("\n\ntrying prime %Zd\n\n", nextq);
        }
    }

    mpz_clear(Qi.X);
    mpz_clear(Qi.Y);
    mpz_clear(Qi.Z);

    mpz_clear(B2);
    mpz_clear(d.d);

    mpz_clear(q);
    mpz_clear(nextq);

    return 0;   //failure

}

void randomECtraditional(struct weirstrassEC * EC, struct ECpoint * Q, struct problemData pd, gmp_randstate_t state)
{
    //random sigma, derive u v anc C from this

    //void gmp_randseed (gmp randstate t state, const mpz t seed)

    //generate random x, y, a in [0, n-1]

    mpz_urandomm(Q->X, state, pd.n);
    //gmp_printf("%Zd\n", Q->X);
    mpz_urandomm(Q->Y, state, pd.n);
    //gmp_printf("%Zd\n", Q->Y);
    mpz_urandomm(EC->a, state, pd.n);
    //gmp_printf("%Zd\n", EC->a);
    //the curve is determined by these values

    //sleep(2);
    mpz_t squareY, cubeX, aX;
    mpz_init(squareY);
    mpz_init(cubeX);
    mpz_init(aX);
    mpz_pow_ui(squareY, Q->Y, 2);
    mpz_pow_ui(cubeX, Q->X, 3);
    mpz_mul(aX, Q->X, EC->a);


    mpz_t temp;
    mpz_init(temp);
    mpz_sub(temp, squareY, cubeX);
    mpz_sub(EC->b, temp, aX);
    mpz_mod(EC->b, EC->b, pd.n);

    mp_bitcnt_t size = mpz_size(pd.n);
    mpz_realloc2(EC->b, size);

    mpz_clear(squareY);
    mpz_clear(cubeX);
    mpz_clear(aX);
    mpz_clear(temp);
    //gmp_printf("%Zd\n", EC->b);

}