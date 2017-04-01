//
// Created by root on 05/03/17.
//


#include <stdlib.h>
#include <stdio.h>
#include "dataStructures.h"
#include <unistd.h>
#include <math.h>
#include <pthread.h>

int classicalECM(struct problemData pd, mpz_t * factor, gmp_randstate_t state, int k, int digits);
int traditionalStageOne(struct weirstrassEC EC, struct ECpoint Q, struct problemData pd, mpz_t  *factor, struct ECpoint * returnQ);
int stageTwo(struct weirstrassEC EC, struct ECpoint Q, struct problemData pd);
void getBBest(int logp, mpz_t * B);
void * loop(void * k);
void optimizedRandomEC(struct weirstrassEC * EC, struct ECpoint * Q, struct problemData pd, gmp_randstate_t state);
//void randomEC(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd);
//void stageOne(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd);
void randomECtraditional(struct weirstrassEC * EC, struct ECpoint  *Q, struct problemData pd, gmp_randstate_t state);
int efficientStageTwo(struct weirstrassEC EC, struct ECpoint Q, struct problemData pd);



struct problemData pd;
int maxdigits = 9;
pthread_mutex_t stage2mtx[8];


int main(int argc, char ** argv)
{


    if(argc != 3)
    {
        perror("missing numbers");
        return EXIT_FAILURE;
    }


    //convert input into internal representation
        //declarations and init

    mpz_init(pd.n);
        //conversion
    printf("%s\n", argv[1]);
    if(mpz_set_str(pd.n, (const char *) argv[1], 10) != 0)
    {
        perror("invalid input type");
        return EXIT_FAILURE;
    }

    mpz_t f1, f2;
    mpz_init(f1);
    mpz_init(f2);
    //mpz_nextprime(f1, )
    //choose criteria
        //init

    for(int j = 0; j < 8; j++)
    {
        pthread_mutex_init(&stage2mtx[j], NULL);
    }

    mpz_init(pd.stageOneB);
    mpz_init(pd.stageTwoB);
    mpz_init(pd.D);
        //assign values
    //mpz_set_ui(pd.stageOneB, 100);


    if(mpz_set_str(pd.stageOneB, (const char *) argv[2], 10) != 0)
    {
        perror("invalid input type");
        return EXIT_FAILURE;
    }

    /*if(mpz_set_str(pd.stageTwoB, (const char *) argv[3], 10) != 0)
    {
        perror("invalid input type");
        return EXIT_FAILURE;
    }*/

    mpz_mul_ui(pd.stageTwoB, pd.stageOneB, 100);


   // mpz_sqrt(pd.stageTwoB, pd.n);
   // gmp_printf("square root of n = %Zd\n", pd.stageTwoB);


    //get BBEst

    /*int maxlog = primeLog+10;

    pid_t pid;
    pid = fork();

    if(pid == 0)
    {
        primeLog = primeLog + 10;
        maxlog = maxlog + 10;
    }
*/

    pthread_t t[7];
    for(int k = 0; k < 7; k++)
    {

        int c =k;
        pthread_create(&t[k], NULL, loop, (void *) &c);
        usleep(100);
    }

    int seven = 7;
    loop(&seven);
}

void * loop(void * k)
{
    int * index = (int *) k;
    int cont = *index;
    //printf("new thread with mtx %d\n", cont);

    pthread_t tid = pthread_self();
    //initialize random state
    gmp_randstate_t state;
    gmp_randinit_mt(state);

    gmp_randseed_ui(state, tid);

    //sleep(1);

    mpz_t factor;
    mpz_init(factor);
    int success = 0;
    int primeLog = 9;   //starting to find primes with 3 digits (this is natural log)
    double ww, w;
    long i;
    long iterations;

    //while((mpz_cmp(pd.stageOneB, pd.stageTwoB) < 0 && primeLog < maxlog) || pid == 0)
    while((mpz_cmp(pd.stageOneB, pd.stageTwoB) < 0))
    {
        //getBBest(primeLog, &pd.stageOneB);        //alpertron
      //  mpz_set_ui(pd.stageOneB, 100);


        w = primeLog/log(mpz_get_ui(pd.stageOneB));


        ww = pow(w, w);
        iterations = lround(ww);
        //iterations = 75;
        printf("process iterating for %ld times\n", iterations);


        //gmp_printf("BBest = %Zd\n\n", pd.stageOneB);
        for(i = 0; i < iterations; i++)       //repeat w^w times
        {
            success = classicalECM(pd, &factor, state, cont, primeLog);
            if(success == 1)
            {
                exit(EXIT_SUCCESS);
            }
            //printf("%ld\n", i);
        }

        primeLog+=3;
        //get new BBest

        printf("failure, incrementing B\n\n");
        //exit(EXIT_FAILURE);
        mpz_mul_ui(pd.stageOneB, pd.stageOneB, 10);
        //sleep(1);
    }
    printf("found no factor\n");



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

int classicalECM(struct problemData pd, mpz_t *factor, gmp_randstate_t state, int k, int digits)
{


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
    //optimizedRandomEC(&EC, &Q, pd, state);

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
        classicalECM(pd, factor, state, k, digits);
    }
    else if (mpz_cmp_ui(g, 1) > 0)
    {
        //return gcdterm as factor
        printf("-------------------------------------\n\n\tfound factor with random choice\n");
        gmp_printf("\t\t%Zd\n-------------------------------\n\n", g);
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
            printf("successful stage 1\n");
            return success;
        }
        else if(success == 0)
        {
            printf("should try stage 2\n");
            //sleep(1);

            if(digits > maxdigits || 1==1)
            {

//                if(pthread_mutex_trylock(&(stage2mtx[k])) == 0)
                if(pthread_mutex_lock(&(stage2mtx[0])) == 0)
                {
                    //printf("thread %ld locked %d mtxfor %d digits \n", pthread_self(), k, digits);

                    maxdigits = digits;
                    //success = stageTwo(EC, result, pd);
                    success = efficientStageTwo(EC, result, pd);
                    if(success)
                    {
                        printf("successful stage two\n");
                    }
                    printf("thread %ld in stage two for %d digits\n", pthread_self(),digits/3);
                    //sleep(1);
//                    pthread_mutex_unlock(&(stage2mtx[k]));
                    pthread_mutex_unlock(&(stage2mtx[0]));
                    //printf("thread %ld leaving stage two for ln = %d\n", pthread_self(),digits);
                    return success;
                }
            }
        }
        else
        {
            printf("found point to infinity, trying another curve\n");
            //sleep(1);
            return 2;
        }
        /*else
        {
            printf("should try stage 2\n");
            //sleep(1);
            success = stageTwo(EC, result, pd);
            if(success)
            {
                printf("successful stage two\n");
            }
            return success;
        }*/

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

void optimizedRandomEC(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd)
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

 /*   mpz_init(res.X);
    mpz_init(res.Y);
    mpz_init(res.Z);*/

    mpz_set(P.X, Q.X);
    mpz_set(P.Y, Q.Y);
    mpz_set(P.Z, Q.Z);
//
//
//   gmp_printf("x = %Zd , y= %Zd , z= %Zd \n", P.X, P.Y, P.Z);
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

            /*gmp_printf("moltiplico per il primo %Zd\n", primen);
            gmp_printf("prep\n\tx = %Zd , y= %Zd , z= %Zd \n", P.X, P.Y, P.Z);
            gmp_printf("pre\n\trx %Zd, ry %Zd, rz %Zd\n", returnQ->X, returnQ->Y, returnQ->Z);*/

            /*gmp_printf("b = %Zd\n\n", EC.b);
            gmp_printf("qx = %Zd\nqy = %Zd\nqz = %Zd\n", P.X, P.Y, P.Z);*/

            *returnQ = ECmultiplyTraditional(&P, primen, EC, pd, &d, returnQ);

            //doubleAndAdd(&P, primen, EC, pd, &d, returnQ);
            //*returnQ = doubleAndAdd2(&P, primen, EC, pd, &d);

            //gmp_printf("post\n\trx %Zd, ry %Zd, rz %Zd\n", returnQ->X, returnQ->Y, returnQ->Z);
            //gmp_printf("postp\n\tx = %Zd , y= %Zd , z= %Zd \n", P.X, P.Y, P.Z);
            if(mpz_cmp_ui(returnQ->X, 0) == 0)
            {
                if(mpz_cmp_ui(returnQ->Y, 1) == 0)
                {
                    if(mpz_cmp_ui(returnQ->Z, 0) == 0)
                    {
                        printf("point to infinity for EZn, bad luck check 3\n");
                        return 2;
                    }
                }
            }
            //sleep(1);

            if(d.flag)
            {
                mpz_t newFactor;
                mpz_init(newFactor);
                mpz_gcd(newFactor, pd.n, d.d);

                if(mpz_cmp(newFactor, pd.n) == 0)
                {
                    printf("bad luck\n");
                    mpz_clear(newFactor);
                    return 0;
                }
                gmp_printf("-------------------------------------\n\n\tfound factor %Zd\n\n-------------------------------------\n\n", newFactor);

                mpz_set(*factor, newFactor);

                return 1;
            }

            /*gmp_printf("qqqqqpost\n\trx %Zd, ry %Zd, rz %Zd\n", returnQ->X, returnQ->Y, returnQ->Z);
            gmp_printf("qqqqpostp\n\tx = %Zd , y= %Zd , z= %Zd \n", P.X, P.Y, P.Z);
*/
            mpz_set(P.Z, returnQ->Z);
            mpz_set(P.Y, returnQ->Y);
            mpz_set(P.X, returnQ->X);

            /*gmp_printf("post\n\trx %Zd, ry %Zd, rz %Zd\n", returnQ->X, returnQ->Y, returnQ->Z);
            gmp_printf("postp\n\tx = %Zd , y= %Zd , z= %Zd \n", P.X, P.Y, P.Z);
*/
            /*mpz_set_ui(returnQ->X, 0);
            mpz_set_ui(returnQ->Y, 0);
            mpz_set_ui(returnQ->Z, 0);
*/

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

    /*mpz_clear(res.X);
    mpz_clear(res.Y);
    mpz_clear(res.Z);
*/


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
    mpz_mul_ui(B2, pd.stageOneB, 100);

    struct ECpoint Qi, term2, res, res2;

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

    mpz_init(res.X);
    mpz_init(res.Y);
    mpz_init(res.Z);

    mpz_init(res2.X);
    mpz_init(res2.Y);
    mpz_init(res2.Z);

    struct nonInvertibleD d;
    mpz_init(d.d);
    d.flag = 0;


    //gmp_printf("stage 2 con qx %Zd, qy %Zd, qz %Zd\n", Q.X, Q.Y, Q.Z);
    //sleep(1);

    mpz_nextprime(nextq, pd.stageOneB);
    mpz_nextprime(q, pd.stageOneB);

    //*previousCalc = ECmultiplyTraditional(previousCalc, q, EC, pd, &d, previousCalc);
    //sostituire con una cosa tipo
    res = ECmultiplyTraditional(previousCalc, q, EC, pd, &d, &res);


    //P = ECmultiplyTraditional(&P, primen, EC, pd, &d, &P);

    //gmp_printf("x %Zd, y %Zd, z %Zd\n", Qi.X, Qi.Y, Qi.Z);
    //sleep(1);

    while(mpz_cmp(q, B2) <= 0)
    {

        mpz_nextprime(nextq, q);
        mpz_sub(delta, nextq, q);

        res2 = ECmultiplyTraditional(&term2, delta, EC, pd, &d, &res2);

        //previousCalc = add(previousCalc, &term2, EC, pd, &d, previousCalc);
        //printf("\n\n\n now step 2\n\n\n");
        //sleep(3);
        add2(&res, &res2, EC, pd, &d, previousCalc);
        //gmp_printf("x %Zd, y %Zd, z %Zd\n", previousCalc->X, previousCalc->Y, previousCalc->Z);
        //sleep(1);
        //devo tenere in memoria la precedente moltiplicazione e fare quella nuova, moltiplicando Q per deltai = qi+1 - qi

        if(d.flag == 1)
        {
            //denominatore non invertibile, found factor
            mpz_t newFactor;
            mpz_init(newFactor);
            mpz_gcd(newFactor, pd.n, d.d);
            if(mpz_cmp(newFactor, pd.n) == 0)
            {
                mpz_clear(newFactor);
                printf("bad luck\n");
                return 0;
            }
            gmp_printf("------------------------------------------------------\n\n\tfound factor during stage two %Zd\n\n-----------------------------------------------------\n", newFactor);

            return 1;

        }
        else
        {
            mpz_swap(q, nextq);

            mpz_set(res.X, previousCalc->X);
            mpz_set(res.Y, previousCalc->Y);
            mpz_set(res.Z, previousCalc->Z);

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

    mpz_set_ui(EC->a, 1);
    //gmp_printf("%Zd\n", Q->Y);
    //mpz_urandomm(EC->a, state, pd.n);       //posso sostituire per avere un A = 1
    //gmp_printf("%Zd\n", EC->a);
    //the curve is determined by these values

    //sleep(2);
    mpz_t squareY, cubeX;
    mpz_init(squareY);
    mpz_init(cubeX);
    //mpz_init(aX);
    mpz_pow_ui(squareY, Q->Y, 2);
    mpz_pow_ui(cubeX, Q->X, 3);
    //mpz_mul(aX, Q->X, EC->a);


    mpz_t temp;
    mpz_init(temp);
    mpz_sub(temp, squareY, cubeX);
    mpz_sub(EC->b, temp, Q->X);
    //mpz_mod(EC->b, EC->b, pd.n);

    //mp_bitcnt_t size = mpz_size(pd.n);
    //mpz_realloc2(EC->b, size);

    mpz_clear(squareY);
    mpz_clear(cubeX);
    //mpz_clear(aX);
    mpz_clear(temp);
    //gmp_printf("b = %Zd\n", EC->b);
    //gmp_printf("a = %Zd\n", EC->a);

    //checkIfCurve(*Q, *EC, pd);
    //printf("peh\n");

}

void optimizedRandomEC(struct weirstrassEC * EC, struct ECpoint * Q, struct problemData pd, gmp_randstate_t state)
{
    //random sigma, derive u v anc C from this

    //initialize random state


    //generate random sigma in [6, n-1]
    unsigned long int six = 6;
    mpz_t maxsigma, sigma;
    mpz_init(maxsigma);
    mpz_init(sigma);
    ;
    mpz_sub_ui(maxsigma, pd.n, six);
    mpz_urandomm(sigma, state, maxsigma);
    mpz_add_ui(sigma, sigma, six);

    //the curve is determined by these values

    //u
    mpz_t squaresigma, nomodval, u;
    mpz_init(squaresigma);
    mpz_init(nomodval);
    mpz_init(u);

    mpz_pow_ui(squaresigma, sigma, 2);
    unsigned long int five = 5;
    mpz_sub_ui(nomodval, squaresigma, five);
    mpz_mod(u, nomodval, pd.n);

    //v
    mpz_t foursigma, v;
    mpz_init(foursigma);
    mpz_init(v);

    unsigned long int four = 4;
    mpz_mul_ui(foursigma, sigma, four);
    mpz_mod(v, foursigma, pd.n);

    //C


    //coordinates for initial point Q (Montgomery representation)
    mpz_powm_ui(Q->X, u, 3, pd.n);
    //y = s^2-1)*s^2-25)*s^4-25

    mpz_t term1, term2, sigmafour, term3;
    mpz_init(term1);
    mpz_init(term2);
    mpz_init(term3);
    mpz_init(sigmafour);

    mpz_sub_ui(term1, squaresigma, 1);
    mpz_sub_ui(term2, squaresigma, 25);
    mpz_pow_ui(sigmafour, squaresigma, 2);
    mpz_sub_ui(term3, sigmafour, 25);

    mpz_mul(Q->Y, term1, term2);
    mpz_mul(Q->Y, Q->Y, term3);
    mpz_mod(Q->Y, Q->Y, pd.n);

    mpz_powm_ui(Q->Z, v, 3, pd.n);

    //coordinates for elliptic curve
    mpz_t u3, partial1, partial2, partial3;
    mpz_init(u3);
    mpz_init(partial1);
    mpz_init(partial2);
    mpz_init(partial3);

    mpz_div(EC->b, u, Q->Z);
    mpz_mod(EC->b, EC->b, pd.n);

    mpz_sub(partial1, v, u);
    mpz_pow_ui(partial1, partial1, 3);

    mpz_mul_ui(u3, u, 3);
    mpz_add(partial2, u3, v);

    mpz_mul(partial3, Q->X, v);
    mpz_mul_ui(partial3, partial3, 4);

    mpz_mul(partial1, partial1, partial2);
    mpz_div(partial1, partial1, partial3);
    mpz_sub_ui(EC->a, partial1, 2);
    mpz_mod(EC->a, EC->a, pd.n);

    mpz_clear(u);
    mpz_clear(v);
    mpz_clear(partial1);
    mpz_clear(partial2);
    mpz_clear(partial3);
    mpz_clear(term1);
    mpz_clear(term2);
    mpz_clear(term3);
    mpz_clear(u3);
    mpz_clear(sigmafour);
    mpz_clear(sigma);
    mpz_clear(squaresigma);
    mpz_clear(nomodval);
    mpz_clear(foursigma);
    mpz_clear(maxsigma);

}

int efficientStageTwo(struct weirstrassEC EC, struct ECpoint Q, struct problemData pd)
{

    /*gmp_printf("Qx = %Zd\n", Q.X);
    gmp_printf("Qy = %Zd\n", Q.Y);
    gmp_printf("Qz = %Zd\n", Q.Z);*/


    printf("precomputation for stage two\n");

    struct JsElem head;

    head.next = NULL;

    mpz_t B2, D, Mmin, Mmax, Dhalf, t1,t2;
    double dD, dMmin, dMmax, dDhalf;

    mpz_init(B2);
    mpz_init(D);
    mpz_init(Mmax);
    mpz_init(Mmin);
    mpz_init(Dhalf);
    mpz_init(t1);
    mpz_init(t2);

    mpz_mul_ui(B2, pd.stageOneB, 100);
    mpz_sqrt(D, B2);

    if(mpz_odd_p(D) != 0)
    {
//        printf("ho un D dispari\n");
        mpz_add_ui(D, D, 1);
    }

    dD = mpz_get_d(D);
    dDhalf = mpz_get_d(D) / 2;

    dMmin = floor((mpz_get_d(pd.stageOneB) + dDhalf) / dD);

    //MMIN ← (B1 +D/2)/D, MMAX ← (B2 −D/2)/D  il primo parte intera inferiore, il secondo parte superiore
//    mpz_div_ui(Dhalf, D, 2);
//    mpz_add(t1, pd.stageOneB, Dhalf);
//    mpz_div(Mmin, t1, D);

//    printf("Mmin = %lf\n", dMmin);

    dMmax = ceil((mpz_get_d(pd.stageTwoB) - dDhalf) / dD);
//
//
//    mpz_sub(t2, B2, Dhalf);
//    mpz_div(Mmax, t2, D);

//    printf("Mmax = %lf\n", dMmax);


    mpz_t gcdVal;
    mpz_init(gcdVal);

    mpz_clear(t1);
    mpz_clear(t2);


//    unsigned long arraylen = mpz_get_ui(Dhalf);
    unsigned long arraylen = (unsigned long) dDhalf;
//    mpz_t GCDtable[arraylen + 1];
    unsigned int GCDtable[arraylen + 1];

    for(unsigned long i = 1; i <= arraylen; i++) //chissà se irene è inclusa
    {
//        mpz_init(GCDtable[i]);
        mpz_gcd_ui(gcdVal, D, i);

        if(mpz_cmp_ui(gcdVal, 1) == 0)
        {
//            mpz_set_ui(GCDtable[i], 1);
            GCDtable[i] = 1;
            //add i to the set Is
            struct JsElem *cell  = malloc(sizeof(struct JsElem));
            cell->index = i;
            cell->next = head.next;
            //el = &cell;
            head.next = cell;
            //printf("aggiungo elemento %lu a lista\n", cell->index);
        }
        else
            GCDtable[i] = 0;
//        printf("| %d ", GCDtable[i]);
    }
//    printf("\n");

//    struct JsElem * in = head.next;
//    while(in != NULL)
//    {
//        printf("found el %ld in list\n", in->index);
//        in = in->next;
//        sleep(1);
//
//    }


//    mpz_t primeTableLen, m, mD, primeCandidate;
    mpz_t m, mD, primeCandidate;
    unsigned int operations;
//    mpz_init(primeTableLen);
    mpz_init(m);
    mpz_init(mD);
    mpz_init(primeCandidate);

//    mpz_sub(primeTableLen, Mmax, Mmin);
//    unsigned long primeLen = mpz_get_ui(primeTableLen);
    unsigned long primeLen = (unsigned long) (dMmax - dMmin);
//    int primaTable[primeLen + 1][arraylen + 1];
    int primaTable[primeLen][arraylen + 1];

    //printf("primelen1 = %ld\nprimelen2 = %ld\n",  primeLen + 1, arraylen);

    int check;
//    printf("%lf\n", dD);
    for(unsigned long k = 0; k < primeLen; k++)     //wrong index
    {
        for(unsigned long j = 1; j <= arraylen; j++)
        {

            operations = (unsigned int) ((dMmin + k) * dD + j);
            mpz_set_ui(primeCandidate, operations);
            // m = Mmin+k
//            mpz_add_ui(m, Mmin, k);
//            mpz_mul(mD, m, D);
//
//            mpz_add_ui(primeCandidate, mD, j);
            primaTable[k][j] = 0;
            if(mpz_probab_prime_p(primeCandidate, 20) == 2)     //the number is prime
            {
                //printf("casella con j %lu a 1 \t", j);
                primaTable[k][j] = 1;
            }
            else
            {
                check = (int) ((dMmin + k) * dD - j);
                if(check > 0)
                {
                    operations = (unsigned int) check;
                    mpz_set_ui(primeCandidate, operations);
//                mpz_sub_ui(primeCandidate, mD, j);
                    if (mpz_probab_prime_p(primeCandidate, 20) == 2) {
                        //printf("casella con j %lu a 1 \t", j);
                        primaTable[k][j] = 1;
                    }
                } else
                    operations=0;
            }
//            printf("|%d %u", primaTable[k][j], operations);
        }
//        printf("\n");
    }

    //Q = Q0
    struct ECpoint P, res;
    mpz_init(P.X);
    mpz_init(P.Z);
    mpz_init(P.Y);

    mpz_init(res.X);
    mpz_init(res.Y);
    mpz_init(res.Z);

    mpz_set(P.X, Q.X);
    mpz_set(P.Y, Q.Y);
    mpz_set(P.Z, Q.Z);

    mpz_t two;
    mpz_init(two);
    mpz_set_ui(two, 2);

    struct ECpoint points[arraylen + 1];    //DHalf
    struct nonInvertibleD d;
    mpz_init(d.d);
    d.flag = 0;

    res = ECmultiplyTraditional(&P, two, EC, pd, &d, &res);

    checkIfCurve(res, EC, pd);

//    gmp_printf("resx = %Zd\n", res.X);      //non trovo più a zero
//    gmp_printf("resy = %Zd\n", res.Y);
//    gmp_printf("resz = %Zd\n", res.Z);


    for(unsigned long k = 1; k <= arraylen; k++)  //arraylen = DHalf
    {
        mpz_init(points[k].X);
        mpz_init(points[k].Y);
        mpz_init(points[k].Z);
    }

    for(unsigned long k = 1; k <= arraylen; k = k+2)  //arraylen = DHalf
    {

//        if(mpz_cmp_ui(GCDtable[k], 1) == 0)
        if(GCDtable[k] == 1)
        {
//            gmp_printf("Px = %Zd\n", P.X);
//            gmp_printf("Py = %Zd\n", P.Y);
//            gmp_printf("Pz = %Zd\n", P.Z);          //no more: always storing the same!!!! BUG

            mpz_set(points[k].X, P.X);
            mpz_set(points[k].Y, P.Y);
            mpz_set(points[k].Z, P.Z);
/*
            gmp_printf("Pointsx = %Zd\n", points[k].X);
            gmp_printf("Pointsy = %Zd\n", points[k].Y);
            gmp_printf("Pointsz = %Zd\n", points[k].Z);*/
        }
        //Q = 2*Q0 + Q


        add2(&res, &P, EC, pd, &d, &P);

        checkIfCurve(P, EC, pd);
//        sleep(1);
    }

    //--------------------------MAIN COMPUTATION-----------------------------------
    printf("main computation for stage two\n");

    //d ← 1, P ← DQ0, R ← MminQ
    mpz_t den, q, mIndex, partial1, partial2, partial3;
    mpz_init(den);
    mpz_init(q);
    mpz_init(mIndex);
    mpz_init(partial1);
    mpz_init(partial2);
    mpz_init(partial3);

    mpz_set_ui(den, 1);

    struct ECpoint S, R, jQ0, result;
    mpz_init(S.X);
    mpz_init(S.Y);
    mpz_init(S.Z);

    mpz_init(R.X);
    mpz_init(R.Y);
    mpz_init(R.Z);

    mpz_init(jQ0.X);
    mpz_init(jQ0.Y);
    mpz_init(jQ0.Z);

    mpz_init(result.X);
    mpz_init(result.Y);
    mpz_init(result.Z);


    //printf("qui\n\n");

    //S = DQ0
    S = ECmultiplyTraditional(&Q, D, EC, pd, &d, &S);
    //R = Mmin*S
    R = ECmultiplyTraditional(&S, Mmin, EC, pd, &d, &R);

//    checkIfCurve(S, EC, pd);
//    checkIfCurve(R, EC, pd);
    //struct JsElem * el;

    volatile struct JsElem * listEl;

    //printf("primelen1 = %ld\nprimelen2 = %ld\n",  primeLen + 1, arraylen);

    mpz_set_ui(Mmin, (unsigned long) dMmin);
    mpz_set_ui(Mmax, (unsigned long) dMmax);

    for(unsigned long k = 0; k <= primeLen; k++)     //m = k+Mmin
    {
        listEl = head.next;
        mpz_add_ui(mIndex, Mmin, k);

        while(listEl != NULL)
        {
            //printf("index k = %ld\nindex i = %ld\nlen1 = %ld\nlen2 = %ld\n", k, listEl->index, primeLen + 1, arraylen);
            if(primaTable[k][listEl->index] == 1)
            {
                //printf("quiccc\n\n");
                mpz_set(jQ0.X, (points[(listEl->index)]).X);
                mpz_set(jQ0.Y, (points[(listEl->index)]).Y);
                mpz_set(jQ0.Z, (points[(listEl->index)]).Z);

                /*printf("retrieving %ld\n", listEl->index);
                gmp_printf("Pointsx = %Zd\n", (points[(listEl->index)]).X);
                gmp_printf("Pointsy = %Zd\n", (points[(listEl->index)]).Y);
                gmp_printf("Pointsz = %Zd\n", (points[(listEl->index)]).Z);

                sleep(2);*/

                //printf("naaaazapegna\n\n");

                //den ← den · (x.R*z.jQ0 − x.jQ0*z.R)
                mpz_mul(partial1, R.X, jQ0.Z);
                mpz_mul(partial2, jQ0.X, R.Z);
                mpz_sub(partial3, partial1, partial2);

                mpz_mul(den, den, partial3);

                //gmp_printf("ecco un den = %Zd\n", den);
                //sleep(1);
            }
            //sleep(1);
            //printf("time to cambiare il listEl!\n\n");
            listEl = listEl->next;

        }
//        printf("sasfasd\n");
        //sleep(1);
        add2(&R, &P, EC, pd, &d, &result);

        mpz_set(R.X, result.X);
        mpz_set(R.Y, result.Y);
        mpz_set(R.Z, result.Z);
//        printf("naaaazapegna\n\n");
        //while(scorro lista collegata di Js)
            //if(primetable[m][js] == 1)
                //faccio cose
        //R = R + Q
    }


    //mpz_mod(den, den, pd.n);


    mpz_gcd(q, den, pd.n);

    //gmp_printf("ecco a voi d = %Zd\n", den);
    gmp_printf("ecco a voi q = %Zd\n", q);
    //sleep(1);

    mpz_clear(R.X);
    mpz_clear(R.Y);
    mpz_clear(R.Z);
    mpz_clear(P.X);
    mpz_clear(P.Y);
    mpz_clear(P.Z);
    mpz_clear(partial1);
    mpz_clear(partial2);
    mpz_clear(partial3);

    for(unsigned long k = 1; k <= arraylen; k++)  //arraylen = DHalf
    {
        mpz_clear(points[k].X);
        mpz_clear(points[k].Y);
        mpz_clear(points[k].Z);
    }

    if(mpz_cmp_ui(q, 1) > 0)
    {
        if(mpz_cmp(q, pd.n) < 0)
        {
            printf("-------------------------------------------success in stage two----------------------------------\n");
            gmp_printf("\t%Zd\n---------------------------------------------------------------------------\n", q);
            //return q
            return 1;
        }
        printf("\n\n\noibo che coincidenza q è proprio n\n\n\n");
        //fail
        return 0;

    }
    else
    {
        printf("failed\n");
        //fail
        return 0;
    }
}