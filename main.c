//
// Created by root on 05/03/17.
//


#include <stdlib.h>
#include <stdio.h>
#include "dataStructures.h"
#include <unistd.h>

void classicalECM(struct problemData pd, mpz_t factor);
void traditionalStageOne(struct weirstrassEC EC, struct ECpoint Q, struct problemData pd, mpz_t factor);


        void randomEC(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd);
void stageOne(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd);
void randomECtraditional(struct weirstrassEC * EC, struct ECpoint  *Q, struct problemData pd);


int main(int argc, char ** argv)
{
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
    mpz_set_ui(pd.stageOneB, 1000);


    mpz_sqrt(pd.stageTwoB, pd.n);



    //sleep(1);

    mpz_t factor;
    mpz_init(factor);

    while(mpz_cmp(pd.stageOneB, pd.stageTwoB) < 0)
    {
        classicalECM(pd, factor);
        printf("failure, possibly increment B\n\n");
        mpz_mul_ui(pd.stageOneB, pd.stageOneB, 10);
        sleep(3);
    }
    printf("the size of B1 exceeded max B");

    //loop through the next steps when there is a significant chance that there are no factors with logB digits





    return EXIT_SUCCESS;
}

void classicalECM(struct problemData pd, mpz_t factor) {
    struct ECpoint Q;
    struct weirstrassEC EC;

    mpz_init(Q.X);
    mpz_init(Q.Y);
    mpz_init(Q.Z);

    mpz_set_ui(Q.Z, 1);

    mpz_init(EC.a);
    mpz_init(EC.b);

    randomECtraditional(&EC, &Q, pd);
    printf("random init terminated\n");

    gmp_printf("%Zd\n", Q.X);

    gmp_printf("%Zd\n", Q.Y);

    gmp_printf("%Zd\n", Q.Z);

    gmp_printf("%Zd\n", EC.a);

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
        classicalECM(pd, factor);
    }
    else if (mpz_cmp_ui(g, 1) > 0)
    {
        //return gcdterm as factor
        printf("found factor with random choice\n");
        gmp_printf("%Zd\n\n", g);
    }
    else
    {
        //prime power multipliers
        printf("starting step 1 in seconds\n");
        //sleep(2);
        traditionalStageOne(EC, Q, pd, factor);

    }
}

void invertionlessECM(struct problemData pd)
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

    /*
     * int mpz_cmp_ui (const mpz t op1, unsigned long int op2) [Macro]
Compare op1 and op2. Return a positive value if op1 > op2, zero if op1 = op2, or a negative
value if op1 < op2.
     * */

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
}

void traditionalStageOne(struct weirstrassEC EC, struct ECpoint Q, struct problemData pd, mpz_t  factor)
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
    printf("i'm here\n\n");

    struct ECpoint P;

    mpz_init(P.X);
    mpz_init(P.Y);
    mpz_init(P.Z);

    mpz_set(P.X, Q.X);
    mpz_set(P.Y, Q.Y);
    mpz_set(P.Z, Q.Z);


    gmp_printf("x = %Zd , y= %Zd , z= %Zd \n", P.X, P.Y, P.Z);
    //sleep(3);
    while(mpz_cmp(primen, pd.stageOneB) <= 0)
    {
        //perform operations
        //find largest integer a such that pi^a <= B1
        unsigned long exp = 1;
        int flag = 1;
        mpz_t power;
        mpz_init(power);
        printf("i'm inside the loop 1\n");
        while(flag)
        {
            mpz_pow_ui(power, primen, exp);
            printf("i'm inside the loop 2\n");
            if(mpz_cmp(power, pd.stageOneB) <= 0)
            {
                exp += 1;
                //looking for a bigger exponent
            }
            else
            {
                flag = 0;
                printf("the exponent is %lu\n", exp);
            }
        }


        unsigned int i;

        struct nonInvertibleD d;
        mpz_init(d.d);
        d.flag = 0;

        printf("starting prime power multipliers\n\n");
        for(i = 1; i <= exp; i++)
        {
            //Q = [pi]Q

            gmp_printf("moltiplico per il primo %Zd\n", primen);
            P = ECmultiplyTraditional(&P, primen, EC, pd, &d);

            if(d.flag)
            {
                mpz_t newFactor;
                mpz_init(newFactor);
                mpz_gcd(newFactor, pd.n, d.d);
                gmp_printf("found factor %Zd\n", newFactor);
                factor = newFactor;
                return;
            }

            //the cycle stops if I find a non invertible denominator in the addition slope
        }
        printf("on with another prime\n");


        //find next prime
        mpz_nextprime(primen, primen);
        gmp_printf("%Zd\n\n", primen);

    }
    printf("failure with this size of B\n");
}

void randomECtraditional(struct weirstrassEC * EC, struct ECpoint * Q, struct problemData pd)
{
    //random sigma, derive u v anc C from this

    //initialize random state
    gmp_randstate_t state;
    gmp_randinit_mt(state);

    //generate random x, y, a in [0, n-1]

    mpz_urandomm(Q->X, state, pd.n);
    gmp_printf("%Zd\n", Q->X);
    mpz_urandomm(Q->Y, state, pd.n);
    gmp_printf("%Zd\n", Q->Y);
    mpz_urandomm(EC->a, state, pd.n);
    gmp_printf("%Zd\n", EC->a);
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
    gmp_printf("%Zd\n", EC->b);

}