//
// Created by root on 05/03/17.
//


#include <stdlib.h>
#include <stdio.h>
#include "dataStructures.h"


void randomEC(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd);
void stageOne(struct ellipticCurve EC, struct ECpoint Q, struct problemData pd);


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


    //loop through the next steps when there is a significant chance that there are no factors with logB digits


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

    return EXIT_SUCCESS;
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

            }
        }

    }
    //gcd
}
