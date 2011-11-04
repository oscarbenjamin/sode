/*
 * weiner.c
 *
 * This file finds numerical solutions of the systems listed in examples but
 * implemented in c for speed comparison.
 */

#include <stdio.h>

/* For generating standard normals. */
#include "randnorm.h"

/* Main program routines */
int print_help(char* argv[], char* msg);
int print_weiner();

int main(int argc, char *argv[])
{
    int systype;
    int i;

    /* Initialise random seed and prepare for rng generation */
    randnorm_seed(RANDNORM_SEED_PID_TIME);

    for (i=0; i < 1000000; i++)
        printf("%.15f\n", RAND_NORMAL);
    return 0;

    if(argc == 1)
        return print_help(argv, NULL);
    else {
        if(!sscanf(argv[1], "%i", &systype))
            return print_help(argv, "Cannot parse SYSNUM\n");
        switch(systype){
            case 0:
                printf("System %d\n", systype);
                return print_weiner();
            default:
                return print_help(argv, "Invalid SYSNUM\n");
        }
    }

    return 0;
}

int print_help(char *argv[], char *errmsg)
{
    if(errmsg != NULL)
        fprintf(stderr, errmsg);
    printf("usage: %s SYSNUM\n", argv[0]);
    printf("\n");
    printf("Numerically solve the stochastic ordinary differential\n");
    printf("equation representing a Weiner process.\n");
    return 0;
}

int print_weiner()
{
    double x1=0;
    double alpha=0, beta=1;
    double t1=0;
    double t2=1;
    double dt=0.01;

    int nsteps=(int) (t2 - t1) / dt + 1;
    int i;
    double x=x1;
    double t=t1;
    printf("t, x\n");
    printf("%f, %f\n", t1, x1);
    for(i=0; i<nsteps; i++) {
        x += alpha * dt + beta * 0;
        t += dt;
        printf("%f, %f\n", t, x);
    }
    return 0;
}
