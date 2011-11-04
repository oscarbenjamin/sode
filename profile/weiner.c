/*
 * weiner.c
 *
 * This file finds numerical solutions of the systems listed in examples but
 * implemented in c for speed comparison.
 */

#include <stdio.h>
#include <math.h>
#include <time.h>

/* For generating standard normals. */
#include "randnorm.h"

/* Main program routines */
int print_help(char* argv[], char* msg);
int print_weiner();

int main(int argc, char *argv[])
{
    int systype;

    /* Initialise random seed and prepare for rng generation */
    randnorm_seed(RANDNORM_SEED_PID_TIME);

    if(argc == 1)
        return print_help(argv, NULL);
    else {
        if(!sscanf(argv[1], "%i", &systype))
            return print_help(argv, "Cannot parse SYSNUM\n");
        switch(systype){
            case 0:
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

void sprint_time(char *timestr, size_t max) {
    time_t raw_time;
    struct tm *timeinfo;
    time(&raw_time);
    timeinfo = localtime(&raw_time);
    strftime(timestr, max, "%a, %d %b %Y %H:%M:%S %z", timeinfo);
}

int print_weiner()
{
    double x1=0;
    double alpha=0, beta=1;
    double t1=0;
    double t2=1;
    double dt=0.01;
    double sqrtdt=sqrt(dt);

    int nsteps=(int) ceil((t2 - t1) / dt);
    int i;
    double x=x1;
    double t=t1;
    char timestr[100];
    sprint_time(timestr, 100);
    printf("# Created by weiner.c at %s\n", timestr);
    printf("time, x\n");
    printf("%f, %f\n", t1, x1);
    for(i=0; i<nsteps; i++) {
        x += alpha * dt + beta * sqrtdt * RANDNORM_NORMAL;
        t += dt;
        printf("%f, %f\n", t, x);
    }
    return 0;
}
