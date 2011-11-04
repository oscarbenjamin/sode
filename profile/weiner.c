/*
 * weiner.c
 *
 * This file finds numerical solutions of the systems listed in examples but
 * implemented in c for speed comparison.
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

/* For generating standard normals. */
#include "randnorm.h"

/* Enumerate the possible systems we can define */
typedef enum {LINEAR2D, LINEARMULT, SINMULT, TANMULT, WEINER, WEINER2} SysType;

/* Struct to hold the data returned by parse_args */
#define INVALID_ARGUMENTS 0
#define GOOD_ARGUMENTS 1
typedef struct InputOptionsStruct {
    SysType systype;
    double t1;
    double t2;
    double dtmax;
    double dtout;
} InputOptions;

/* Function used to process the command line args. */
int parse_args(int argc, char *argv[], InputOptions *opts);

int main(int argc, char *argv[])
{
    InputOptions opts;

    /* Initialise random seed and prepare for rng generation */
    randnorm_seed(RANDNORM_SEED_PID_TIME);

    if(parse_args(argc, argv, &opts) == INVALID_ARGUMENTS)
        return 1;
    else
        return 0;
}

int print_help(int argc, char *argv[], char *errmsg)
{
    if(errmsg != NULL)
        fprintf(stderr, "%s\n", errmsg);
    printf("usage: %s SYSNUM\n", argv[0]);
    printf("\n");
    printf("Numerically solve the stochastic ordinary differential\n");
    printf("equation representing a Weiner process.\n");
    return 0;
}

int invalid_arguments(int argc, char *argv[], char *errmsg) {
    print_help(argc, argv, errmsg);
    return INVALID_ARGUMENTS;
}

/* Used by main to process the input arguments */
int parse_args(int argc, char *argv[], InputOptions *opts) {
    /* Check the number of arguments first */
    if(argc != 1 + 1)
        return print_help(argc, argv, "Wrong number of arguments");

    /* Check the arguments in order */
    argv++;
    if(!strcmp(*argv, "weiner"))
        opts->systype = WEINER;
    else
        return invalid_arguments(argc, argv, "Unrecognised sysname");
    return GOOD_ARGUMENTS;
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
    double t2=100;
    double dt=0.0001;
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
