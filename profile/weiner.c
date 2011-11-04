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

/* Struct to define the features of the system */
typedef struct SysInfoStruct {
    char *sysname;
} SysInfo;

/* Set of systems to choose from */
SysInfo example_systems[] = {
    {"weiner"     },
    {"linear2d"   },
    {"linearmult" },
    {"sinmult"    },
    {"tanmult"    },
    {"weiner"     },
    {"weiner2"    },
};
const unsigned int num_systems =
        sizeof(example_systems) / sizeof(example_systems[0]);

/* Struct to hold the data returned by parse_args */
#define INVALID_ARGUMENTS 0
#define GOOD_ARGUMENTS 1
typedef struct InputOptionsStruct {
    SysInfo *sysinfo;
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
    printf("usage: %s SYSNUM T1 T2 DTMAX DTOUT\n", argv[0]);
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
    int n;
    char *sysstr, *t1str, *t2str, *dtmaxstr, *dtoutstr;

    /* Check the number of arguments first */
    if(argc != 5 + 1)
        return print_help(argc, argv, "Wrong number of arguments");

    /* Separate into named strings */
    sysstr = argv[1];
    t1str = argv[2];
    t2str = argv[3];
    dtmaxstr = argv[4];
    dtoutstr = argv[5];

    /* Check for the system name */
    opts->sysinfo = NULL;
    for(n=0; n<num_systems; n++)
        if(!strcmp(example_systems[n].sysname, sysstr))
            opts->sysinfo = &example_systems[n];
    if(opts->sysinfo == NULL)
        return invalid_arguments(argc, argv, "Unrecognised SYSNAME");

    /* Parse numeric arguments */

    if(!sscanf(t1str, "%lf", &opts->t1))
        return invalid_arguments(argc, argv, "Unable to parse T1");

    if(!sscanf(t2str, "%lf", &opts->t2))
        return invalid_arguments(argc, argv, "Unable to parse T2");

    if(!sscanf(dtmaxstr, "%lf", &opts->dtmax))
        return invalid_arguments(argc, argv, "Unable to parse DTMAX");

    if(!sscanf(dtoutstr, "%lf", &opts->dtout))
        return invalid_arguments(argc, argv, "Unable to parse DTOUT");

    /* Indicate success */
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
