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
#include "examples.h"
#include "solvers.h"

/* Struct to hold the data returned by parse_args */
typedef struct InputOptionsStruct {
    SysInfo *sysinfo;
    double t1;
    double t2;
    double dtmax;
    double dtout;
} InputOptions;


/* Function used to process the command line args. */
#define ARGUMENTS_INVALID 0
#define ARGUMENTS_GOOD 1
int parse_args(int argc, char *argv[], InputOptions *opts);

#define SOLVE_FAILED 0
#define SOLVE_SUCCEEDED 1
int print_sode(InputOptions *opts);

int main(int argc, char *argv[])
{
    InputOptions opts;

    /* Initialise random seed and prepare for rng generation */
    randnorm_seed(RANDNORM_SEED_PID_TIME);

    /* Parse the input arguments */
    if(parse_args(argc, argv, &opts) == ARGUMENTS_INVALID)
        return 1;

    /* Actually solve and write the result to stdout */
    print_sode(&opts);

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
    return ARGUMENTS_INVALID;
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
    return ARGUMENTS_GOOD;
}

void outfunc(const double *x, const double t, const size_t nvars)
{
    size_t v;
    printf("%f", t);
    for(v=0; v<nvars; v++)
        printf(", %f", x[v]);
    printf("\n");
}

/* Solve the equations specified by the user */
int print_sode(InputOptions *opts) {

   /*         void solve_sode(vecfield drift, vecfield diffusion, size_t nvars,
                            const double *x0, const double t0,
                            const double dtmax, const double dtout,
                            outputvec outfunc) {*/

    SysInfo *sys = opts->sysinfo;
    solve_sode(sys->drift, sys->diffusion, sys->nvars, sys->x0,
               opts->t1, opts->t2, opts->dtmax, opts->dtout, outfunc);
    return SOLVE_SUCCEEDED;
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
    double sqrtdt = sqrt(dt);

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
        x += alpha * dt + beta * sqrtdt * RANDNORM_NORMAL();
        t += dt;
        printf("%f, %f\n", t, x);
    }
    return 0;
}
