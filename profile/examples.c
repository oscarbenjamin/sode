/*
 * examples.c
 *
 * This file gives example stochastic ordinary differential equations to be
 * tested by weiner.c
 */

#include "examples.h"

/* To determine the size of an array */
#define COUNT_OF(x) ((sizeof(x)/sizeof(0[x])) / ((size_t)(!(sizeof(x) % sizeof(0[x])))))

const double weiner_x0[] = {0};
void weiner_drift(const double *x1, const double t1, double *a) {
    a[0] = 0;
}
void weiner_diffusion(const double *x1, const double t1, double *b) {
    b[0] = 1;
}


const double weiner2_x0[] = {0};

const double lin2d_x0[] = {0};

const double linmult_x0[] = {0};

const double sinmult_x0[] = {0};

const double tanmult_x0[] = {0};


/* Set of systems to choose from */
SysInfo example_systems[] = {
    {"weiner" ,  weiner_x0, COUNT_OF( weiner_x0), &weiner_drift, &weiner_diffusion},
    {"weiner2", weiner2_x0, COUNT_OF(weiner2_x0), NULL, NULL},
    {  "lin2d",   lin2d_x0, COUNT_OF(  lin2d_x0), NULL, NULL},
    {"linmult", linmult_x0, COUNT_OF(linmult_x0), NULL, NULL},
    {"sinmult", sinmult_x0, COUNT_OF(sinmult_x0), NULL, NULL},
    {"tanmult", tanmult_x0, COUNT_OF(tanmult_x0), NULL, NULL},
};
const unsigned int num_systems =
        sizeof(example_systems) / sizeof(example_systems[0]);

