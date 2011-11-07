/*
 * examples.c
 *
 * This file gives example stochastic ordinary differential equations to be
 * tested by weiner.c
 */

#include <math.h>

#include "examples.h"

/* To determine the size of an array */
#define COUNT_OF(x) ((sizeof(x)/sizeof(0[x])) / ((size_t)(!(sizeof(x) % sizeof(0[x])))))

/* Plain weiner process */
const double weiner_x0[] = {0};
void weiner_drift(const double *x, const double t1, double *a) {
    a[0] = 0;
}
void weiner_diffusion(const double *x, const double t1, double *b) {
    b[0] = 1;
}

/* 2d linear SODE */
const double lin2d_x0[] = {0, 0};
#define LIN2D_BETA 0.01
#define LIN2D_ALPHA 1
void lin2d_drift(const double *x, const double t1, double *a) {
    a[0] = - LIN2D_BETA * x[0] + LIN2D_ALPHA * x[1];
    a[1] = - LIN2D_BETA * x[1] - LIN2D_ALPHA * x[0];
}
void lin2d_diffusion(const double *x, const double t1, double *b) {
    b[0] = b[1] = 1;
}

/* linear SODE with multiplicative noise */
const double linmult_x0[] = {1};
#define LINMULT_ALPHA -0.0001
#define LINMULT_BETA 1
void linmult_drift(const double *x, const double t1, double *a) {
    a[0] = LINMULT_ALPHA * x[0];
}
void linmult_diffusion(const double *x, const double t1, double *b) {
    b[0] = LINMULT_BETA * x[0];
}

/* nonlinear SODE with multiplicative noise */
const double sinmult_x0[] = {0};
#define SINMULT_ALPHA 0.01
void sinmult_drift(const double *x, const double t1, double *a) {
    a[0] = - SINMULT_ALPHA * SINMULT_ALPHA * x[0] / 2;
}
void sinmult_diffusion(const double *x, const double t1, double *b) {
    b[0] = SINMULT_ALPHA * sqrt(1 - x[0]*x[0]);
}

/* nonlinear SODE with multiplicative noise */
const double tanmult_x0[] = {0};
#define TANMULT_ALPHA 0.01
void tanmult_drift(const double *x, const double t1, double *a) {
    a[0] = TANMULT_ALPHA * TANMULT_ALPHA * x[0] * (1 - x[0]*x[0]);
}
void tanmult_diffusion(const double *x, const double t1, double *b) {
    b[0] = TANMULT_ALPHA * (1 + x[0]*x[0]);
}


/* Set of systems to choose from */
SysInfo example_systems[] = {
    { "weiner",  weiner_x0, COUNT_OF( weiner_x0),  &weiner_drift,  &weiner_diffusion},
    {"weiner2",  weiner_x0, COUNT_OF( weiner_x0),  &weiner_drift,  &weiner_diffusion},
    {  "lin2d",   lin2d_x0, COUNT_OF(  lin2d_x0),   &lin2d_drift,   &lin2d_diffusion},
    {"linmult", linmult_x0, COUNT_OF(linmult_x0), &linmult_drift, &linmult_diffusion},
    {"sinmult", sinmult_x0, COUNT_OF(sinmult_x0), &sinmult_drift, &sinmult_diffusion},
    {"tanmult", tanmult_x0, COUNT_OF(tanmult_x0), &tanmult_drift, &tanmult_diffusion},
};
const unsigned int num_systems =
        sizeof(example_systems) / sizeof(example_systems[0]);

