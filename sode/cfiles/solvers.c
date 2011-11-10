/*
 * solvers.c
 *
 * This file gives example stochastic ordinary differential equations to be
 * tested by weiner.c
 */

#include "solvers.h"
#include "randnorm.h"
#include <stdio.h>

/*
 * Euler-Mariyama integration step.
 */
void solve_EM(vecfield drift, vecfield diffusion, size_t nvars, /* SODEs */
              const double *x1, const double t1, const double dt, /* Input */
              double *x2) /* Output */
{
    size_t i;
    double temp1[nvars];
    double temp2[nvars];
    double sqrtdt = sqrt(dt);
    (*drift)(x1, t1, temp1);
    (*diffusion)(x1, t1, temp2);
    for(i=0; i<nvars; i++)
        x2[i] = x1[i] + dt * temp1[i] + sqrtdt * temp2[i] * RANDNORM_NORMAL();
}

/* Solve the equations specified by the user */
void solve_sode(vecfield drift, vecfield diffusion, const size_t nvars,
                const double *x0, const double t1, const double t2,
                const double dtmax, const double dtout,
                outputvec outfunc) {

    /* Initial conditions */
    double t = t1, tnext = t1, tout;
    double x[nvars];
    size_t v;
    for(v=0; v<nvars; v++)
        x[v] = x0[v];

    while(t < t2) {
        tout = t + dtout;
        tout = (tout >= t2) ? t2 : tout;
        while(t < tout) {
            tnext = t + dtmax;
            tnext = (tnext <= tout) ? tnext : tout;
            solve_EM(drift, diffusion, nvars, x, t, tnext-t, x);
            if(t == tnext) {
                fprintf(stderr, "Stepsize underflow\n");
                return;
            }
            t = tnext;
        }
        (*outfunc)(x, t, nvars);
    }
}
