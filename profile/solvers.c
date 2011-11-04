/*
 * solvers.c
 *
 * This file gives example stochastic ordinary differential equations to be
 * tested by weiner.c
 */

#include "solvers.h"
#include "randnorm.h"


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
        x2[i] = x1[i] + dt * temp1[i] + sqrtdt * temp2[i] * RANDNORM_NORMAL;
}
