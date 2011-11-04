/*
 * solvers.h
 *
 * Header file exporting the Stochastic ODE solvers used in weiner.c
 */
#ifndef _SOLVERS_H
#define _SOLVERS_H

/* for size_t */
#include <stddef.h>

typedef void (*vecfield) (const double*, const double, double*);

void solve_EM(vecfield drift, vecfield diffusion, size_t nvars, /* SODEs */
              const double *x1, const double t1, const double dt, /* Input */
              double *x2); /* Output */

#endif /* _SOLVERS_H */


