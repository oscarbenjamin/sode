/*
 * solvers.h
 *
 * Header file exporting the Stochastic ODE solvers used in weiner.c
 */
#ifndef _SOLVERS_H
#define _SOLVERS_H

/* for size_t */
#include <stddef.h>

/* signature for vector field functions used by the solvers */
typedef void (*vecfield) (const double*, const double, double*);

/* a single Euler Mariyama step */
void solve_EM(vecfield drift, vecfield diffusion, size_t nvars, /* SODEs */
              const double *x1, const double t1, const double dt, /* Input */
              double *x2); /* Output */

/* signature for output functions passed in to solver */
typedef void (*outputvec) (const double*, const double, const size_t);

/* Solve the equations specified by the user */
void solve_sode(vecfield drift, vecfield diffusion, const size_t nvars,
                const double *x0, const double t1, const double t2,
                const double dtmax, const double dtout,
                outputvec outfunc);

#endif /* _SOLVERS_H */


