/*
 * examples.h
 *
 * Header file exporting the Stochastic ODEs listed in examples.c
 */
#ifndef _EXAMPLES_H
#define _EXAMPLES_H

/* for size_t */
#include <stddef.h>

/* for vecfield */
#include "solvers.h"

/* Struct to define the features of the system */
typedef struct SysInfoStruct {
    char *sysname;
    const double *x0;
    const size_t nvars;
    vecfield drift;
    vecfield diffusion;
} SysInfo;

/* Set of systems to choose from */
extern SysInfo example_systems[];
extern const unsigned int num_systems;

#endif /* _EXAMPLES_H */

