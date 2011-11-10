/*
 * For generating standard normal random variables. Adapted from
 * From Marsgalia and Tsang
 * http://www.jstatsoft.org/v05/i08/supp/1
 *
 * randnorm exports the MACROS:
 *  (1) RANDNORM_SHR3 - to generate uniform random integers in [0, 2^32]
 *  (2) RANDNORM_UNIF - to generate uniform random doubles in [0, 1]
 *  (3) RANDNORM_NORMAL - to generate (standard) normal doubles.
 *
 * Before these MACROS can be used the function randnorm_seed() must be called to
 * generate the random seed and initialise the data structures used for
 * RANDNORM_NORMAL. Otherwise RANDNORM_NORMAL will produce all zeros.
 */
#ifndef _RANDNORM_H
#define _RANDNORM_H

#include <math.h>
#include <limits.h>

/* Global variables needed by the RANDNORM_ macros */
extern unsigned long randnorm_jz, randnorm_jsr;
extern long randnorm_hz;
extern unsigned long randnorm_iz, randnorm_kn[128];
extern double randnorm_wn[128];

/* Exported functions */
void randnorm_seed_ziggurat(unsigned long jsrseed);
double randnorm_nfix(void); /* Needed by RANDNORM_NORMAL */
double randnorm_boxmuller(void);

/* Use this so that randnorm_seed chosses its own seed */
#define RANDNORM_SEED_PID_TIME 0
void randnorm_seed(unsigned int seed);

/* Macros to generate random numbers.*/
#define RANDNORM_SHR3 (\
            randnorm_jz=randnorm_jsr,\
            randnorm_jsr^=(randnorm_jsr<<13),\
            randnorm_jsr^=(randnorm_jsr>>17),\
            randnorm_jsr^=(randnorm_jsr<<5),\
          randnorm_jz+randnorm_jsr\
        )
#define RANDNORM_SHR3_MAX ULONG_MAX

#define RANDNORM_UNIF (.5 + (signed) RANDNORM_SHR3*.2328306e-9)
#define RANDNORM_NORMAL() (\
            randnorm_hz=RANDNORM_SHR3,\
            randnorm_iz=randnorm_hz&127,\
          (fabs(randnorm_hz)<randnorm_kn[randnorm_iz])\
                ? randnorm_hz*randnorm_wn[randnorm_iz] : randnorm_nfix()\
        )

#ifdef RANDNORM_DEBUG
void randnorm_debug_test();
#endif /* RANDNORM_DEBUG */

#endif /* _RANDNORM_H */
