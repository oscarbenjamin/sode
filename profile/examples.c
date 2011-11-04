/*
 * examples.c
 *
 * This file gives example stochastic ordinary differential equations to be
 * tested by weiner.c
 */

#include "examples.h"

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

