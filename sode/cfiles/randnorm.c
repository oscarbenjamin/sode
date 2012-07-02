/*
 * For generating random normals. The Box Muller and Ziggurat algorithms
 * adapted from Jeremy Lea and Marsaglia and Tsang.
 */

/* This is to use OS-specific stuff when initialising the seed */
#if WIN32 | _WIN32 | __WIN32__
    #define WINDOWS 1
#elif __APPLE__ & __MACH__
    #define OSX 1
#elif __unix__ | __posix__ | __linux__
    #define NIX 1
#else
    /* Unrecognised system so just use the standard c time() function */
#endif

#include <math.h>
#include <time.h>

/* Debug function needs printf */
#if RANDNORM_DEBUG
#include <stdio.h>
#endif

/* Random seed uses pid and time in us or ns */
#if WINDOWS
    #include <process.h>
    #include <windows.h>
#elif NIX
    #include <unistd.h>
#elif OSX
    #include <unistd.h>
    #include <mach/clock.h>
    #include <mach/mach.h>
#endif


/* To get the MACROs defined there */
#include "randnorm.h"

/* Used in seed initialisation */
#define RANDNORM_CONG(n) (69069*n+1234567)

/* Global variables used by RANDNORM_ macros */
unsigned long randnorm_jz, randnorm_jsr;
long randnorm_hz;
unsigned long randnorm_iz, randnorm_kn[128];
double randnorm_wn[128],randnorm_fn[128];

/*
 * First, get a good initial seed. If possible randnorm_seed attempts to use
 * the time in microseconds or nanoseconds and process id in conjunction with
 * the time in seconds. These are put through a congruential generator and
 * then XOR'ed. If not possible time() is used.
 */

/* Used by randnorm_seed, but not exported by randnorm.h */
unsigned int initial_seed(void);
void randnorm_seed_ziggurat(unsigned long jsrseed);

/* Must be called by programs using randnorm.h */
void randnorm_seed(unsigned int seed) {
    seed = (seed != RANDNORM_SEED_PID_TIME) ? seed : initial_seed();
    randnorm_seed_ziggurat(seed);
}

/* Use OS-specific functions to initialise seed if possible */
unsigned int initial_seed(void) {
    unsigned int seed, pid, tprec;
#if NIX
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    pid = getpid();
    tprec = ts.tv_nsec;
#elif OSX
    /* OSX doesn't have clock_gettime(): use clock_get_time() */
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    pid = getpid();
    tprec = mts.tv_nsec;
#elif WINDOWS
    double microseconds;
    LARGE_INTEGER clockfreq, t;
    if(!QueryPerformanceFrequency(&clockfreq))
        tprec = 0;
    else {
        QueryPerformanceCounter(&t);
        microseconds = (double)t.LowPart / ((double)clockfreq.QuadPart / 1e6);
        tprec = ((unsigned int) microseconds) % 1000000;
    }
    pid = _getpid();
#endif
    seed = RANDNORM_CONG(time(NULL));
#if OSX | NIX | WINDOWS
    seed ^= RANDNORM_CONG(pid);
    seed ^= RANDNORM_CONG(tprec);
#endif
    return seed;
}

/*
 * Box Muller method
 *
 * From Jeremy Lea:
 * http://home.online.no/~pjacklam/notes/invnorm/impl/lea/lea.c
 *
 * A normally distributed random number generator. We avoid
 * the uniform rv's being 0.0 since this will result in infinte
 * values, and double count the 0 == 2pi.
 */

double randnorm_boxmuller() {
    static int i = 1;
    static double u[2] = {0.0, 0.0};
    register double r1, r2, t;

    if (i == 1) {
        t = RANDNORM_UNIF;
        t = t ? t : 6.2831853071796;
        r1 = sqrt(-2 * log((double)(t)));
        r2 = 6.2831853071796 * (double)RANDNORM_UNIF;
        u[0] = r1 * sin(r2);
        u[1] = r1 * cos(r2);
        i = 0;
    } else {
        i = 1;
    }

    return u[i];
}

/*
 * Ziggurat algorithm
 *
 * From Marsgalia and Tsang:
 * http://www.jstatsoft.org/v05/i08/supp/1
 *
 * Combine the code below with the main program in which you want normal or
 * exponential variates. Then use of RANDNORM_NORMAL in any expression will
 * provide a standard normal variate with mean zero, variance 1. Before using
 * RANDNORM_NORMAL in your main, insert a command such as
 * randnorm_seed_ziggurat(86947731); with your own choice of seed value>0,
 * rather than 86947731. (If you do not invoke randnorm_seed_ziggurat(...) you
 * will get all zeros for RANDNORM_NORMAL.) For details of the method, see
 * Marsaglia and Tsang, "The ziggurat method for generating random variables",
 * Journ.  Statistical Software.
 */

/*
 * randnorm_nfix() generates variates from the residue when rejection in
 * RANDNORM_NORMAL occurs.
 */

double randnorm_nfix(void)
{
    const double r = 3.442620;     /* The start of the right tail */
    static double x, y;
    for(;;)
    {
        x=randnorm_hz*randnorm_wn[randnorm_iz];
        /* randnorm_iz==0, handles the base strip */
        if(randnorm_iz==0)
        {
            do{
                x=-log(RANDNORM_UNIF)*0.2904764; y=-log(RANDNORM_UNIF);
            } while(y+y<x*x);  /* .2904764 is 1/r */
            return (randnorm_hz>0)? r+x : -r-x;
        }
        /* randnorm_iz>0, handle the wedges of other strips */
        if( randnorm_fn[randnorm_iz]
            + RANDNORM_UNIF*(randnorm_fn[randnorm_iz-1]-randnorm_fn[randnorm_iz])
                < exp(-.5*x*x) )
            return x;

        /* initiate, try to exit for(;;) for loop*/
        randnorm_hz=RANDNORM_SHR3;
        randnorm_iz=randnorm_hz&127;
        if(fabs(randnorm_hz)<randnorm_kn[randnorm_iz])
            return (randnorm_hz*randnorm_wn[randnorm_iz]);
    }
}

/*--------This procedure sets the seed and creates the tables------*/

void randnorm_seed_ziggurat(unsigned long jsrseed)
{
    const double m1 = RANDNORM_SHR3_MAX / 2;
    double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
    int i;
    randnorm_jsr = 123456789 ^ jsrseed;

    /* Set up tables for RANDNORM_NORMAL */
    q=vn/exp(-.5*dn*dn);
    randnorm_kn[0]=(dn/q)*m1;
    randnorm_kn[1]=0;

    randnorm_wn[0]=q/m1;
    randnorm_wn[127]=dn/m1;

    randnorm_fn[0]=1.;
    randnorm_fn[127]=exp(-.5*dn*dn);

    for(i=126;i>=1;i--) {
        dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
        randnorm_kn[i+1]=(dn/tn)*m1;
        tn=dn;
        randnorm_fn[i]=exp(-.5*dn*dn);
        randnorm_wn[i]=dn/m1;
    }
#ifdef RANDNORM_DEBUG
    randnorm_debug_test();
#endif /* RANDNORM_DEBUG */
}

#ifdef RANDNORM_DEBUG

#define PRINT_SIZEOF(x) printf("sizeof(" #x "): %d bytes\n", (int)sizeof(x))
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define MIN(a, b) ((a) < (b)) ? (a) : (b)

#ifndef RANDNORM_NTIMES
#define RANDNORM_NTIMES 1000000
#endif /* RANDNORM_NTIMES */

/* MACRO to test the statistics of a larg number of random numbers */
#define MINMAXMEANVAR(type, fmt, name, expr, description, details)      \
    type name##_min, name##_max, name##_current;                        \
    double name##_mean, name##_var;                                     \
    name##_min = name##_max = expr;                                     \
    name##_mean = (double)name##_min;                                   \
    name##_var = 0;                                                     \
    for(i=2; i<=RANDNORM_NTIMES; i++) {                                 \
        id = (double) i;                                                \
        name##_current = expr;                                          \
        name##_min = MIN(name##_min, name##_current);                   \
        name##_max = MAX(name##_max, name##_current);                   \
        name##_mean = ((id-1.)/id)*name##_mean + (1./id)*(double)name##_current;\
        dev = (double)name##_current - name##_mean;                     \
        name##_var = ((id-1.)/id)*name##_var + (1./(id - 1.))*dev*dev;  \
    }                                                                   \
    printf("\tMin: " fmt " Max: " fmt " Mean: %le Var: %le\n",          \
            name##_min, name##_max, name##_mean, name##_var);           \
    printf("True values for %s:\n\t%s\n", description, details)

/* MACRO to time the generation of random numbers using RANDNORM macros */
#define TIMEREPEAT(type, name, expr)                                    \
    clock_t name##_t1, name##_t2;                                       \
    name##_t1 = clock();                                                \
    type name##_total = 0;                                              \
    for(i=1; i<=RANDNORM_NTIMES; ++i)                                   \
        name##_total += expr;                                           \
    name##_t2 = clock();                                                \
    double name##_tsecs = ((double)(clock() - name##_t1))               \
                                        / ((double)CLOCKS_PER_SEC);     \
    printf("%.4lf seconds to make and add %d numbers\n",                \
            name##_tsecs, RANDNORM_NTIMES);                             \
    printf("%.1lf nanoseconds per eval\n",                              \
            name##_tsecs * 1.0e9 / ((double)RANDNORM_NTIMES) )          \

/* MACRO to combine the above and print a header */
#define TEST_RANDNORM_MACRO(type, fmt, name, expr, description, details)\
    printf("\n\n-------%s   %s\n\n", #name, description);               \
    printf("3 examples:");                                              \
    for(i=1; i<=3; ++i)                                                 \
        printf(" " fmt, expr);                                          \
    printf("\nStats from %u examples of %s:\n", RANDNORM_NTIMES, #name);\
    MINMAXMEANVAR(type, fmt, name, expr, description, details);         \
    TIMEREPEAT(type, name, expr);


void randnorm_debug_test()
{
    /* Separate from earlier output */
    printf("\n\n\nRANDNORM_DEBUG_TEST: ------------BEGIN---------\n\n");

    /* print the size of all standard types */
    PRINT_SIZEOF(char);
    PRINT_SIZEOF(unsigned char);
    PRINT_SIZEOF(short);
    PRINT_SIZEOF(unsigned short);
    PRINT_SIZEOF(int);
    PRINT_SIZEOF(unsigned int);
    PRINT_SIZEOF(long);
    PRINT_SIZEOF(unsigned long);
    PRINT_SIZEOF(long long);
    PRINT_SIZEOF(unsigned long long);
    PRINT_SIZEOF(float);
    PRINT_SIZEOF(double);

    printf("RANDNORM_SHR3_MAX: %lu", RANDNORM_SHR3_MAX);

    /* collect stats on a large number of random variables */
    int i; double id, dev;/* Used by MINMAXMEANVAR */

    char shr3_details[100];
    sprintf(shr3_details, "Min: 0 Max: %lu Mean: %le Var: %le",
            RANDNORM_SHR3_MAX, (double)RANDNORM_SHR3_MAX / 2.,
            (((double)RANDNORM_SHR3_MAX * (double)RANDNORM_SHR3_MAX) - 1.) / 12.);

    TEST_RANDNORM_MACRO(unsigned long, "%lu", shr3, RANDNORM_SHR3,
                        "uniform integers",
                        shr3_details);

    TEST_RANDNORM_MACRO(double, "%lf", unif, RANDNORM_UNIF,
                        "uniform reals",
                        "Min: 0 Max: 1 Mean: 0.5 Var: 0.83333333");

    TEST_RANDNORM_MACRO(double, "%lf", normal, RANDNORM_NORMAL,
                        "standard normals",
                        "Min: -Inf Max: +Inf Mean: 0 Var: 1");

    TEST_RANDNORM_MACRO(double, "%lf", boxm, randnorm_boxmuller(),
                        "standard normals",
                        "Min: -Inf Max: +Inf Mean: 0 Var: 1");

    /* Separate from later output */
    printf("\n\n\nRANDNORM_DEBUG_TEST: -------------END----------\n\n");
}

#endif /* RANDNORM_DEBUG */
