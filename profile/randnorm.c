/*
 * For generating random normals. The Box Muller and Ziggurat algorithms
 * adapted from Jeremy Lea and Marsaglis and Tsang.
 */

/* This is to use OS-specific stuff when initialising the seed */
#if WIN32 | _WIN32 | __WIN32__
    #define WINDOWS 1
#elif __unix__ | __posix__ | __linux__ | __APPLE__
    #define NIX 1
#else
    /* Unrecognised system so just use the standard c time() function */
#endif

#include <math.h>
#include <time.h>

/* Random seed uses pid and time in us or ns */
#if WINDOWS
    #include <process.h>
    #include <windows.h>
#else
    #include <unistd.h>
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
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    pid = getpid();
    tprec = s.tv_nsec;
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
#if NIX | WINDOWS
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
    const double m1 = 2147483648.0;
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
}
