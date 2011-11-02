/*
 * weiner.c
 *
 * This file finds numerical solutions of the systems listed in examples but
 * implemented in c for spedd comparison.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Number of iterations after seeding. */
#define RANDITERMIN 20

int print_help(char* argv[], char* msg);
int print_weiner();

void setup_rngs(void);
double random_mod(void);

int main(int argc, char *argv[])
{
    int systype;
    int i;

    /* Initialise random seed and prepare for rng generation */
    setup_rngs();

    for (i=0; i < 20; i++)
        printf("%lf\n", random_mod());
    return 0;

    if(argc == 1)
        return print_help(argv, NULL);
    else {
        if(!sscanf(argv[1], "%i", &systype))
            return print_help(argv, "Cannot parse SYSNUM\n");
        switch(systype){
            case 0:
                printf("System %d\n", systype);
                return print_weiner();
            default:
                return print_help(argv, "Invalid SYSNUM\n");
        }
    }

    return 0;
}

int print_help(char *argv[], char *errmsg)
{
    if(errmsg != NULL)
        fprintf(stderr, errmsg);
    printf("usage: %s SYSNUM\n", argv[0]);
    printf("\n");
    printf("Numerically solve the stochastic ordinary differential\n");
    printf("equation representing a Weiner process.\n");
    return 0;
}

int print_weiner()
{
    double x1=0;
    double alpha=0, beta=1;
    double t1=0;
    double t2=1;
    double dt=0.01;

    int nsteps=(int) (t2 - t1) / dt + 1;
    int i;
    double x=x1;
    double t=t1;
    printf("t, x\n");
    printf("%lf, %lf\n", t1, x1);
    for(i=0; i<nsteps; i++) {
        x += alpha * dt + beta * 0;
        t += dt;
        printf("%lf, %lf\n", t, x);
    }
    return 0;
}

/*
 * Random number generators. First we need to generate random numbers on the
 * interval [0, 1)
 */

void setup_rngs(void) {
    /* Initialise with system time and iterate a few times */
    int i;
    srand(time(NULL));
    for(i =0; i<RANDITERMIN; i++)
        rand();

    /* Initialise cong_seed */
    cong_seed = rand();
    rng_cong();

    /* Initialise */
}

static unsigned long rng_cong_seed=123456789;
unsigned long rng_cong(void) {
    return (cong_seed=69069*cong_seed + 362437);
}

static unsigned long xorshift_x=123456789,xorshift_y=362436069,xorshift_z=521288629,xorshift_w=88675123,xor_shift_v=886756453; 
/* replace defaults with five random seed values in calling program */ 
unsigned long xorshift(void) 
{
    unsigned long t; 
    t=(x^(x>>7)); x=y; y=z; z=w; w=v; 
    v=(v^(v<<6))^(t^(t<<13)); return (y+y+1)*v;
} 

double uniform_simple(){
    return (double)rand() / ((double)RAND_MAX + 1.0);
}

double uniform_better()
