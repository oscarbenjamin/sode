/*
 * weiner.c
 *
 * This file finds numerical solutions of the systems listed in examples but
 * implemented in c for spedd comparison.
 */

#include <stdio.h>


char *progname;
int print_help(char* msg);


int main(int argc, char *argv[])
{
    int systype;

    progname = argv[0];

    if(argc == 1)
        return print_help(NULL);
    else {
        if(!sscanf(argv[1], "%i", &systype))
            return print_help("Cannot parse SYSNUM\n");
        switch(systype){
            case 0:
                printf("You have selected 0\n");
                return 0;
            default:
                return print_help("Invalid SYSNUM\n");
        }
    }

    return 0;
}

int print_help(char *errmsg)
{
    if(errmsg != NULL)
        fprintf(stderr, errmsg);
    printf("usage: %s SYSNUM\n", progname);
    printf("\n");
    printf("Numerically solve the stochastic ordinary differential\n");
    printf("equation representing a Weiner process.\n");
    return 0;
}
