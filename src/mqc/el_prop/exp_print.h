#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

// Routine to print xf debug info for exponential propagator
static void exp_xf_print_coef(int nst, double complex *coef, double **dec, double *dotpopdec){
    int ist, jst;
 
    for(ist = 0; ist < nst; ist++){
        dotpopdec[ist] = 0.0;
        for(jst = 0; jst < nst; jst++){
            dotpopdec[ist] -= 2 * creal(dec[jst][ist] * conj(coef[ist]) * coef[ist] * conj(coef[jst]) * coef[jst]);
        }
    }
}


