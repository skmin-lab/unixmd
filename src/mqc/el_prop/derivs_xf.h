#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

// Routine to calculate cdot contribution originated from XF term
// CAUTION!!! different from the k_lk in CT!!!
static void xf_cdot(int nst, double **k_lk, double complex *c, double complex *xfcdot){

    double *rho = malloc(nst * sizeof(double));

    int lst, mst;

    // Calculate densities from current coefficients
    for(lst = 0; lst < nst; lst++){
        rho[lst] = creal(conj(c[lst]) * c[lst]);
    }

    // Get cdot contribution from decoherence term, state-pair expression
    for(lst = 0; lst < nst; lst++){
        xfcdot[lst] = 0.0 + 0.0 * I;
        for(mst = 0; mst < nst; mst++){
            xfcdot[lst] -= rho[mst] * k_lk[mst][lst] * c[lst];
        }
    }

}

// Routine to print xf debug info 
static void xf_print_coef(int nst, double complex *coef, double complex *xfcdot, double *dotpopdec){
    int lst;
    
    for(lst = 0; lst < nst; lst++){
        dotpopdec[lst] = 2.0 * creal(xfcdot[lst] * conj(coef[lst]));
    }
}

// CAUTION!!! different from the k_lk in CT!!!
// Routine to calculate rhodot contribution originated from XF term
static void xf_rhodot(int nst, double **k_lk, double complex **rho, double complex **xfrhodot){

    int lst, kst, mst;

    // Get rhodot contribution from decoherence term
    for(lst = 0; lst < nst; lst++){
        // Diagonal components
        xfrhodot[lst][lst] = 0.0 + 0.0 * I;
        for(mst = 0; mst < nst; mst++){
            xfrhodot[lst][lst] -= 2.0 * k_lk[mst][lst] * rho[lst][mst] * rho[mst][lst];
        }
        // Off-diagonal components
        for(kst = lst + 1; kst < nst; kst++){
            xfrhodot[lst][kst] = 0.0 + 0.0 * I;
            for(mst = 0; mst < nst; mst++){
                xfrhodot[lst][kst] -= (k_lk[mst][lst] + k_lk[mst][kst]) * rho[lst][mst] * rho[mst][kst];
            }
            xfrhodot[kst][lst] = conj(xfrhodot[lst][kst]);
        }
    }

}

// Routine to print xf debug info 
static void xf_print_rho(int nst, double complex **xfrhodot, double *dotpopdec){
    int lst;
    
    for(lst = 0; lst < nst; lst++){
        dotpopdec[lst] = creal(xfrhodot[lst][lst]);
    }
}

