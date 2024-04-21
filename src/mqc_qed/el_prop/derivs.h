#ifndef DERIVS_H
#define DERIVS_H

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

//// Routine to calculate dot product from two temporary arrays
//static double dot(int nst, double complex *u, double complex *v){
//
//    double complex sum;
//    double norm;
//    int ist;
//
//    sum = 0.0 + 0.0 * I;
//    for(ist = 0; ist < nst; ist++){
//        sum += conj(u[ist]) * v[ist];
//    }
//
//    norm = creal(sum);
//    return norm;
//
//}

// Routine to calculate cdot contribution originated from Ehrenfest term
static void cdot(int pst, double complex **prop_mat, double complex *c, double complex *c_dot){

    int ist, jst;

    for(ist = 0; ist < pst; ist++){
        c_dot[ist] = 0.0 + 0.0 * I;
    }

    for(ist = 0; ist < pst; ist++){
        for(jst = 0; jst < pst; jst++){
            c_dot[ist] += prop_mat[ist][jst] * c[jst];
        }
    }

}
#endif

//// Routine to calculate rhodot contribution originated from Ehrenfest term
//static void rhodot(int nst, double *e, double **dv, double complex **rho, double complex **rho_dot){
//
//    int ist, jst, kst;
//
//    for(ist = 0; ist < nst; ist++){
//        for(jst = 0; jst < nst; jst++){
//            rho_dot[ist][jst] = 0.0 + 0.0 * I;
//        }
//    }
//
//    for(ist = 0; ist < nst; ist++){
//        for(jst = 0; jst < nst; jst++){
//            if(ist != jst){
//                rho_dot[ist][ist] -= dv[ist][jst] * 2.0 * creal(rho[ist][jst]);
//            }
//        }
//    }
//
//    for(ist = 0; ist < nst; ist++){
//        for(jst = ist + 1; jst < nst; jst++){
//            rho_dot[ist][jst] -=  1.0 * I * (e[jst] - e[ist]) * rho[ist][jst];
//            for(kst = 0; kst < nst; kst++){
//                rho_dot[ist][jst] -= dv[ist][kst] * rho[kst][jst] + dv[jst][kst] * rho[ist][kst];
//            }
//            rho_dot[jst][ist] = conj(rho_dot[ist][jst]);
//        }
//    }
//
//}


