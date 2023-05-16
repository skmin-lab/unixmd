#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

// Routine to calculate cdot contribution originating from XF term
static void xf_cdot(int pst, double **unitary, double **dec_mat, double complex *c, double complex *xfcdot){

    double *rho = malloc(pst * sizeof(double));

    int ast, ist, jst;

    // Calculate densities from current coefficients
    for(ist = 0; ist < pst; ist++){
        rho[ist] = creal(conj(c[ist]) * c[ist]);
    }

    for(ast = 0; ast < pst; ast++){
        xfcdot[ast] = 0.0 + 0.0 * I;
    }

    for(ast = 0; ast < pst; ast++){

        for(ist = 0; ist < pst; ist++){
            for(jst = 0; jst < pst; jst++){
                xfcdot[ast] -= unitary[ast][ist] * rho[jst] * dec_mat[jst][ist] * c[ist];
            }
        }

    }

    free(rho);

}

// Routine to print xf debug info 
static void xf_print_coef(int pst, double **unitary, double **dec_mat, double complex *coef_d,
    double complex *coef_a, double *dotpopdec){

    double *rho = malloc(pst * sizeof(double));

    int ast, ist, jst;

    for(ist = 0; ist < pst; ist++){
        rho[ist] = creal(conj(coef_a[ist]) * coef_a[ist]);
    }

    for(ast = 0; ast < pst; ast++){
        dotpopdec[ast] = 0.0;
    }

    for(ast = 0; ast < pst; ast++){

        for(ist = 0; ist < pst; ist++){
            for(jst = 0; jst < pst; jst++){
                dotpopdec[ast] -= unitary[ast][ist] * rho[jst] * dec_mat[jst][ist]
                    * creal( conj(coef_a[ist]) * coef_d[ast] + coef_a[ist] * conj(coef_d[ast]) );
            }
        }

    }

    free(rho);

}

//// Routine to calculate rhodot contribution originated from XF term
//static void xf_rhodot(int nat, int ndim, int nst, int *l_coh, double *mass, double *sigma,
//    double **pos, double **qmom, double ***aux_pos, double ***phase, double complex **rho, double complex **xfrhodot){
//
//    double **dec = malloc(nst * sizeof(double*));
//
//    int ist, jst, kst, iat, isp;
//
//    // Initialize variables related to decoherence
//    for(iat = 0; iat < nat; iat++){
//        for(isp = 0; isp < ndim; isp++){
//            qmom[iat][isp] = 0.0;
//        }
//    }
//
//    for(ist = 0; ist < nst; ist++){
//        dec[ist] = malloc(nst * sizeof(double));
//        for(jst = 0; jst < nst; jst++){
//            dec[ist][jst] = 0.0;
//        }
//    }
//
//    // Get quantum momentum from auxiliary positions and sigma values
//    for(ist = 0; ist < nst; ist++){
//
//        if(l_coh[ist] == 1){
//            for(iat = 0; iat < nat; iat++){
//                for(isp = 0; isp < ndim; isp++){
//                    qmom[iat][isp] += 0.5 * creal(rho[ist][ist]) * (pos[iat][isp] - aux_pos[ist][iat][isp])
//                        / pow(sigma[iat], 2.0) / mass[iat];
//                }
//            }
//        }
//
//    }
//
//    // Get decoherence term from quantum momentum and phase
//    for(ist = 0; ist < nst; ist++){
//        for(jst = ist + 1; jst < nst; jst++){
//
//            if(l_coh[ist] == 1 && l_coh[jst] == 1){
//                for(iat = 0; iat < nat; iat++){
//                    for(isp = 0; isp < ndim; isp++){
//                        dec[ist][jst] += qmom[iat][isp] * (phase[ist][iat][isp] - phase[jst][iat][isp]);
//                    }
//                }
//            }
//            dec[jst][ist] = - 1.0 * dec[ist][jst];
//
//        }
//    }
//
//    // Get rhodot contribution from decoherence term
//    for(ist = 0; ist < nst; ist++){
//        // Diagonal components
//        xfrhodot[ist][ist] = 0.0 + 0.0 * I;
//        for(kst = 0; kst < nst; kst++){
//            xfrhodot[ist][ist] -= 2.0 * dec[kst][ist] * rho[ist][kst] * rho[kst][ist];
//        }
//        // Off-diagonal components
//        for(jst = ist + 1; jst < nst; jst++){
//            xfrhodot[ist][jst] = 0.0 + 0.0 * I;
//            for(kst = 0; kst < nst; kst++){
//                xfrhodot[ist][jst] -= (dec[kst][ist] + dec[kst][jst]) * rho[ist][kst] * rho[kst][jst];
//            }
//            xfrhodot[jst][ist] = conj(xfrhodot[ist][jst]);
//        }
//    }
//
//    // Deallocate temporary arrays
//    for(ist = 0; ist < nst; ist++){
//        free(dec[ist]);
//    }
//
//    free(dec);
//
//}

//// Routine to print xf debug info 
//static void xf_print_rho(int nst, double complex **xfrhodot, double *dotpopdec){
//    int ist;
//    
//    for(ist = 0; ist < nst; ist++){
//        dotpopdec[ist] = creal(xfrhodot[ist][ist]);
//    }
//}


