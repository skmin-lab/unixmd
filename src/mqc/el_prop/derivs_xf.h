#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

// Routine to calculate cdot contribution originated from XF term
static void xf_cdot(int nat, int nsp, int nst, int *l_coh, double *mass, double *wsigma,
    double **pos, double ***aux_pos, double ***phase, double complex *c, double complex *xfcdot){

    double **qmom = malloc(nat * sizeof(double*));
    double **dec = malloc(nst * sizeof(double*));
    double *rho = malloc(nst * sizeof(double));

    int ist, jst, iat, isp;

    // Initialize variables related to decoherence
    for(iat = 0; iat < nat; iat++){
        qmom[iat] = malloc(nsp * sizeof(double));
        for(isp = 0; isp < nsp; isp++){
            qmom[iat][isp] = 0.0;
        }
    }

    for(ist = 0; ist < nst; ist++){
        dec[ist] = malloc(nst * sizeof(double));
        for(jst = 0; jst < nst; jst++){
            dec[ist][jst] = 0.0;
        }
    }

    // Calculate densities from current coefficients
    for(ist = 0; ist < nst; ist++){
        rho[ist] = creal(conj(c[ist]) * c[ist]);
    }

    // Get quantum momentum from auxiliary positions and sigma values
    for(ist = 0; ist < nst; ist++){

        if(l_coh[ist] == 1){
            for(iat = 0; iat < nat; iat++){
                for(isp = 0; isp < nsp; isp++){
                    qmom[iat][isp] += 0.5 * rho[ist] * (pos[iat][isp] - aux_pos[ist][iat][isp])
                        / pow(wsigma[iat],2.0) / mass[iat];
                }
            }
        }

    }

    // Get decoherence term from quantum momentum and phase
    for(ist = 0; ist < nst; ist++){
        for(jst = ist + 1; jst < nst; jst++){

            if(l_coh[ist] == 1 && l_coh[jst] == 1){
                for(iat = 0; iat < nat; iat++){
                    for(isp = 0; isp < nsp; isp++){
                        dec[ist][jst] += qmom[iat][isp] * (phase[ist][iat][isp] - phase[jst][iat][isp]);
                    }
                }
            }
            dec[jst][ist] = - 1.0 * dec[ist][jst];

        }
    }

    // Get cdot contribution from decoherence term
    for(ist = 0; ist < nst; ist++){
        xfcdot[ist] = 0.0 + 0.0 * I;
        for(jst = 0; jst < nst; jst++){
            xfcdot[ist] -= rho[jst] * dec[jst][ist] * c[ist];
        }
    }

    // Deallocate temporary arrays
    for(iat = 0; iat < nat; iat++){
        free(qmom[iat]);
    }
    for(ist = 0; ist < nst; ist++){
        free(dec[ist]);
    }

    free(qmom);
    free(dec);
    free(rho);

}

// Routine to calculate rhodot contribution originated from XF term
static void xf_rhodot(int nat, int nsp, int nst, int *l_coh, double *mass, double *wsigma,
    double **pos, double ***aux_pos, double ***phase, double complex **rho, double complex **xfrhodot){

    double **qmom = malloc(nat * sizeof(double*));
    double **dec = malloc(nst * sizeof(double*));

    int ist, jst, kst, iat, isp;

    // Initialize variables related to decoherence
    for(iat = 0; iat < nat; iat++){
        qmom[iat] = malloc(nsp * sizeof(double));
        for(isp = 0; isp < nsp; isp++){
            qmom[iat][isp] = 0.0;
        }
    }

    for(ist = 0; ist < nst; ist++){
        dec[ist] = malloc(nst * sizeof(double));
        for(jst = 0; jst < nst; jst++){
            dec[ist][jst] = 0.0;
        }
    }

    // Get quantum momentum from auxiliary positions and sigma values
    for(ist = 0; ist < nst; ist++){

        if(l_coh[ist] == 1){
            for(iat = 0; iat < nat; iat++){
                for(isp = 0; isp < nsp; isp++){
                    qmom[iat][isp] += 0.5 * creal(rho[ist][ist]) * (pos[iat][isp] - aux_pos[ist][iat][isp])
                        / pow(wsigma[iat],2.0) / mass[iat];
                }
            }
        }

    }

    // Get decoherence term from quantum momentum and phase
    for(ist = 0; ist < nst; ist++){
        for(jst = ist + 1; jst < nst; jst++){

            if(l_coh[ist] == 1 && l_coh[jst] == 1){
                for(iat = 0; iat < nat; iat++){
                    for(isp = 0; isp < nsp; isp++){
                        dec[ist][jst] += qmom[iat][isp] * (phase[ist][iat][isp] - phase[jst][iat][isp]);
                    }
                }
            }
            dec[jst][ist] = - 1.0 * dec[ist][jst];

        }
    }

    // Get rhodot contribution from decoherence term
    for(ist = 0; ist < nst; ist++){
        // Diagonal components
        xfrhodot[ist][ist] = 0.0 + 0.0 * I;
        for(kst = 0; kst < nst; kst++){
            xfrhodot[ist][ist] -= 2.0 * dec[kst][ist] * rho[ist][kst] * rho[kst][ist];
        }
        // Off-diagonal components
        for(jst = ist + 1; jst < nst; jst++){
            xfrhodot[ist][jst] = 0.0 + 0.0 * I;
            for(kst = 0; kst < nst; kst++){
                xfrhodot[ist][jst] -= (dec[kst][ist] + dec[kst][jst]) * rho[ist][kst] * rho[kst][jst];
            }
            xfrhodot[jst][ist] = conj(xfrhodot[ist][jst]);
        }
    }

    // Deallocate temporary arrays
    for(iat = 0; iat < nat; iat++){
        free(qmom[iat]);
    }
    for(ist = 0; ist < nst; ist++){
        free(dec[ist]);
    }

    free(qmom);
    free(dec);

}


