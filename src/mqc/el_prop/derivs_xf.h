#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

static void xf_cdot(int nst, int nat, int nsp, int *l_coh, double *wsigma, 
                        double *mass, double complex *c, double **pos, double ***aux_pos, double ***phase, 
                        double complex *xfcdot){
    
    int ist, jst, iat, isp;
    double **qmom = malloc(nat * sizeof(double*));
    double **dec = malloc(nst * sizeof(double*));
    double *rho = malloc(nst * sizeof(double)); 
    
    for(iat = 0; iat < nat; iat++){
        qmom[iat] = malloc(nsp * sizeof(double));
        for(isp = 0; isp < nsp; isp++){
            qmom[iat][isp] = 0.0;
        }
    }
    for(ist = 0; ist < nst; ist++){
        dec[ist] = malloc(nst * sizeof(double));
        rho[ist] = creal(conj(c[ist]) * c[ist]);
    }
    
    // get qmom
    for(ist = 0; ist < nst; ist++){
        if(l_coh[ist] == 1){
            for(iat = 0; iat < nat; iat++){
                for(isp = 0; isp < nsp; isp++){
                    qmom[iat][isp] += 0.5 * rho[ist] * (pos[iat][isp] - aux_pos[ist][iat][isp]) / pow(wsigma[iat],2.0) / mass[iat];
                }
            }
        }
    }
    for(ist = 0; ist < nst; ist++){
        dec[ist][ist] = 0.0;
        for(jst = ist + 1; jst < nst; jst++){
            dec[ist][jst] = 0.0;
            if(l_coh[ist] == 1 && l_coh[jst] == 1){
                for(iat = 0; iat < nat; iat++){
                    for(isp = 0; isp < nsp; isp++){
                        dec[ist][jst] += qmom[iat][isp]*(phase[ist][iat][isp]-phase[jst][iat][isp]);
                    }
                }
            }
            dec[jst][ist] = -1.0 * dec[ist][jst];
        }
    }
    for(ist = 0; ist < nst; ist++){
        xfcdot[ist] = 0.0 + 0.0 * I;
        for(jst = 0; jst < nst; jst++){
            xfcdot[ist] -= rho[jst] * dec[jst][ist] * c[ist];
        }
    }
    free(rho);
    for(ist = 0; ist < nst; ist++){
        free(dec[ist]);
    }
    for(iat = 0; iat < nat; iat++){
        free(qmom[iat]);
    }
    free(dec);
    free(qmom);
}

static void xf_rhodot(int nst, int nat, int nsp, int *l_coh, double *wsigma, 
                      double *mass, double complex **rho, double **pos, double ***aux_pos, double ***phase, 
                      double complex **xfrhodot){
    int ist, jst, kst, iat, isp;
    double **qmom = malloc(nat * sizeof(double*));
    double **dec = malloc(nst * sizeof(double*));
    
    for(iat = 0; iat < nat; iat++){
        qmom[iat] = malloc(nsp * sizeof(double));
        for(isp = 0; isp < nsp; isp++){
            qmom[iat][isp] = 0.0;
        }
    }
    for(ist = 0; ist < nst; ist++){
        dec[ist] = malloc(nst * sizeof(double));
    }
    
    // get qmom
    for(ist = 0; ist < nst; ist++){
        if(l_coh[ist] == 1){
            for(iat = 0; iat < nat; iat++){
                for(isp = 0; isp < nsp; isp++){
                    qmom[iat][isp] += 0.5 * creal(rho[ist][ist]) * (pos[iat][isp] - aux_pos[ist][iat][isp]) / pow(wsigma[iat],2.0) / mass[iat];
                }
            }
        }
    }
    for(ist = 0; ist < nst; ist++){
        dec[ist][ist] = 0.0;
        for(jst = ist + 1; jst < nst; jst++){
            dec[ist][jst] = 0.0;
            if(l_coh[ist] == 1 && l_coh[jst] == 1){
                for(iat = 0; iat < nat; iat++){
                    for(isp = 0; isp < nsp; isp++){
                        dec[ist][jst] += qmom[iat][isp]*(phase[ist][iat][isp]-phase[jst][iat][isp]);
                    }
                }
            }
            dec[jst][ist] = -1.0 * dec[ist][jst];
        }
    }
    for(ist = 0; ist < nst; ist++){
        xfrhodot[ist][ist] = 0.0 + 0.0 * I; 
        for(kst = 0; kst < nst; kst++){
            xfrhodot[ist][ist] -= 2.0 * dec[kst][ist] * rho[ist][kst] * rho[kst][ist];
        }
        for(jst = ist + 1; jst < nst; jst++){
            xfrhodot[ist][jst] = 0.0 + 0.0 * I;
            for(kst = 0; kst < nst; kst++){
                xfrhodot[ist][jst] -= (dec[kst][ist]+dec[kst][jst]) * rho[ist][kst] * rho[kst][jst];
            }
            xfrhodot[jst][ist] = conj(xfrhodot[ist][jst]);
        }
    }
    for(ist = 0; ist < nst; ist++){
        free(dec[ist]);
    }
    for(iat = 0; iat < nat; iat++){
        free(qmom[iat]);
    }
    free(dec);
    free(qmom);
}
