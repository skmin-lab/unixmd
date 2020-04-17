#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

//double dot(double complex* u, double complex* v, int nst);
//void cdot(int nst, double complex* c, double* e, double **dv, double complex *c_dot);

/*static void rhodot(int nst, double complex **rho, double *e, double **dv, double complex **rho_dot){
    
    int ist, jst, kst;
    for(ist = 0; ist < nst; ist++){
        for(jst = 0; jst < nst; jst++){
            rho_dot[ist][jst] = 0.0 + 0.0 * I;
        }
    }

    for(ist = 0; ist < nst; ist++){
        for(jst = 0; jst < nst; jst++){
            if(ist != jst){
                rho_dot[ist][ist] -= dv[ist][jst] * 2.0 * creal(rho[ist][jst]);
            }
        }
    }

    for(ist = 0; ist < nst; ist++){
        for(jst = ist + 1; jst < nst; jst++){
            rho_dot[ist][jst] -=  1.0 * I * (e[jst] - e[ist]) * rho[ist][jst];

            for(kst = 0; kst < nst; kst++){
                rho_dot[ist][jst] -= dv[ist][kst] * rho[kst][jst] + dv[jst][kst] * rho[ist][kst];
            }
            rho_dot[jst][ist] = conj(rho_dot[ist][jst]);
        }
    }
}*/

//static void TD_NAC(int nst, int nesteps, double dt, double complex **rho, double *energy, double *energy_old, double **nacme, double **nacme_old){
static void TD_NAC(int nst){
//    double complex **rho_dot = malloc(nst * sizeof(double complex *));
//    double complex **rho_new = malloc(nst * sizeof(double complex *));
//    double complex **k1 = malloc(nst * sizeof(double complex *));
//    double complex **k2 = malloc(nst * sizeof(double complex *));
//    double complex **k3 = malloc(nst * sizeof(double complex *));
//    double complex **k4 = malloc(nst * sizeof(double complex *));
//    double complex **kfunction = malloc(nst * sizeof(double complex *));
//    double complex **variation = malloc(nst * sizeof(double complex *));
//    double *eenergy = malloc(nst * sizeof(double));
//    double **dv = malloc(nst * sizeof(double*));
//    double *na_term = malloc(nst * sizeof(double));
//    int iestep, ist, jst;
//    double frac, edt, norm;

//    printf("TD_NAC : nst = %8d\n", nst);

/*    frac = 1.0 / (double)nesteps;
    edt = dt * frac;

    for(ist = 0; ist < nst; ist++){
        dv[ist] = malloc(nst * sizeof(double));
        rho_dot[ist] = malloc(nst * sizeof(double complex));
        rho_new[ist] = malloc(nst * sizeof(double complex));
        k1[ist] = malloc(nst * sizeof(double complex));
        k2[ist] = malloc(nst * sizeof(double complex));
        k3[ist] = malloc(nst * sizeof(double complex));
        k4[ist] = malloc(nst * sizeof(double complex));
        kfunction[ist] = malloc(nst * sizeof(double complex));
        variation[ist] = malloc(nst * sizeof(double complex));
    }
    
    for(iestep = 0; iestep < nesteps; iestep++){
        
        for(ist = 0; ist < nst; ist++){
            eenergy[ist] = energy_old[ist] + (energy[ist] - energy_old[ist]) *(double) iestep * frac;
            for(jst = 0; jst < nst; jst++){
                dv[ist][jst] = nacme_old[ist][jst] + (nacme[ist][jst] - nacme_old[ist][jst]) *(double) iestep * frac;
            }
        }

        rhodot(nst, rho, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k1[ist][jst] = edt * rho_dot[ist][jst];
                kfunction[ist][jst] = 0.5 * k1[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        rhodot(nst, rho_new, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k2[ist][jst] = edt * rho_dot[ist][jst];
                kfunction[ist][jst] = 0.5 * (-1.0 + sqrt(2.0)) *  k1[ist][jst] + (1.0 - 0.5 * sqrt(2.0)) * k2[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        rhodot(nst, rho_new, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k3[ist][jst] = edt * rho_dot[ist][jst];
                kfunction[ist][jst] = - 0.5 * sqrt(2.0) * k2[ist][jst] + (1.0 + 0.5 * sqrt(2.0)) * k3[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        rhodot(nst, rho_new, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k4[ist][jst] = edt * rho_dot[ist][jst];
                variation[ist][jst] = (k1[ist][jst] + (2.0 - sqrt(2.0)) * k2[ist][jst] + (2.0 + sqrt(2.0)) * k3[ist][jst] + k4[ist][jst]) / 6.0;
                rho[ist][jst] += variation[ist][jst];
            }
        }
    }
    
    for(ist = 0; ist < nst; ist++){
        na_term[ist] = 0.0;
        for(jst = 0; jst < nst; jst++){
            if(jst != ist){
                na_term[ist] -= 2.0 * dv[ist][jst] * creal(rho[ist][jst]);
            }
        }
    }
    norm = 0.0;
    for(ist = 0; ist < nst; ist++){
        norm += creal(rho[ist][ist]);
    }

    free(eenergy);
    free(na_term);
    for(ist = 0; ist < nst; ist++){
        free(dv[ist]);
        free(k1[ist]);
        free(k2[ist]);
        free(k3[ist]);
        free(k4[ist]);
        free(kfunction[ist]);
        free(variation[ist]);
        free(rho_dot[ist]);
        free(rho_new[ist]);
    }
    free(dv);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(kfunction);
    free(variation);
    free(rho_dot);
    free(rho_new);
    */
}


