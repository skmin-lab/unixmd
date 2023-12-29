#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "derivs.h"

// Routine for coefficient propagation scheme in rk4 propagator
static void rk4_coef(int pst, int nesteps, double dt, int **get_d_ind, double **ham_d,
    double **ham_d_old, double **nacme, double **nacme_old, double complex *coef_d);

// Routine for density propagation scheme in rk4 propagator
//static void rk4_rho(int nst, int nesteps, double dt, double *energy, double *energy_old,
//    double **nacme, double **nacme_old, double complex **rho);

// Interface routine for propagation scheme in rk4 propagator
static void rk4(int pst, int nesteps, double dt, char *elec_object, int **get_d_ind, double **ham_d,
    double **ham_d_old, double **nacme, double **nacme_old, double complex *coef_d){

    if(strcmp(elec_object, "coefficient") == 0){
        rk4_coef(pst, nesteps, dt, get_d_ind, ham_d, ham_d_old, nacme, nacme_old, coef_d);
    }
//    else if(strcmp(elec_object, "density") == 0){
//        rk4_rho(nst, nesteps, dt, energy, energy_old, nacme, nacme_old, rho);
//    }

}

// Routine for coefficient propagation scheme in rk4 propagator
static void rk4_coef(int pst, int nesteps, double dt, int **get_d_ind, double **ham_d,
    double **ham_d_old, double **nacme, double **nacme_old, double complex *coef_d){

    double complex *k1 = malloc(pst * sizeof(double complex));
    double complex *k2 = malloc(pst * sizeof(double complex));
    double complex *k3 = malloc(pst * sizeof(double complex));
    double complex *k4 = malloc(pst * sizeof(double complex));
    double complex *kfunction = malloc(pst * sizeof(double complex));
    double complex *variation = malloc(pst * sizeof(double complex));
    double complex *c_dot = malloc(pst * sizeof(double complex));
    double complex *coef_new = malloc(pst * sizeof(double complex));
    double complex **prop_mat = malloc(pst * sizeof(double complex*));

    int ist, jst, iestep, ind_mol1, ind_mol2, ind_photon1, ind_photon2;
    // TODO : Is norm necessary?
    double frac, edt, erel, tmp1, tmp2;//, norm;

    for(ist = 0; ist < pst; ist++){
        prop_mat[ist] = malloc(pst * sizeof(double complex));
    }

    frac = 1.0 / (double)nesteps;
    edt = dt * frac;

    for(iestep = 0; iestep < nesteps; iestep++){

        // Interpolate ham_d and NACME terms between time t and t + dt
        for(ist = 0; ist < pst; ist++){
            ind_mol1 = get_d_ind[ist][0];
            ind_photon1 = get_d_ind[ist][1];
            for(jst = 0; jst < pst; jst++){
                ind_mol2 = get_d_ind[jst][0];
                ind_photon2 = get_d_ind[jst][1];

                tmp1 = ham_d_old[ist][jst] + (ham_d[ist][jst] - ham_d_old[ist][jst])
                    * (double)iestep * frac;
                // Save the lowest energy at current electronic step
                if(ist == 0 && jst == 0){
                    erel = tmp1;
                }
                // To increase the stability of electronic propagation, subtract the relative energy
                if(ist == jst){
                    tmp1 -= erel;
                }
                tmp2 = 0.0;
                if(ind_photon1 == ind_photon2){
                    tmp2 = nacme_old[ind_mol1][ind_mol2] + (nacme[ind_mol1][ind_mol2] - nacme_old[ind_mol1][ind_mol2])
                        * (double)iestep * frac;
                }
                prop_mat[ist][jst] = - 1.0 * tmp1 * I - tmp2;

            }
        }

        // Calculate k1
        cdot(pst, prop_mat, coef_d, c_dot);

        for(ist = 0; ist < pst; ist++){
            k1[ist] = edt * c_dot[ist];
            kfunction[ist] = 0.5 * k1[ist];
            coef_new[ist] = coef_d[ist] + kfunction[ist];
        }

        // Calculate k2
        cdot(pst, prop_mat, coef_new, c_dot);

        for(ist = 0; ist < pst; ist++){
            k2[ist] = edt * c_dot[ist];
            kfunction[ist] = 0.5 * (- 1.0 + sqrt(2.0)) * k1[ist] + (1.0 - 0.5 * sqrt(2.0)) * k2[ist];
            coef_new[ist] = coef_d[ist] + kfunction[ist];
        }

        // Calculate k3
        cdot(pst, prop_mat, coef_new, c_dot);

        for(ist = 0; ist < pst; ist++){
            k3[ist] = edt * c_dot[ist];
            kfunction[ist] = - 0.5 * sqrt(2.0) * k2[ist] + (1.0 + 0.5 * sqrt(2.0)) * k3[ist];
            coef_new[ist] = coef_d[ist] + kfunction[ist];
        }

        // Calculate k4
        cdot(pst, prop_mat, coef_new, c_dot);

        for(ist = 0; ist < pst; ist++){
            k4[ist] = edt * c_dot[ist];
            variation[ist] = (k1[ist] + (2.0 - sqrt(2.0)) * k2[ist] + (2.0 + sqrt(2.0))
                * k3[ist] + k4[ist]) / 6.0;
            coef_new[ist] = coef_d[ist] + variation[ist];
        }

        // TODO : Is this part necessary?
        // Renormalize the coefficients
//        norm = dot(nst, coef_new, coef_new);
        for(ist = 0; ist < pst; ist++){
//            coef_new[ist] /= sqrt(norm);
            coef_d[ist] = coef_new[ist];
        }

    }

    for(ist = 0; ist < pst; ist++){
        free(prop_mat[ist]);
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(kfunction);
    free(variation);
    free(c_dot);
    free(coef_new);
    free(prop_mat);

}

//// Routine for density propagation scheme in rk4 propagator
//static void rk4_rho(int nst, int nesteps, double dt, double *energy, double *energy_old,
//    double **nacme, double **nacme_old, double complex **rho){
//
//    double complex **k1 = malloc(nst * sizeof(double complex*));
//    double complex **k2 = malloc(nst * sizeof(double complex*));
//    double complex **k3 = malloc(nst * sizeof(double complex*));
//    double complex **k4 = malloc(nst * sizeof(double complex*));
//    double complex **kfunction = malloc(nst * sizeof(double complex*));
//    double complex **variation = malloc(nst * sizeof(double complex*));
//    double complex **rho_dot = malloc(nst * sizeof(double complex*));
//    double complex **rho_new = malloc(nst * sizeof(double complex*));
//    double *eenergy = malloc(nst * sizeof(double));
////    double *na_term = malloc(nst * sizeof(double));
//    double **dv = malloc(nst * sizeof(double*));
//
//    int ist, jst, iestep;
//    double frac, edt;
//
//    for(ist = 0; ist < nst; ist++){
//        k1[ist] = malloc(nst * sizeof(double complex));
//        k2[ist] = malloc(nst * sizeof(double complex));
//        k3[ist] = malloc(nst * sizeof(double complex));
//        k4[ist] = malloc(nst * sizeof(double complex));
//        kfunction[ist] = malloc(nst * sizeof(double complex));
//        variation[ist] = malloc(nst * sizeof(double complex));
//        rho_dot[ist] = malloc(nst * sizeof(double complex));
//        rho_new[ist] = malloc(nst * sizeof(double complex));
//        dv[ist] = malloc(nst * sizeof(double));
//    }
//
//    frac = 1.0 / (double)nesteps;
//    edt = dt * frac;
//
//    for(iestep = 0; iestep < nesteps; iestep++){
//
//        // Interpolate energy and NACME terms between time t and t + dt
//        for(ist = 0; ist < nst; ist++){
//            eenergy[ist] = energy_old[ist] + (energy[ist] - energy_old[ist]) * (double)iestep * frac;
//            for(jst = 0; jst < nst; jst++){
//                dv[ist][jst] = nacme_old[ist][jst] + (nacme[ist][jst] - nacme_old[ist][jst])
//                    * (double)iestep * frac;
//            }
//        }
//
//        // Calculate k1
//        rhodot(nst, eenergy, dv, rho, rho_dot);
//
//        for(ist = 0; ist < nst; ist++){
//            for(jst = 0; jst < nst; jst++){
//                k1[ist][jst] = edt * rho_dot[ist][jst];
//                kfunction[ist][jst] = 0.5 * k1[ist][jst];
//                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
//            }
//        }
//
//        // Calculate k2
//        rhodot(nst, eenergy, dv, rho_new, rho_dot);
//
//        for(ist = 0; ist < nst; ist++){
//            for(jst = 0; jst < nst; jst++){
//                k2[ist][jst] = edt * rho_dot[ist][jst];
//                kfunction[ist][jst] = 0.5 * (- 1.0 + sqrt(2.0)) * k1[ist][jst]
//                    + (1.0 - 0.5 * sqrt(2.0)) * k2[ist][jst];
//                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
//            }
//        }
//
//        // Calculate k3
//        rhodot(nst, eenergy, dv, rho_new, rho_dot);
//
//        for(ist = 0; ist < nst; ist++){
//            for(jst = 0; jst < nst; jst++){
//                k3[ist][jst] = edt * rho_dot[ist][jst];
//                kfunction[ist][jst] = - 0.5 * sqrt(2.0) * k2[ist][jst]
//                    + (1.0 + 0.5 * sqrt(2.0)) * k3[ist][jst];
//                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
//            }
//        }
//
//        // Calculate k4
//        rhodot(nst, eenergy, dv, rho_new, rho_dot);
//
//        for(ist = 0; ist < nst; ist++){
//            for(jst = 0; jst < nst; jst++){
//                k4[ist][jst] = edt * rho_dot[ist][jst];
//                variation[ist][jst] = (k1[ist][jst] + (2.0 - sqrt(2.0)) * k2[ist][jst]
//                    + (2.0 + sqrt(2.0)) * k3[ist][jst] + k4[ist][jst]) / 6.0;
//                rho[ist][jst] += variation[ist][jst];
//            }
//        }
//
//    }
//
//    /*
//    for(ist = 0; ist < nst; ist++){
//        na_term[ist] = 0.0;
//        for(jst = 0; jst < nst; jst++){
//            if(jst != ist){
//                na_term[ist] -= 2.0 * dv[ist][jst] * creal(rho[ist][jst]);
//            }
//        }
//        printf("RHODOT_NAC %d %f\n",ist+1, na_term[ist]);
//    }
//    norm = 0.0;
//    for(ist = 0; ist < nst; ist++){
//        norm += creal(rho[ist][ist]);
//    }
//    printf("RK4_COEF : NORM = %15.8f\n", norm);
//    */
//
//    for(ist = 0; ist < nst; ist++){
//        free(k1[ist]);
//        free(k2[ist]);
//        free(k3[ist]);
//        free(k4[ist]);
//        free(kfunction[ist]);
//        free(variation[ist]);
//        free(rho_dot[ist]);
//        free(rho_new[ist]);
//        free(dv[ist]);
//    }
//
//    free(k1);
//    free(k2);
//    free(k3);
//    free(k4);
//    free(kfunction);
//    free(variation);
//    free(rho_dot);
//    free(rho_new);
//    free(eenergy);
////    free(na_term);
//    free(dv);
//
//}


