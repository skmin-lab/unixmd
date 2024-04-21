#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "derivs.h"
#include "derivs_xf.h"
#include "transform.h"

// Routine for coefficient propagation scheme in rk4 propagator
static void rk4_coef(int nat, int ndim, int pst, int nesteps, int verbosity, double dt,
    int *l_coh, double *mass, double *sigma, int **get_d_ind, double **unitary, double **ham_d,
    double **ham_d_old, double **nacme, double **nacme_old, double **pos, double ***aux_pos,
    double ***phase, double *dotpopdec_d, double complex *coef_d, double **qmom);

// Routine for density propagation scheme in rk4 propagator
//static void rk4_rho(int nat, int ndim, int nst, int nesteps, double dt, int *l_coh,
//    double *mass, double *energy, double *energy_old, double *sigma, double **nacme,
//    double **nacme_old, double **pos, double **qmom, double ***aux_pos, double ***phase, double complex **rho,
//    int verbosity, double *dotpopdec);

// Interface routine for propagation scheme in rk4 propagator
static void rk4(int nat, int ndim, int pst, int nesteps, int verbosity, double dt, char *elec_object,
    int *l_coh, double *mass, double *sigma, int **get_d_ind, double **unitary, double **ham_d,
    double **ham_d_old, double **nacme, double **nacme_old, double **pos, double ***aux_pos,
    double ***phase, double *dotpopdec_d, double complex *coef_d, double **qmom){

    if(strcmp(elec_object, "coefficient") == 0){
        rk4_coef(nat, ndim, pst, nesteps, verbosity, dt, l_coh, mass, sigma, get_d_ind, unitary,
            ham_d, ham_d_old, nacme, nacme_old, pos, aux_pos, phase, dotpopdec_d, coef_d, qmom);
    }
//    else if(strcmp(elec_object, "density") == 0){
//        rk4_rho(nat, ndim, nst, nesteps, dt, l_coh, mass, energy, energy_old, sigma,
//            nacme, nacme_old, pos, qmom, aux_pos, phase, rho, verbosity, dotpopdec);
//    }

}

// Routine for coefficient propagation scheme in rk4 propagator
static void rk4_coef(int nat, int ndim, int pst, int nesteps, int verbosity, double dt,
    int *l_coh, double *mass, double *sigma, int **get_d_ind, double **unitary, double **ham_d,
    double **ham_d_old, double **nacme, double **nacme_old, double **pos, double ***aux_pos,
    double ***phase, double *dotpopdec_d, double complex *coef_d, double **qmom){

    double complex *k1 = malloc(pst * sizeof(double complex));
    double complex *k2 = malloc(pst * sizeof(double complex));
    double complex *k3 = malloc(pst * sizeof(double complex));
    double complex *k4 = malloc(pst * sizeof(double complex));
    double complex *kfunction = malloc(pst * sizeof(double complex));
    double complex *variation = malloc(pst * sizeof(double complex));
    double complex *c_dot = malloc(pst * sizeof(double complex));
    double complex *xf_c_dot = malloc(pst * sizeof(double complex));
    double complex *coef_a = malloc(pst * sizeof(double complex));
    double complex *coef_new = malloc(pst * sizeof(double complex));
    double complex **prop_mat_d = malloc(pst * sizeof(double complex*));
    double **dec_mat = malloc(pst * sizeof(double*));

    int ist, jst, iestep, ind_mol1, ind_mol2, ind_photon1, ind_photon2;
    int iat, isp;
    // TODO : Is norm necessary?
    double frac, edt, erel, tmp1, tmp2;//, norm;
    double rho;

    for(ist = 0; ist < pst; ist++){
        prop_mat_d[ist] = malloc(pst * sizeof(double complex));
    }

    for(ist = 0; ist < pst; ist++){
        dec_mat[ist] = malloc(pst * sizeof(double));
    }

    frac = 1.0 / (double)nesteps;
    edt = dt * frac;

    for(iestep = 0; iestep < nesteps; iestep++){

        // Get polaritonic state coefficients for calculation of decoherence term
        transform_d2a(pst, unitary, coef_d, coef_a);

        // Get quantum momentum from auxiliary positions and sigma values
        for(iat = 0; iat < nat; iat++){
            for(isp = 0; isp < ndim; isp++){
                qmom[iat][isp] = 0.0;
            }
        }
        for(ist = 0; ist < pst; ist++){

            if(l_coh[ist] == 1){
                rho = creal(conj(coef_a[ist]) * coef_a[ist]);
                for(iat = 0; iat < nat; iat++){
                    for(isp = 0; isp < ndim; isp++){
                        qmom[iat][isp] += 0.5 * rho * (pos[iat][isp] - aux_pos[ist][iat][isp])
                            / pow(sigma[iat], 2.0) / mass[iat];
                    }
                }
            }

        }

        // Get decoherence term from quantum momentum and phase
        for(ist = 0; ist < pst; ist++){
            for(jst = 0; jst < pst; jst++){
                dec_mat[ist][jst] = 0.0;
            }
        }
        for(ist = 0; ist < pst; ist++){
            for(jst = ist + 1; jst < pst; jst++){

                if(l_coh[ist] == 1 && l_coh[jst] == 1){
                    for(iat = 0; iat < nat; iat++){
                        for(isp = 0; isp < ndim; isp++){
                            dec_mat[ist][jst] += qmom[iat][isp] * (phase[ist][iat][isp] - phase[jst][iat][isp]);
                        }
                    }
                }
                dec_mat[jst][ist] = - 1.0 * dec_mat[ist][jst];

            }
        }

        // Calculate cdot contribution originating from XF term
        xf_cdot(pst, unitary, dec_mat, coef_a, xf_c_dot);

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
                prop_mat_d[ist][jst] = - 1.0 * tmp1 * I - tmp2;

            }
        }

        // Calculate k1
        cdot(pst, prop_mat_d, coef_d, c_dot);

        for(ist = 0; ist < pst; ist++){
            k1[ist] = edt * (c_dot[ist] + xf_c_dot[ist]);
            kfunction[ist] = 0.5 * k1[ist];
            coef_new[ist] = coef_d[ist] + kfunction[ist];
        }

        // Calculate k2
        cdot(pst, prop_mat_d, coef_new, c_dot);

        for(ist = 0; ist < pst; ist++){
            k2[ist] = edt * (c_dot[ist] + xf_c_dot[ist]);
            kfunction[ist] = 0.5 * (- 1.0 + sqrt(2.0)) * k1[ist] + (1.0 - 0.5 * sqrt(2.0)) * k2[ist];
            coef_new[ist] = coef_d[ist] + kfunction[ist];
        }

        // Calculate k3
        cdot(pst, prop_mat_d, coef_new, c_dot);

        for(ist = 0; ist < pst; ist++){
            k3[ist] = edt * (c_dot[ist] + xf_c_dot[ist]);
            kfunction[ist] = - 0.5 * sqrt(2.0) * k2[ist] + (1.0 + 0.5 * sqrt(2.0)) * k3[ist];
            coef_new[ist] = coef_d[ist] + kfunction[ist];
        }

        // Calculate k4
        cdot(pst, prop_mat_d, coef_new, c_dot);

        for(ist = 0; ist < pst; ist++){
            k4[ist] = edt * (c_dot[ist] + xf_c_dot[ist]);
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

    if(verbosity >= 1){
        transform_d2a(pst, unitary, coef_d, coef_a);
        xf_print_coef(pst, unitary, dec_mat, coef_d, coef_a, dotpopdec_d);
    }

    for(ist = 0; ist < pst; ist++){
        free(prop_mat_d[ist]);
    }
    for(ist = 0; ist < pst; ist++){
        free(dec_mat[ist]);
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(kfunction);
    free(variation);
    free(c_dot);
    free(xf_c_dot);
    free(coef_a);
    free(coef_new);
    free(prop_mat_d);
    free(dec_mat);

}

//// Routine for density propagation scheme in rk4 propagator
//static void rk4_rho(int nat, int ndim, int nst, int nesteps, double dt, int *l_coh,
//    double *mass, double *energy, double *energy_old, double *sigma, double **nacme,
//    double **nacme_old, double **pos, double **qmom, double ***aux_pos, double ***phase, double complex **rho,
//    int verbosity, double *dotpopdec){
//
//    double complex **k1 = malloc(nst * sizeof(double complex*));
//    double complex **k2 = malloc(nst * sizeof(double complex*));
//    double complex **k3 = malloc(nst * sizeof(double complex*));
//    double complex **k4 = malloc(nst * sizeof(double complex*));
//    double complex **kfunction = malloc(nst * sizeof(double complex*));
//    double complex **variation = malloc(nst * sizeof(double complex*));
//    double complex **rho_dot = malloc(nst * sizeof(double complex*));
//    double complex **xf_rho_dot = malloc(nst * sizeof(double complex*));
//    double complex **rho_new = malloc(nst * sizeof(double complex*));
//    double *eenergy = malloc(nst * sizeof(double));
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
//        xf_rho_dot[ist] = malloc(nst * sizeof(double complex));
//        rho_new[ist] = malloc(nst * sizeof(double complex));
//        dv[ist] = malloc(nst * sizeof(double));
//    }
//
//    frac = 1.0 / (double)nesteps;
//    edt = dt * frac;
//
//    for(iestep = 0; iestep < nesteps; iestep++){
//
//        // Calculate rhodot contribution originated from XF term
//        xf_rhodot(nat, ndim, nst, l_coh, mass, sigma, pos, qmom, aux_pos, phase, rho, xf_rho_dot);
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
//                k1[ist][jst] = edt * (rho_dot[ist][jst] + xf_rho_dot[ist][jst]);
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
//                k2[ist][jst] = edt * (rho_dot[ist][jst] + xf_rho_dot[ist][jst]);
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
//                k3[ist][jst] = edt * (rho_dot[ist][jst] + xf_rho_dot[ist][jst]);
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
//                k4[ist][jst] = edt * (rho_dot[ist][jst] + xf_rho_dot[ist][jst]);
//                variation[ist][jst] = (k1[ist][jst] + (2.0 - sqrt(2.0)) * k2[ist][jst]
//                    + (2.0 + sqrt(2.0)) * k3[ist][jst] + k4[ist][jst]) / 6.0;
//                rho[ist][jst] += variation[ist][jst];
//            }
//        }
//
//    }
//
//    if(verbosity >= 1){
//        xf_print_rho(nst, xf_rho_dot, dotpopdec); 
//    }
//
//    for(ist = 0; ist < nst; ist++){
//        free(k1[ist]);
//        free(k2[ist]);
//        free(k3[ist]);
//        free(k4[ist]);
//        free(kfunction[ist]);
//        free(variation[ist]);
//        free(rho_dot[ist]);
//        free(xf_rho_dot[ist]);
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
//    free(xf_rho_dot);
//    free(rho_new);
//    free(eenergy);
//    free(dv);
//
//}


