#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "derivs.h"

// Complex datatype
struct _dcomplex {double real, imag;};
typedef struct _dcomplex dcomplex;

// Import external heev and gemm routines
extern void zheev_(char *jobz, char *uplo, int *n, dcomplex *a, int *lda, double *w,
    dcomplex *work, int *lwork, double *rwork, int *info);
extern void zgemm_(char *transa, char *transb, int *m, int *n, int *k, dcomplex *alpha,
    dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *beta, dcomplex *c, int *ldc);

// Routine for coefficient propagation scheme in exponential propagator
static void exponential_coef(int pst, int nesteps, double dt, int **get_d_ind, double **ham_d,
    double **ham_d_old, double **nacme, double **nacme_old, double complex *coef_d);

// Routine for density propagation scheme in rk4 propagator
//static void rk4_rho(int nat, int ndim, int nst, int nesteps, double dt, int *l_coh,
//    double *mass, double *energy, double *energy_old, double *sigma, double **nacme,
//    double **nacme_old, double **pos, double **qmom, double ***aux_pos, double ***phase, double complex **rho,
//    int verbosity, double *dotpopdec);

// Interface routine for propagation scheme in exponential propagator
static void exponential(int pst, int nesteps, double dt, char *elec_object, int **get_d_ind, double **ham_d,
    double **ham_d_old, double **nacme, double **nacme_old, double complex *coef_d){

    if(strcmp(elec_object, "coefficient") == 0){
        exponential_coef(pst, nesteps, dt, get_d_ind, ham_d, ham_d_old, nacme, nacme_old, coef_d);
    }
//    else if(strcmp(elec_object, "density") == 0){
//        rk4_rho(nat, ndim, nst, nesteps, dt, l_coh, mass, energy, energy_old, sigma,
//            nacme, nacme_old, pos, qmom, aux_pos, phase, rho, verbosity, dotpopdec);
//    }

}

// Routine for coefficient propagation scheme in exponential propagator
static void exponential_coef(int pst, int nesteps, double dt, int **get_d_ind, double **ham_d,
    double **ham_d_old, double **nacme, double **nacme_old, double complex *coef_d){

    double complex *coef_new = malloc(pst * sizeof(double complex));
    double complex **prop_mat_d = malloc(pst * sizeof(double complex*));

    // (Hamiltonian - i * (NACME + decoherence)) * dt
    double complex **exponent = malloc((pst) * sizeof(double complex*));
    // eigenvectors of (Hamiltonian - i * (NACME + decoherence)) * dt, P
    dcomplex *eigenvectors = malloc((pst * pst) * sizeof(dcomplex));
    // eigenvalues of (Hamiltonian - i * (NACME + decoherence)) * dt, D
    double *eigenvalues = malloc(pst * sizeof(double));

    // diagonal matrix using eigenvalues, exp(- i * D)
    dcomplex *exp_idiag = malloc((pst * pst) * sizeof(dcomplex));
    // exp(- i * exponent) = P * exp(- i * D) * P^-1
    dcomplex *exp_iexponent = malloc((pst * pst) * sizeof(dcomplex));
    // product of (P * exp(- i * D) * P^-1) until previous step
    dcomplex *product_old = malloc((pst * pst) * sizeof(dcomplex));
    // product of (P * exp(- i * D) * P^-1) until current step
    dcomplex *product_new = malloc((pst * pst) * sizeof(dcomplex));
    // final product of exp(- i * exponent)
    double complex **propagator = malloc((pst) * sizeof(double complex*));

    dcomplex *tmp_mat = malloc((pst * pst) * sizeof(dcomplex));
    dcomplex *identity = malloc((pst * pst) * sizeof(dcomplex));

    dcomplex wkopt;
    dcomplex *work;
    double *rwork = malloc((3 * pst - 2) * sizeof(double));
    int lwork, info;

    dcomplex dcone = {1.0, 0.0};
    dcomplex dczero = {0.0, 0.0};

    int ist, jst, kst, iestep, ind_mol1, ind_mol2, ind_photon1, ind_photon2;
    // TODO : Is norm necessary?
    double frac, edt, erel, tmp1, tmp2;//, norm;
    double complex tmp_coef;

    for(ist = 0; ist < pst; ist++){
        prop_mat_d[ist] = malloc(pst * sizeof(double complex));
    }

    for(ist = 0; ist < pst; ist++){
        exponent[ist] = malloc(pst * sizeof(double complex));
        propagator[ist] = malloc(pst * sizeof(double complex));
    }

    for(ist = 0; ist < pst; ist++){
        // diagonal elements
        identity[pst * ist + ist].real = 1.0;
        identity[pst * ist + ist].imag = 0.0;
        product_old[pst * ist + ist].real = 1.0;
        product_old[pst * ist + ist].imag = 0.0;
        for(jst = ist + 1; jst < pst; jst++){
            // off-diagonal elements
            // upper triangle
            identity[pst * ist + jst].real = 0.0;
            identity[pst * ist + jst].imag = 0.0;
            product_old[pst * ist + jst].real = 0.0;
            product_old[pst * ist + jst].imag = 0.0;
            // lower triangle
            identity[pst * jst + ist].real = 0.0;
            identity[pst * jst + ist].imag = 0.0;
            product_old[pst * jst + ist].real = 0.0;
            product_old[pst * jst + ist].imag = 0.0;
        }
    }

    for(ist = 0; ist < pst * pst; ist++){
        exp_idiag[ist].real = 0.0;
        exp_idiag[ist].imag = 0.0;
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
                prop_mat_d[ist][jst] = - 1.0 * tmp1 * I - tmp2;

            }
        }

        // Construct (i * total propagation matrix) to make hermitian matrix
        // exponent = (Hamiltonian - i * (NACME + decoherence)) * dt
        for(ist = 0; ist < pst; ist++){
            for (jst = 0; jst < pst; jst++){
                exponent[ist][jst] = prop_mat_d[ist][jst] * I * edt;
            }
        }

        // Convert the data type for exponent to dcomplex to exploit external math libraries
        // eigenvectors: column-major matrix
        kst = 0;
        for(ist = 0; ist < pst; ist++){
            for (jst = 0; jst < pst; jst++){
                eigenvectors[kst].real = creal(exponent[jst][ist]);
                eigenvectors[kst].imag = cimag(exponent[jst][ist]);
                kst += 1;
            }
        }

        // Diagonalize the matrix (exponent) to obtain eigenvectors and eigenvalues.
        // After this operation, eigenvectors becomes the eigenvectors defined as P.
        lwork = - 1;
        zheev_("Vectors", "Lower", &pst, eigenvectors, &pst, eigenvalues, &wkopt, &lwork, rwork, &info);
        lwork = (int)wkopt.real;
        work = (dcomplex*)malloc(lwork * sizeof(dcomplex));
        zheev_("Vectors", "Lower", &pst, eigenvectors, &pst, eigenvalues, work, &lwork, rwork, &info);
        free(work);

        // Create the diagonal matrix (exp(- i * eigenvalues))
        kst = 0;
        for(ist = 0; ist < pst; ist++){
            for (jst = 0; jst < pst; jst++){
                if(ist == jst){
                    exp_idiag[kst].real = creal(cexp(- 1.0 * eigenvalues[ist] * I));
                    exp_idiag[kst].imag = cimag(cexp(- 1.0 * eigenvalues[ist] * I));
                }
                kst += 1;
            }
        }

        // Compute the product using eigenvectors and eigenvalues
        // exp(- i * exponent) = P * exp(- i * D) * P^-1
        zgemm_("N", "N", &pst, &pst, &pst, &dcone, eigenvectors, &pst, exp_idiag, &pst, &dczero, tmp_mat, &pst);
        zgemm_("N", "C", &pst, &pst, &pst, &dcone, tmp_mat, &pst, eigenvectors, &pst, &dczero, exp_iexponent, &pst);
        // Update the product
        zgemm_("N", "N", &pst, &pst, &pst, &dcone, exp_iexponent, &pst, product_old, &pst, &dczero, product_new, &pst);
        // Backup the product
        zgemm_("N", "N", &pst, &pst, &pst, &dcone, identity, &pst, product_new, &pst, &dczero, product_old, &pst);

    }

    // Convert the data type for the product to double complex to make total propagation matrix
    kst = 0;
    for(ist = 0; ist < pst; ist++){
        for (jst = 0; jst < pst; jst++){
            // product_new: column-major matrix
            propagator[jst][ist] = product_new[kst].real + product_new[kst].imag * I;
            kst += 1;
        }
    }

    // Update the coefficients using total propagation matrix
//    zgemv_("N", &nst, &nst, &dcone, product_old, &nst, coef, 1, &dczero, tmp_coef, 1)
    for(ist = 0; ist < pst; ist++){
        tmp_coef = 0.0 + 0.0 * I;
        for (jst = 0; jst < pst; jst++){
            tmp_coef += propagator[ist][jst] * coef_d[jst];
        }
        coef_new[ist] = tmp_coef;
    }

    for(ist = 0; ist < pst; ist++){
        coef_d[ist] = coef_new[ist];
    }

    for(ist = 0; ist < pst; ist++){
        free(prop_mat_d[ist]);
    }
    for(ist = 0; ist < pst; ist++){
        free(exponent[ist]);
        free(propagator[ist]);
    }

    free(coef_new);
    free(prop_mat_d);

    free(exponent);
    free(propagator);

    free(eigenvectors);
    free(eigenvalues);

    free(exp_idiag);
    free(exp_iexponent);
    free(product_old);
    free(product_new);

    free(tmp_mat);
    free(identity);
    free(rwork);

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


