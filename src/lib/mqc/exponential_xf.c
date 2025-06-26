#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

// Complex datatype
struct _dcomplex {double real, imag;};
typedef struct _dcomplex dcomplex;

// Importing heev and gemm
extern void zheev_(char *jobz, char *uplo, int *n, dcomplex *a, int *lda, double *w, dcomplex *work, int *lwork, double *rwork, int *info);
extern void zgemm_(char *transa, char *transb, int *m, int *n, int *k, dcomplex *alpha, dcomplex *a, int *lda,
    dcomplex *b, int *ldb, dcomplex *beta, dcomplex *c, int *ldc);

// Routine for coefficient propagation scheme in exponential propagator
static void exponential_coef(int nat, int ndim, int nst, int nesteps, double dt, int *l_coh,
    double *mass, double *energy, double *energy_old, double **sigma, double **nacme,
    double **nacme_old, double **pos, double **qmom, double ***aux_pos, double ***phase, double complex *coef,
    int verbosity, double *dotpopdec);

// Interface routine for propagation scheme in exponential propagator
static void exponential(int nat, int ndim, int nst, int nesteps, double dt, char *elec_object, int *l_coh,
    double *mass, double *energy, double *energy_old, double **sigma, double **nacme, double **nacme_old,
    double **pos, double **qmom, double ***aux_pos, double ***phase, double complex *coef, int verbosity,
    double *dotpopdec){

    if(strcmp(elec_object, "coefficient") == 0){
        exponential_coef(nat, ndim, nst, nesteps, dt, l_coh, mass, energy, energy_old, sigma,
            nacme, nacme_old, pos, qmom, aux_pos, phase, coef, verbosity, dotpopdec);
    }
//    else if(strcmp(elec_object, "density") == 0){
//        exponential_rho(nat, ndim, nst, nesteps, dt, l_coh, mass, energy, energy_old, sigma,
//            nacme, nacme_old, pos, qmom, aux_pos, phase, rho, verbosity, dotpopdec);
//    }

}

static void exponential_coef(int nat, int ndim, int nst, int nesteps, double dt, int *l_coh,
    double *mass, double *energy, double *energy_old, double **sigma, double **nacme,
    double **nacme_old, double **pos, double **qmom, double ***aux_pos, double ***phase, double complex *coef,
    int verbosity, double *dotpopdec){

    double *eenergy = malloc(nst * sizeof(double));
    double **dv = malloc(nst * sizeof(double*));
    double **dec = malloc(nst * sizeof(double*));
    double complex **rho = malloc(nst * sizeof(double complex*));
    double complex **h_dec = malloc(nst * sizeof(double complex*));
    double complex *coef_new = malloc(nst * sizeof(double complex));

    // (energy - i * (NACME + decoherence)) * dt
    double complex **exponent = malloc(nst * sizeof(double complex*)); 
    // eigenvectors of (energy - i * (NACME + decoherence)) * dt, P
    dcomplex *eigenvectors = malloc((nst * nst) * sizeof(dcomplex));
    // eigenvalues of (energy - i * (NACME + decoherence)) * dt, D
    double *eigenvalues = malloc(nst * sizeof(double));

    // diagonal matrix using eigenvalues, exp(- i * D)
    dcomplex *exp_idiag = malloc((nst * nst) * sizeof(dcomplex));
    // exp(- i * exponent) = P * exp(- i * D) * P^-1
    dcomplex *exp_iexponent = malloc((nst * nst) * sizeof(dcomplex));
    // product of (P * exp(- i * D) * P^-1) until previous step
    dcomplex *product_old = malloc((nst * nst) * sizeof(dcomplex));
    // product of (P * exp(- i * D) * P^-1) until current step
    dcomplex *product_new = malloc((nst * nst) * sizeof(dcomplex));
    // final product of exp(- i * exponent)
    double complex **propagator = malloc(nst * sizeof(double complex*));

    dcomplex *tmp_mat = malloc((nst * nst) * sizeof(dcomplex));
    dcomplex *identity = malloc((nst * nst) * sizeof(dcomplex));

    dcomplex wkopt;
    dcomplex *work;
    double *rwork = malloc((3 * nst - 2) * sizeof(double));
    int lwork, info;

    int ist, jst, isp, iat, iestep;
    double frac, edt;
    double complex tmp_coef;

    for(ist = 0; ist < nst; ist++){
        dv[ist] = malloc(nst * sizeof(double));
        dec[ist] = malloc(nst * sizeof(double));
        rho[ist] = malloc(nst * sizeof(double complex));
        h_dec[ist] = malloc(nst * sizeof(double complex));
        exponent[ist] = malloc(nst * sizeof(double complex));
        propagator[ist] = malloc(nst * sizeof(double complex));
    }

    dcomplex dcone = {1.0, 0.0};
    dcomplex dczero = {0.0, 0.0};

    for(ist = 0; ist < nst; ist++){
        // Diagonal element
        identity[nst * ist + ist].real = 1.0;
        identity[nst * ist + ist].imag = 0.0;
        product_old[nst * ist + ist].real = 1.0;
        product_old[nst * ist + ist].imag = 0.0;
        for(jst = ist + 1; jst < nst; jst++){
            // Off-diagonal elements
            // Upper triangle
            identity[nst * ist + jst].real = 0.0;
            identity[nst * ist + jst].imag = 0.0;
            product_old[nst * ist + jst].real = 0.0;
            product_old[nst * ist + jst].imag = 0.0;
            // Lower triangle
            identity[nst * jst + ist].real = 0.0;
            identity[nst * jst + ist].imag = 0.0;
            product_old[nst * jst + ist].real = 0.0;
            product_old[nst * jst + ist].imag = 0.0;
        }
    }

    for(ist = 0; ist < nst * nst; ist++){
        exp_idiag[ist].real = 0.0;
        exp_idiag[ist].imag = 0.0;
    }

    /*
    // TODO : Use memset
    memset(identity, 0, (nst * nst)*sizeof(identity[0]));
    memset(product_old, 0, (nst * nst)*sizeof(product_old[0]));
    memset(exp_idiag, 0, (nst * nst)*sizeof(exp_idiag[0]));

    for(ist = 0; ist < nst; ist++){
        identity[nst * ist + ist].real = 1.0;
        product_old[nst * ist + ist].real = 1.0;
    }
    */

    frac = 1.0 / (double)nesteps;
    edt = dt * frac;

    for(iestep = 0; iestep < nesteps; iestep++){

        // Initialize variables related to decoherence
        for(iat = 0; iat < nat; iat++){
            for(isp = 0; isp < ndim; isp++){
                qmom[iat][isp] = 0.0;
            }
        }

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                dec[ist][jst] = 0.0;
                h_dec[ist][jst] = 0.0 + 0.0 * I;
            }
        }

        // Calculate densities from current coefficients
        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                rho[jst][ist] = (conj(coef[jst]) * coef[ist]);
            }
        }

        // Get quantum momentum from auxiliary positions and sigma values
        for(ist = 0; ist < nst; ist++){ 

            if(l_coh[ist] == 1){
                for(iat = 0; iat < nat; iat++){
                    for(isp = 0; isp < ndim; isp++){
                        qmom[iat][isp] += 0.5 * rho[ist][ist] * (pos[iat][isp] - aux_pos[ist][iat][isp])
                            / pow(sigma[iat][isp], 2.0) / mass[iat];
                    }
                }
            }
        }

        // Get decoherence term from quantum momentum and phase
        for(ist = 0; ist < nst; ist++){
            for(jst = ist + 1; jst < nst; jst++){ 

                if(l_coh[ist] == 1 && l_coh[jst] == 1){
                    for(iat = 0; iat < nat; iat++){
                        for(isp = 0; isp < ndim; isp++){
                            dec[ist][jst] += qmom[iat][isp] * (phase[ist][iat][isp] - phase[jst][iat][isp]);
                        }
                    }
                }
                dec[jst][ist] = - 1.0 * dec[ist][jst];

            }
        }

        // Get hamiltonian contribution from decoherence term
        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                h_dec[ist][jst] -= rho[jst][ist] * dec[jst][ist] * I;
            }
        }

        // Interpolate energy and NACME terms between time t and t + dt
        for(ist = 0; ist < nst; ist++){
            eenergy[ist] = energy_old[ist] + (energy[ist] - energy_old[ist]) * (double)iestep * frac;
            for(jst = 0; jst < nst; jst++){
                dv[ist][jst] = nacme_old[ist][jst] + (nacme[ist][jst] - nacme_old[ist][jst])
                    * (double)iestep * frac;
            }
        }

        // Construct (i * propagation matrix) to make hermitian matrix
        for(ist = 0; ist < nst; ist++){
            for (jst = 0; jst < nst; jst++){
                if (ist == jst){
                    exponent[ist][jst] = (eenergy[ist] - eenergy[0]) * edt;
                }
                else{
                    exponent[ist][jst] = (- 1.0 * I * dv[ist][jst] + h_dec[ist][jst]) * edt;
                }
            }
        }

        // Convert the data type for exponent to dcomplex to exploit external math libraries
        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                eigenvectors[nst * ist + jst].real = creal(exponent[jst][ist]);
                eigenvectors[nst * ist + jst].imag = cimag(exponent[jst][ist]);
            }
        }

        // Diagonalize the matrix (exponent) to obtain eigenvectors and eigenvalues
        // After this operation, eigenvectors becomes the eigenvectors defined as P
        lwork = -1;
        zheev_("Vectors", "Lower", &nst, eigenvectors, &nst, eigenvalues, &wkopt, &lwork, rwork, &info);
        lwork = (int)wkopt.real;
        work = (dcomplex*)malloc(lwork * sizeof(dcomplex));
        zheev_("Vectors", "Lower", &nst, eigenvectors, &nst, eigenvalues, work, &lwork, rwork, &info);
        free(work);

        // Create the diagonal matrix (exp_idiag = exp(- i * D))
        for(ist = 0; ist < nst; ist++){
            exp_idiag[nst * ist + ist].real = creal(cexp(- 1.0 * eigenvalues[ist] * I));
            exp_idiag[nst * ist + ist].imag = cimag(cexp(- 1.0 * eigenvalues[ist] * I));
        }

        // Compute the product (P * exp(- i * D) * P^-1) and update the product for every electronic step
        zgemm_("N", "N", &nst, &nst, &nst, &dcone, eigenvectors, &nst, exp_idiag, &nst, &dczero, tmp_mat, &nst);
        zgemm_("N", "C", &nst, &nst, &nst, &dcone, tmp_mat, &nst, eigenvectors, &nst, &dczero, exp_iexponent, &nst);
        // Update the product
        zgemm_("N", "N", &nst, &nst, &nst, &dcone, exp_iexponent, &nst, product_old, &nst, &dczero, product_new, &nst);
        // Backup the product
        zgemm_("N", "N", &nst, &nst, &nst, &dcone, identity, &nst, product_new, &nst, &dczero, product_old, &nst);
    }

    // Convert the data type for the term (exp(- i * exponent)) to double complex to make propagation matrix
    for(ist = 0; ist < nst; ist++){
        for(jst = 0; jst < nst; jst++){
            propagator[ist][jst] = product_new[nst * jst + ist].real + product_new[nst * jst + ist].imag * I;
        }
    }

    // Update the coefficients using the propagation matrix
    // TODO Is it necessary to change this to zgemv?
//    zgemv_("N", &nst, &nst, &dcone, product_old, &nst, coef, 1, &dczero, tmp_coef, 1)
    for(ist = 0; ist < nst; ist++){
        tmp_coef = 0.0 + 0.0 * I;
        for(jst = 0; jst < nst; jst++){
            tmp_coef += propagator[ist][jst] * coef[jst];
        }
        coef_new[ist] = tmp_coef;
    }

    for(ist = 0; ist < nst; ist++){
        coef[ist] = coef_new[ist];
    }

    for(ist = 0; ist < nst; ist++){
        free(dv[ist]);
        free(dec[ist]);
        free(rho[ist]);
        free(h_dec[ist]);
        free(exponent[ist]);
        free(propagator[ist]);
    }

    free(eenergy);
    free(dv);
    free(dec);
    free(rho);
    free(h_dec);
    free(coef_new);

    free(exponent);
    free(eigenvectors);
    free(eigenvalues);

    free(exp_idiag);
    free(exp_iexponent);
    free(product_old);
    free(product_new);
    free(propagator);

    free(tmp_mat);
    free(identity);
    free(rwork);

}
