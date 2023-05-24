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
static void expon_coef(int nst, int nesteps, double dt, double *energy, double *energy_old,
    double **nacme, double **nacme_old, double complex *coef);

// Interface routine for propagation scheme in exponential propagator
static void expon(int nst, int nesteps, double dt, char *elec_object, double *energy, double *energy_old,
    double **nacme, double **nacme_old, double complex *coef){

    if(strcmp(elec_object, "coefficient") == 0){
        expon_coef(nst, nesteps, dt, energy, energy_old, nacme, nacme_old, coef);
    }

//    else if(strcmp(elec_object, "density") == 0){
//        expon_rho(nst, nesteps, dt, energy, energy_old, nacme, nacme_old, rho);
//    }

}

static void expon_coef(int nst, int nesteps, double dt, double *energy, double *energy_old,
    double **nacme, double **nacme_old, double complex *coef){

    // (E - i \sigma) dt = PDP^{-1} (diagonalization)
    // \exp( - i (E - i \sigma) dt) = P \exp( - i D) P^{-1}
    // \exp( - i (E_1 - i \sigma_1) dt ) \times \exp( - i (E_2 - i \sigma_2) dt) \times \cdots \times \exp( - i (E_n - i \sigma_n) dt) =
    // P_1(\exp( - i D_1))P_1^{-1} \times P_2(\exp( - i D_2))P_2^{-1} \times \cdots \times P_n(\exp( - i D_n))P_n^{-1}

    // where E is energy, sigma is nonadiabatic coupling matrix elements(NACME), and n is the interpolation step number,
    // P is a matrix consisting of eigenvectors of (E - i * sigma) * dt
    // D is a diagonal matrix consisting of eigenvalues of (E - i * sigma) * dt
    double *eenergy = malloc(nst * sizeof(double));
    // Eigenvalues of (energy - i * nacme) * dt
    double *eigenvalues = malloc(nst * sizeof(double));
    double *rwork = malloc((3 * nst - 2) * sizeof(double));
    // (energy - i * nacme) * dt
    double complex *exponent = malloc((nst * nst) * sizeof(double complex));
    double complex *coef_new = malloc(nst * sizeof(double complex));
    // Double complex type of exp( - i * exponent)
    double complex *propagator = malloc((nst * nst) * sizeof(double complex));
    // Diagonal matrix using eigenvalues, exp( - i * D)
    dcomplex *exp_idiag = malloc((nst * nst) * sizeof(dcomplex));
    dcomplex *eigenvectors = malloc((nst * nst) * sizeof(dcomplex));
    dcomplex *tmp_mat = malloc((nst * nst) * sizeof(dcomplex));
    dcomplex *exp_iexponent = malloc((nst * nst) * sizeof(dcomplex));
    // Product of Pexp( - i * D)P^-1 until previous step (old)
    dcomplex *product_old = malloc((nst * nst) * sizeof(dcomplex));
    // Product of Pexp( - i * D)P^-1 until current step (new)
    dcomplex *product_new = malloc((nst * nst) * sizeof(dcomplex));
    dcomplex *identity = malloc((nst * nst) * sizeof(dcomplex));
    double **dv = malloc(nst * sizeof(double*));

    int ist, jst, iestep, lwork, info;
    double frac, edt;
    double complex tmp_coef;

    for(ist = 0; ist < nst; ist++){
        dv[ist] = malloc(nst * sizeof(double));
    }

    dcomplex dcone = {1.0, 0.0};
    dcomplex dczero = {0.0, 0.0};
    dcomplex wkopt;
    dcomplex *work;

    // Set identity matrix
    for(ist = 0; ist < nst; ist++){
        // Diagonal element to one
        identity[nst * ist + ist].real = 1.0;
        identity[nst * ist + ist].imag = 0.0;
        product_old[nst * ist + ist].real = 1.0;
        product_old[nst * ist + ist].imag = 0.0;
        for(jst = ist + 1; jst < nst; jst++){
            // Off-diagonal elements to zero
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

    // Set zero matrix
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
                    exponent[nst * ist + jst] = (eenergy[ist] - eenergy[0]) * edt;
                }
                else{
                    exponent[nst * ist + jst] = - 1.0 * I * dv[jst][ist] * edt;
                }
            }
        }

        // Convert the data type for exponent to dcomplex to exploit external math libraries
        for(ist = 0; ist < nst * nst; ist++){
            eigenvectors[ist].real = creal(exponent[ist]);
            eigenvectors[ist].imag = cimag(exponent[ist]);
        }

        // Diagonalize the matrix (exponent) to obtain eigenvectors and eigenvalues
        // After this operation, eigenvectors becomes the eigenvectors defined as P
        lwork = -1;
        zheev_("Vectors", "Lower", &nst, eigenvectors, &nst, eigenvalues, &wkopt, &lwork, rwork, &info);
        lwork = (int)wkopt.real;
        work = (dcomplex*)malloc(lwork * sizeof(dcomplex));
        zheev_("Vectors", "Lower", &nst, eigenvectors, &nst, eigenvalues, work, &lwork, rwork, &info);
        free(work);

        // Create the diagonal matrix (exp_idiag = exp( - i * D)) where D is a matrix consisting of eigenvalues obtained from upper operation
        for(ist = 0; ist < nst; ist++){
                exp_idiag[nst * ist + ist].real = creal(cexp( - 1.0 * eigenvalues[ist] * I));
                exp_idiag[nst * ist + ist].imag = cimag(cexp( - 1.0 * eigenvalues[ist] * I));
        }

        // Compute the product (P * exp( - i * D) * P^-1) and update the product for every electronic step
        zgemm_("N", "N", &nst, &nst, &nst, &dcone, eigenvectors, &nst, exp_idiag, &nst, &dczero, tmp_mat, &nst);
        zgemm_("N", "C", &nst, &nst, &nst, &dcone, tmp_mat, &nst, eigenvectors, &nst, &dczero, exp_iexponent, &nst);
        zgemm_("N", "N", &nst, &nst, &nst, &dcone, exp_iexponent, &nst, product_old, &nst, &dczero, product_new, &nst);

        // Update coefficent
        zgemm_("N", "N", &nst, &nst, &nst, &dcone, identity, &nst, product_new, &nst, &dczero, product_old, &nst);
    }

    // Convert the data type for the term (exp( - i * exponent)) to original double complex to make propagation matrix
    for(ist = 0; ist < nst * nst; ist++){
        propagator[ist] = product_new[ist].real + product_new[ist].imag * I;
    }

    // Matrix - vector multiplication
    // TODO Is it necessary to change this to zgemv?
//    zgemv_("N", &nst, &nst, &dcone, product_old, &nst, coef, 1, &dczero, tmp_coef, 1)
    for(ist = 0; ist < nst; ist++){
        tmp_coef = 0.0 + 0.0 * I;
        for (jst = 0; jst < nst; jst++){
            tmp_coef += propagator[nst * jst + ist] * coef[jst];
        }
        coef_new[ist] = tmp_coef;
    }

    for(ist = 0; ist < nst; ist++){
        coef[ist] = coef_new[ist];
    }

    for(ist = 0; ist < nst; ist++){
        free(dv[ist]);
    }

    free(coef_new);
    free(propagator);
    free(exp_iexponent);
    free(identity);
    free(product_new);
    free(product_old);
    free(eigenvectors);
    free(exp_idiag);
    free(tmp_mat);
    free(eigenvalues);
    free(rwork);
    free(exponent);
    free(eenergy);
    free(dv);

}


