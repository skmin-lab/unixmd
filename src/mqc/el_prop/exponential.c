#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

// Complex datatype 
struct _dcomplex { double real, imag; };
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
/*    else if(strcmp(elec_object, "density") == 0){
        expon_rho(nst, nesteps, dt, energy, energy_old, nacme, nacme_old, rho);
    }*/

}
  
static void expon_coef(int nst, int nesteps, double dt, double *energy, double *energy_old, 
    double **nacme, double **nacme_old, double complex *coef){

    // (( energy - i * NACME ) * dt) = PDP^-1 (diagonalization), exp(( - i * energy - NACME ) * dt ) = P( exp( -iD ) )P^-1,  
    // exp( ( - i * energy_1 - NACME_1 ) * dt ) * exp(( - i * energy_2 - NACME_2 ) * dt )···exp(( - i * energy_n - NACME_n ) * dt ) 
    // = P_1( exp( - i * D_1 ) )P_1^-1 * P_2( exp( - i * D_2 ) )P_2^-1 ··· P_n( exp( - i * D_n ) )P_n^-1
    // n is interpolation step number
    
    // \[(energy - i\cdot NACME)\cdot dt = PDP^{-1} ,\](diagonalization)  \[exp(-i\cdot (energy-i\cdot NACME)\cdot dt)= P exp(-i\cdot D) P^{-1}\]
    // \[exp( -i\cdot (energy_1 -i\cdot NACME_1 ) \cdot dt ) \cdot exp( - i \cdot (energy_2 -i\cdot NACME_2)  \cdot dt) \cdots exp( -i\cdot (energy_n -i\cdot NACME_n ) \cdot dt )\] 
    // \[= P_1( exp( - i \cdot D_1 ) )P_1^{-1} \cdot P_2( exp( - i \cdot D_2 ) )P_2^{-1} \cdots P_n( exp( - i \cdot D_n ) )P_n^{-1} \]
    // energy and NACME are interpolated. Because of that, diagonal matrix and eigenvector is also changed. They are different value at different time step.
    double *eenergy = malloc(nst * sizeof(double));
    double *eigenvalue = malloc(nst * sizeof(double)); // eigenvalue of (energy - i * NACME) * dt
    double *rwork = malloc((3 * nst - 2) * sizeof(double));  // temporary value for zheev
    double complex *exponent = malloc((nst * nst) * sizeof(double complex)); // (energy - i * NACME) * dt  
    double complex *coef_new = malloc(nst * sizeof(double complex));  // need to calculate coef
    double complex *exp_iexponent = malloc((nst * nst) * sizeof(double complex)); // double complex type of exp(-i * exponent * dt)
    dcomplex *diag_dcom = malloc((nst * nst) * sizeof(dcomplex));     // diagonal matrix using eigenvalue, exp(-iD), D is diagonal matrix and diagonal elements are eigenvalue of (energy - i * NACME) * dt
    dcomplex *eigenvectors = malloc((nst * nst) * sizeof(dcomplex));  // it is eigenvectors
    dcomplex *tmp_mat_dcom= malloc((nst * nst) * sizeof(dcomplex));
    dcomplex *dg_exp_iexponent = malloc((nst * nst) * sizeof(dcomplex)); // dg_exp_iexponent is diagonalized exp(-i*exponent), Pexp(-iD)P^-1 
    dcomplex *product_old = malloc((nst * nst) * sizeof(dcomplex)); // product_old is product of Pexp(-iD)P^-1 until previous step (old)
    dcomplex *identity_dcom = malloc((nst * nst) * sizeof(dcomplex)); // identity matrix
    dcomplex *product_new = malloc((nst * nst) * sizeof(dcomplex)); // product_new is product of Pexp(-iD)P^-1 until current step (new)
    double **dv = malloc(nst * sizeof(double*));

    int ist, jst, iestep, lwork, info;  // lwork : The length of the array WORK, info : confirmation that heev is working
    double frac, edt; 
    double complex tmp_coef;

    for(ist = 0; ist < nst; ist++){
        dv[ist] = malloc(nst * sizeof(double));
    }

    dcomplex dcone = {1.0, 0.0}; 
    dcomplex dczero = {0.0, 0.0};
    dcomplex wkopt;  // need to get optimized lwork
    dcomplex *work;  // length of lwork

    // Set zero matrix
    for(ist = 0; ist < nst; ist++){
        for(jst = 0; jst < nst; jst++){
            identity_dcom[ist * nst + jst].real = 0.0;
            identity_dcom[ist * nst + jst].imag = 0.0;
            product_old[ist * nst + jst].real = 0.0;
            product_old[ist * nst + jst].imag = 0.0;
        }
    }

    // TODO : Use memset
    // memset(identity_dcom, 0, (nst*nst)*sizeof(identity_dcom[0]));
    // memset(product_old, 0, (nst*nst)*sizeof(product_old[0]));
    
    // Set identity matrix
    for(ist = 0; ist < nst; ist++){
        identity_dcom[nst * ist + ist].real = 1.0;
        product_old[nst * ist + ist].real = 1.0;
    }
   
    frac = 1.0 / (double)nesteps;
    edt = dt * frac;

    for(iestep = 0; iestep < nesteps; iestep++){

        // TODO : Use memset
        // memset(diag_dcom, 0, (nst*nst)*sizeof(diag_dcom[0]));

        // Reset the diagonal matrix
        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                diag_dcom[ist * nst + jst].real = 0.0;
                diag_dcom[ist * nst + jst].imag = 0.0;
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
        // exponent = energy - i * NACME 
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
        
        // Diagonalize the matrix (exponent) to obtain eigenvectors and eigenvalues.
        // After this operation, eigenvectors becomes the eigenvectors defined as P.
        lwork = -1;
        zheev_("Vectors", "Lower", &nst, eigenvectors, &nst, eigenvalue, &wkopt, &lwork, rwork, &info);
        lwork = (int)wkopt.real;
        work = (dcomplex*)malloc(lwork * sizeof(dcomplex));
        zheev_("Vectors", "Lower", &nst, eigenvectors, &nst, eigenvalue, work, &lwork, rwork, &info); // eigenvectors -> P(unitary matrix is consisting of eigenvector)
        free(work);

        // Create the diagonal matrix (diag_dcom = exp( -i*D )) where D is a matrix consisting of eigenvalues obtained from upper operation.
        for(ist = 0; ist < nst; ist++){ 
                diag_dcom[nst * ist + ist].real = creal(cexp(- 1.0 * eigenvalue[ist] * I));
                diag_dcom[nst * ist + ist].imag = cimag(cexp(- 1.0 * eigenvalue[ist] * I));
        }

        // Compute the product ( P*exp( -i*D )*P^-1  ) and update the product for every electronic step
        zgemm_("N", "N", &nst, &nst, &nst, &dcone, eigenvectors, &nst, diag_dcom, &nst, &dczero, tmp_mat_dcom, &nst); // P*exp(-iD)
        zgemm_("N", "C", &nst, &nst, &nst, &dcone, tmp_mat_dcom, &nst, eigenvectors, &nst, &dczero, dg_exp_iexponent, &nst); // Pexp(-iD) * P^-1 
        zgemm_("N", "N", &nst, &nst, &nst, &dcone, dg_exp_iexponent, &nst, product_old, &nst, &dczero, product_new, &nst); // Pexp(-iD)P^-1  * (old Pexp(-iD)P^-1 )

        // Update coefficent 
        zgemm_("N", "N", &nst, &nst, &nst, &dcone, identity_dcom, &nst, product_new, &nst, &dczero, product_old, &nst); // to keep Pexp(-iD)P^-1 value in total_coef_dcom 
    }

    // Convert the data type for the term ( exp( - i * exponent ) ) to original double complex to make propagation matrix
    for(ist = 0; ist < nst * nst; ist++){
        exp_iexponent[ist] = product_old[ist].real + product_old[ist].imag * I;
    }

    // matrix - vector multiplication //TODO Is it necessary to change this to zgemv?
    //zgemv_("N", &nst, &nst, &dcone, product_old, &nst, coef, 1, &dczero, tmp_coef, 1)
    for(ist = 0; ist < nst; ist++){
        tmp_coef = 0.0 + 0.0 * I;
        for (jst = 0; jst < nst; jst++){
            tmp_coef += exp_iexponent[nst * jst + ist] * coef[jst];
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
    free(exp_iexponent);
    free(dg_exp_iexponent);
    free(identity_dcom);
    free(product_new);
    free(product_old);
    free(eigenvectors);
    free(diag_dcom);
    free(tmp_mat_dcom);
    free(eigenvalue);
    free(rwork);
    free(exponent);
    free(eenergy);
    free(dv);

}
