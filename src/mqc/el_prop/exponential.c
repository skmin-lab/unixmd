#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

// Complex datatype 
struct _dcomplex { double real, imag; };
typedef struct _dcomplex dcomplex;

// importing heev and gemm
extern void zheev_(char* jobz, char* uplo, int* n, dcomplex* a, int* lda, double* eigenvalue, dcomplex* work, int* lwork, double* rwork, int* info);
extern void zgemm_(char* transpa, char* transpb, int* m, int* n, int* k, dcomplex* alpha, dcomplex* a, int* lda, 
    dcomplex* b, int* ldb, dcomplex* beta, dcomplex* c, int* ldc);

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

    // ((energy-i*nacme)*dt) = PDP^-1 (diagonalization), e^((-i*energy-nacme)*dt) = P(exp(-iD))P^-1,  
    // e^((-i*energy-nacme)*dt_1)*e^((-i*energy-nacme)*dt_2)~~~ = P_1(exp(-iD_1))P_1^-1*P_2(exp(-iD_2))P_2^-1~~~
    double *eenergy = malloc(nst * sizeof(double));
    double *eigenvalue = malloc(nst * sizeof(double)); // eigenvalue
    double *rwork = malloc((3 * nst - 2) * sizeof(double));  // temperary value for zheev
    double complex *emt = malloc((nst * nst) * sizeof(double complex)); // energy - tou(nacme)
    double complex *coef_new = malloc(nst * sizeof(double complex));  // need to calculate coef
    double complex *exp_iemt = malloc((nst *nst) * sizeof(double complex)); // double complex type of exp(i*emt)
    dcomplex *emt_dcom = malloc((nst * nst) * sizeof(dcomplex));  // dcomplex type of emt
    dcomplex *p_dcom = malloc((nst * nst) * sizeof(dcomplex));  // p_dcom is eigenvector
    dcomplex *pd_dcom= malloc((nst * nst) * sizeof(dcomplex));      // pd_dcom is PDP^-1 > (PD) part
    dcomplex *pdp_dcom = malloc((nst *nst) * sizeof(dcomplex)); // pdp_dcom is pdp 
    dcomplex *diag_dcom = malloc((nst * nst) * sizeof(dcomplex));     // diagonal matrix using eigenvalue 
    dcomplex *product_pdps_dcom = malloc((nst * nst) * sizeof(dcomplex)); // product_pdps_dcom is product of PDPs (PDP^-1)
    dcomplex *identity_temp_dcom = malloc((nst *nst) * sizeof(dcomplex)); // temperary value for zgemm (identity matrix)
    dcomplex *product_pdps_temp_dcom = malloc((nst *nst) * sizeof(dcomplex)); // temperary value for zgemm (product of pdp)
    double **dv = malloc(nst * sizeof(double*));

    int ist, jst, iestep, lwork, info;  // lwork : The length of the array WORK, info : confirmation that heev is working
    double frac, edt; 
    double complex tem_coef = 0.0;

    for(ist = 0; ist < nst; ist++){
        dv[ist] = malloc(nst * sizeof(double));
    }

    dcomplex dcone = {1.0, 0.0}; 
    dcomplex dczero = {0.0, 0.0};
    dcomplex wkopt;  // need to get optimized lwork
    dcomplex* work;  // length of lwork

    // Set zero matrix
    for(ist = 0; ist < nst; ist++){
        for(jst = 0; jst < nst; jst++){
            diag_dcom[ist * nst + jst].real = 0.0;
            diag_dcom[ist * nst + jst].imag = 0.0;
            product_pdps_dcom[ist * nst + jst].real = 0.0;
            product_pdps_dcom[ist * nst + jst].imag = 0.0;
        }
    }

    // TODO : use memset
    // memset(diag_dcom, 0, (nst*nst)*sizeof(diag_dcom[0]));
    // memset(product_pdps_dcom, 0, (nst*nst)*sizeof(product_pdps_dcom[0]));
    
    // Set identity matrix
    for(ist = 0; ist < nst; ist++){
        product_pdps_dcom[nst * ist + ist].real = 1.0;
    }
   
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

        // Make emt, emt means Energy - Tou(nacme)
        for(ist = 0; ist < nst; ist++){
            for (jst = 0; jst < nst; jst++){
                if (ist == jst){
                    emt[nst * ist + jst] = (eenergy[ist] - eenergy[0]) * edt;
                }
                else{
                    emt[nst * ist + jst] = - 1.0 * I * dv[jst][ist] * edt;
                }         
            }
        }

        // Change complex type
        for(ist = 0; ist < nst * nst; ist++){
            emt_dcom[ist].real = creal(emt[ist]);
            emt_dcom[ist].imag = cimag(emt[ist]);
            p_dcom[ist].real = creal(emt[ist]);
            p_dcom[ist].imag = cimag(emt[ist]);
        }    
        
        // diagonaliztion 
        lwork = -1;
        zheev_("Vectors", "Lower", &nst, p_dcom, &nst, eigenvalue, &wkopt, &lwork, rwork, &info);
        lwork = (int)wkopt.real;
        work = (dcomplex*)malloc(lwork * sizeof(dcomplex));
        zheev_("Vectors", "Lower", &nst, p_dcom, &nst, eigenvalue, work, &lwork, rwork, &info); // p_dcom -> P(unitary matrix which is consisting of eigen vector)

        // Make diagonal matrix D
        for(ist = 0; ist < nst; ist++){ 
                diag_dcom[nst * ist + ist].real = creal(cexp(- 1.0 * eigenvalue[ist] * I));
                diag_dcom[nst * ist + ist].imag = cimag(cexp(- 1.0 * eigenvalue[ist] * I));
        }

        // Matrix multiplication PDP^-1 and 
        zgemm_("N","N",&nst, &nst, &nst, &dcone, p_dcom, &nst, diag_dcom, &nst, &dczero, pd_dcom, &nst); // P*D
        zgemm_("N","C",&nst, &nst, &nst, &dcone, pd_dcom, &nst, p_dcom, &nst, &dczero, pdp_dcom, &nst); // PD * P^-1
        zgemm_("N","N",&nst, &nst, &nst, &dcone, pdp_dcom, &nst, product_pdps_dcom, &nst, &dczero, product_pdps_temp_dcom, &nst); // PDP^-1 * (old PDP^-1)

        //reset the diag_dcom
        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                diag_dcom[ist * nst + jst].real = 0.0;
                diag_dcom[ist * nst + jst].imag = 0.0;
            }
        }

        //reset the identity_temp_dcom
        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                identity_temp_dcom[ist * nst + jst].real = 0.0;
                identity_temp_dcom[ist * nst + jst].imag = 0.0;
            }
        }

        // TODO : use memset
        // memset(diag_dcom, 0, (nst*nst)*sizeof(diag_dcom[0]));
        // memset(identity_temp_dcom, 0, (nst*nst)*sizeof(identity_temp_dcom[0]));

        for(ist = 0; ist < nst; ist++){
            identity_temp_dcom[nst * ist + ist].real = 1.0;
        }

        // update coefficent 
        zgemm_("N","N", &nst, &nst, &nst, &dcone, identity_temp_dcom, &nst, product_pdps_temp_dcom, &nst, &dczero, product_pdps_dcom, &nst); // to keep PDP^-1 value in total_coef_dcom 
    }

    //change complex type
    for(ist = 0; ist < nst * nst; ist++){
        exp_iemt[ist] = product_pdps_dcom[ist].real + product_pdps_dcom[ist].imag * I;
    }

    // matrix - vector multiplication //TODO Is it necessary to change this to zgemv?
    //zgemv_("N", &nst, &nst, &dcone, product_pdps_dcom, &nst, coef, 1, &dczero, tem_coef, 1)
    for(ist = 0; ist < nst; ist++){
        for (jst = 0; jst < nst; jst++){
            tem_coef += exp_iemt[nst * jst + ist] * coef[jst];
        }
        coef_new[ist] = tem_coef;
        tem_coef = 0.0;
    }

    for(ist = 0; ist < nst; ist++){
        coef[ist] = coef_new[ist];        
    }

    for(ist = 0; ist < nst; ist++){
            free(dv[ist]);
    }
  
    free(coef_new);
    free(exp_iemt);
    free(pdp_dcom);
    free(identity_temp_dcom);
    free(product_pdps_temp_dcom);
    free(product_pdps_dcom);
    free(p_dcom);
    free(emt_dcom);
    free(diag_dcom);
    free(pd_dcom);
    free(eigenvalue);
    free(rwork);
    free(emt);
    free(eenergy);
    free(dv);

}