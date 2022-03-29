#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

// Complex datatype 
struct _dcomplex { double real, imag; };
typedef struct _dcomplex dcomplex;

// importing heev and gemm
extern void zheev_(char* jobz, char* uplo, int* n, dcomplex* a, int* lda, double* w, dcomplex* work, int* lwork, double* rwork, int* info );
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

    double *eenergy = malloc(nst * sizeof(double));
    double *w = malloc(nst * sizeof(double)); // eigenvalue
    double *rwork = malloc((3 * nst - 2) * sizeof(double));  // DOUBLE PRECISION array
    double complex *emt = malloc((nst * nst) * sizeof(double complex)); // energy - tou(nacme)
    double complex *coef_new = malloc(nst * sizeof(double complex));  // need to calculate coef
    double complex *total_coef = malloc((nst *nst) * sizeof(double complex)); // double complex type of total coefficient
    dcomplex *emt_dcom = malloc((nst * nst) * sizeof(dcomplex));  // dcomplex type of emt
    dcomplex *pd = malloc((nst * nst) * sizeof(dcomplex));      // pd is PDP^-1 > (PD) part
    dcomplex *dia = malloc((nst * nst) * sizeof(dcomplex));     // diagonal matrix using eigenvalue 
    dcomplex *total_coef_dcom = malloc((nst * nst) * sizeof(dcomplex)); // total coefficient (PDP^-1)
    dcomplex *temporary_value = malloc((nst * nst) * sizeof(dcomplex)); // need to move coefficient
     
    int ist, jst, iestep, lwork, info;  // lwork : The length of the array WORK, info : confirmation that heev is working
    double frac, edt; 

    dcomplex C1 = {1.0, 0.0}; 
    dcomplex C0 = {0.0, 0.0};
    dcomplex wkopt;
    dcomplex* work;

    // Set zero matrix
    for(ist = 0; ist < nst; ist++){
        for(jst = 0; jst < nst; jst++){
            dia[ist * nst + jst].real = 0;
            dia[ist * nst + jst].imag = 0;
            total_coef_dcom[ist * nst + jst].real = 0;
            total_coef_dcom[ist * nst + jst].imag = 0;
        }
    }

    // TODO : use memset
    // memset(dia, 0, (nst*nst)*sizeof(dia[0]));
    // memset(total_coef_dcom, 0, (nst*nst)*sizeof(total_coef_dcom[0]));
    // memset(temporary_value, 0, (nst*nst)*sizeof(temporay_value[0])); 
    
    // Set identity matrix
    for(ist = 0; ist < nst; ist++){
        total_coef_dcom[nst * ist + ist].real = 1;
    }

    double **dv = malloc(nst * sizeof(double*));
    for(ist = 0; ist < nst; ist++){
        dv[ist] = malloc(nst * sizeof(double));
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
                    emt[nst*ist+jst] = (eenergy[ist] - eenergy[0]) * edt;
                }
                else{
                    emt[nst*ist+jst] = -I * dv[jst][ist] * edt;
                }         
            }
        }

        // Change complex type
        for(ist = 0; ist < nst * nst; ist++){
            emt_dcom[ist].real = creal(emt[ist]);
            emt_dcom[ist].imag = cimag(emt[ist]);
        }    
        
        // diagonaliztion 
        lwork = -1;
        zheev_( "Vectors", "Lower", &nst, emt_dcom, &nst, w, &wkopt, &lwork, rwork, &info );
        lwork = (int)wkopt.real;
        work = (dcomplex*)malloc(lwork*sizeof(dcomplex));
        zheev_( "Vectors", "Lower", &nst, emt_dcom, &nst, w, work, &lwork, rwork, &info );

        // Make diagonal matrix D
        for(ist = 0; ist < nst; ist++){ 
                dia[nst*ist+ist].real = creal(cexp(-w[ist]*I));
                dia[nst*ist+ist].imag = cimag(cexp(-w[ist]*I));
        }

        // Matrix multiplication PDP^-1 and 
        zgemm_("N","N",&nst, &nst, &nst, &C1, emt_dcom, &nst, dia, &nst, &C0, pd, &nst);
        zgemm_("N","C",&nst, &nst, &nst, &C1, pd, &nst, emt_dcom, &nst, &C0, dia, &nst);
        zgemm_("N","N",&nst, &nst, &nst, &C1, dia, &nst, total_coef_dcom, &nst, &C0, emt_dcom, &nst);

        //reset the temporary value
        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                temporary_value[ist * nst + jst].real = 0;
                temporary_value[ist * nst + jst].imag = 0;
            }
        }

        for(ist = 0; ist < nst; ist++){
            temporary_value[nst * ist + ist].real = 1;
        }

        // update coefficent 
        zgemm_("N","N", &nst, &nst, &nst, &C1, temporary_value, &nst, emt_dcom, &nst, &C0, total_coef_dcom, &nst);
    }

    //change complex type
    for(ist = 0; ist < nst * nst; ist++){
        total_coef[ist] = total_coef_dcom[ist].real + total_coef_dcom[ist].imag * I;
    }

    // matrix - vector multiplication //TODO Is it necessary to change this to zgemv?
    //zgemv_("N", &nst, &nst, &C1, total_coef_dcom, &nst, coef, 1, &C0, tem_coef, 1)
    double complex tem_coef = 0;
    for(ist = 0; ist < nst; ist++){
        for (jst = 0; jst < nst; jst++){
            tem_coef += total_coef[nst*jst+ist] * coef[jst];
        }
        coef_new[ist] = tem_coef;
        tem_coef = 0;
    }

    for(ist = 0; ist < nst; ist++){
        coef[ist] = coef_new[ist];        
    }

    for(ist = 0; ist < nst; ist++){
            free(dv[ist]);
    }
  
    free(coef_new);
    free(total_coef);
    free(total_coef_dcom);
    free(emt_dcom);
    free(dia);
    free(pd);
    free(w);
    free(rwork);
    free(emt);
    free(eenergy);
    free(dv);

}