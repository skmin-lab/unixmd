#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* Complex datatype */
struct _dcomplex { double real, imag; };
typedef struct _dcomplex dcomplex;

/* ZHEEV prototype */
extern void zheev_( char* jobz, char* uplo, int* n, dcomplex* a, int* lda, double* w, dcomplex* work, int* lwork, double* rwork, int* info );
extern void zgemm_(char* , char*, int*, int*, int*, dcomplex*, dcomplex*, int*, dcomplex*, int*, dcomplex*, dcomplex*, int*);

// Routine for coefficient propagation scheme in exponential propagator
static void expon_coef(int nst, int nesteps, double dt, double *energy, double *energy_old,
    double **nacme, double **nacme_old, double complex *coef);

// Routine for density propagation scheme in exponential propagator
/*static void expon_rho(int nst, int nesteps, double dt, double *energy, double *energy_old,
    double **nacme, double **nacme_old, double complex **rho);
*/
// Interface routine for propagation scheme in exponential propagator
static void expon(int nst, int nesteps, double dt, char *elec_object, double *energy, double *energy_old,
    double **nacme, double **nacme_old, double complex *coef, double complex **rho){

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
    double *w = malloc(nst * sizeof(double));
    double *rwork = malloc((3*nst-2) * sizeof(double));
    double complex *emt = malloc((nst*nst) * sizeof(double complex)); 
    double complex *coef_new = malloc(nst * sizeof(double complex));
    double complex *total_coef_com = malloc((nst *nst) * sizeof(double complex));
    dcomplex *matrix = malloc((nst*nst) * sizeof(dcomplex));
    dcomplex *emtt = malloc((nst*nst) * sizeof(dcomplex));
    dcomplex *Dia = malloc((nst*nst) * sizeof(dcomplex)); 
    dcomplex *total_coef= malloc((nst*nst) * sizeof(dcomplex)); 
     
    int ist, jst, iestep,lwork, info;
    double frac, edt; 

    dcomplex Alpha = {1.0,0.0};
    dcomplex Beta = {0.0,0.0};
    dcomplex wkopt;
    dcomplex* work;

    // Set zero matrix
    memset(Dia, 0, (nst*nst)*sizeof(Dia[0]));
    memset(total_coef, 0, (nst*nst)*sizeof(total_coef[0])); 
    
    // Set identity matrix
    for(ist =0; ist<nst; ist++){
        total_coef[nst*ist+ist].real = 1;
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

        // emt means Energy - Tou(nacme)
        for(ist = 0; ist < nst; ist++){
            for (jst = 0; jst < nst; jst++){
                if (ist == jst){
                    emt[nst*ist+jst] = (eenergy[ist]- eenergy[0])*edt;
                }else{
                    emt[nst*ist+jst] = (-I*dv[jst][ist])*edt;
                }         
            }
        }

        // change complex type
        for(ist=0; ist<nst*nst; ist++){
            matrix[ist].real = creal(emt[ist]);
            matrix[ist].imag = cimag(emt[ist]);
        }    
        
        // Diagonaliztion 
        lwork = -1;
        zheev_( "Vectors", "Lower", &nst, matrix, &nst, w, &wkopt, &lwork, rwork, &info );
        lwork = (int)wkopt.real;
        work = (dcomplex*)malloc( lwork*sizeof(dcomplex) );
        zheev_( "Vectors", "Lower", &nst, matrix, &nst, w, work, &lwork, rwork, &info );

        // Make diagonal matrix D
        for(ist = 0; ist<nst; ist++){ 
                Dia[nst*ist+ist].real = creal(cexp(-w[ist]*I));
                Dia[nst*ist+ist].imag = cimag(cexp(-w[ist]*I));
        }

        // Matrix multiplication PDP^-1 and 
        zgemm_("N","N",&nst,&nst,&nst,&Alpha,matrix, &nst,Dia, &nst, &Beta, emtt, &nst);
        zgemm_("N","C",&nst,&nst,&nst,&Alpha,emtt, &nst,matrix, &nst, &Beta, Dia, &nst);
        zgemm_("N","N",&nst,&nst,&nst,&Alpha,Dia, &nst,total_coef, &nst, &Beta, matrix, &nst);

        // reused matrix by removing date 
        memset(Dia, 0, (nst*nst)*sizeof(Dia[0]));

        for(ist =0; ist<nst; ist++){
                Dia[nst*ist+ist].real = 1;
        }

        // update coefficent 
        zgemm_("N","N", &nst, &nst, &nst, &Alpha, Dia,&nst, matrix, &nst,&Beta,total_coef,&nst);
    }

    //zgemv_("N", &nst, &nst, &Alpha, total_coef, &nst, coef, 1, &Beta, tem_coef, 1)

    //change complex type
    for(ist=0; ist<nst*nst; ist++){
        total_coef_com[ist] = total_coef[ist].real + total_coef[ist].imag * I;
    }
    // matrix - vector multiplication //TODO Is it necessary to change this to zgemv?
    double complex tem_coef = 0;
    for(ist = 0; ist < nst; ist++){
        for (jst = 0; jst <nst; jst++){
            tem_coef += total_coef_com[nst*(jst)+ist]*coef[jst];
    
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
    free(total_coef_com);
    free(total_coef);
    free(matrix);
    free(Dia);
    free(emtt);
    free(w);
    free(rwork);
    free(emt);
    free(eenergy);
    free(dv);

}