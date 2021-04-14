#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "derivs_ct.h"

// Routine for coefficient elec_object scheme in rk4 solver
static void rk4_coef(int nst, int nesteps, double dt, double *energy, double *energy_old,
    double **nacme, double **nacme_old, double **k_lk, double complex *coef);

// Routine for density elec_object scheme in rk4 solver
/*
static void rk4_rho(int nst, int nesteps, double dt, double *energy, double *energy_old,
    double **nacme, double **nacme_old, double **k_lk, double *dotpopd, double complex **rho);
*/

// Interface routine for elec_object scheme in rk4 solver
static void rk4(int nst, int nesteps, double dt, char *elec_object, double *energy, double *energy_old,
    double **nacme, double **nacme_old, double **k_lk, double complex *coef, double complex **rho){

    if(strcmp(elec_object, "coefficient") == 0){
        rk4_coef(nst, nesteps, dt, energy, energy_old, nacme, nacme_old, k_lk, coef);
    }
    /*
    else if(strcmp(elec_object, "density") == 0){
        rk4_rho(nst, nesteps, dt, energy, energy_old, nacme, nacme_old, k_lk, rho);
    }
    */
}

// Routine for coefficient elec_object scheme in rk4 solver
static void rk4_coef(int nst, int nesteps, double dt, double *energy, double *energy_old,
    double **nacme, double **nacme_old, double **k_lk, double complex *coef){

    double complex *k1 = malloc(nst * sizeof(double complex));
    double complex *k2 = malloc(nst * sizeof(double complex));
    double complex *k3 = malloc(nst * sizeof(double complex));
    double complex *k4 = malloc(nst * sizeof(double complex));
    double complex *kfunction = malloc(nst * sizeof(double complex));
    double complex *variation = malloc(nst * sizeof(double complex));
    double complex *c_dot = malloc(nst * sizeof(double complex));
    double complex *coef_new = malloc(nst * sizeof(double complex));
    double *eenergy = malloc(nst * sizeof(double));
//    double *na_term = malloc(nst * sizeof(double));
    double **dv = malloc(nst * sizeof(double*));

    int ist, jst, iestep;
    double frac, edt, norm;

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

        // Calculate k1
        ct_cdot(nst, eenergy, dv, k_lk, coef, c_dot);

        for(ist = 0; ist < nst; ist++){
            k1[ist] = edt * c_dot[ist];
            kfunction[ist] = 0.5 * k1[ist];
            coef_new[ist] = coef[ist] + kfunction[ist];
        }

        // Calculate k2
        ct_cdot(nst, eenergy, dv, k_lk, coef_new, c_dot);

        for(ist = 0; ist < nst; ist++){
            k2[ist] = edt * c_dot[ist];
            kfunction[ist] = 0.5 * (- 1.0 + sqrt(2.0)) * k1[ist] + (1.0 - 0.5 * sqrt(2.0)) * k2[ist];
            coef_new[ist] = coef[ist] + kfunction[ist];
        }

        // Calculate k3
        ct_cdot(nst, eenergy, dv, k_lk, coef_new, c_dot);

        for(ist = 0; ist < nst; ist++){
            k3[ist] = edt * c_dot[ist];
            kfunction[ist] = - 0.5 * sqrt(2.0) * k2[ist] + (1.0 + 0.5 * sqrt(2.0)) * k3[ist];
            coef_new[ist] = coef[ist] + kfunction[ist];
        }

        // Calculate k4
        ct_cdot(nst, eenergy, dv, k_lk, coef_new, c_dot);

        for(ist = 0; ist < nst; ist++){
            k4[ist] = edt * c_dot[ist];
            variation[ist] = (k1[ist] + (2.0 - sqrt(2.0)) * k2[ist] + (2.0 + sqrt(2.0))
                * k3[ist] + k4[ist]) / 6.0;
            coef_new[ist] = coef[ist] + variation[ist];
        }

        // TODO : Is this part necessary?
        // Renormalize the coefficients
        norm = dot(nst, coef_new, coef_new);
        for(ist = 0; ist < nst; ist++){
            coef_new[ist] /= sqrt(norm);
            coef[ist] = coef_new[ist];
        }

    }

    /* 
    for(ist = 0; ist < nst; ist++){
        na_term[ist] = 0.0;
        for(jst = 0; jst < nst; jst++){
            if(jst != ist){
                na_term[ist] -= dv[ist][jst] * coef[jst];
            }
        }
        printf("RHODOT_NAC %d %f\n",ist+1,2.0*creal(conj(coef[ist])) * na_term[ist]);
    }
    printf("RK4_COEF : NORM = %15.8f\n", creal(norm));
    */

    for(ist = 0; ist < nst; ist++){
        free(dv[ist]);
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(kfunction);
    free(variation);
    free(c_dot);
    free(coef_new);
    free(eenergy);
//    free(na_term);
    free(dv);

}

/*
// Routine for density elec_object scheme in rk4 solver
static void rk4_rho(int nst, int nesteps, double dt, double *energy, double *energy_old,
    double **nacme, double **nacme_old, double **k_lk, double *dotpopd, double complex **rho){

    double complex **k1 = malloc(nst * sizeof(double complex*));
    double complex **k2 = malloc(nst * sizeof(double complex*));
    double complex **k3 = malloc(nst * sizeof(double complex*));
    double complex **k4 = malloc(nst * sizeof(double complex*));
    double complex **kfunction = malloc(nst * sizeof(double complex*));
    double complex **variation = malloc(nst * sizeof(double complex*));
    double complex **rho_dot = malloc(nst * sizeof(double complex*));
    double complex **rho_new = malloc(nst * sizeof(double complex*));
    double *eenergy = malloc(nst * sizeof(double));
//    double *na_term = malloc(nst * sizeof(double));
    double **dv = malloc(nst * sizeof(double*));

    int ist, jst, iestep;
    double frac, edt;

    for(ist = 0; ist < nst; ist++){
        k1[ist] = malloc(nst * sizeof(double complex));
        k2[ist] = malloc(nst * sizeof(double complex));
        k3[ist] = malloc(nst * sizeof(double complex));
        k4[ist] = malloc(nst * sizeof(double complex));
        kfunction[ist] = malloc(nst * sizeof(double complex));
        variation[ist] = malloc(nst * sizeof(double complex));
        rho_dot[ist] = malloc(nst * sizeof(double complex));
        rho_new[ist] = malloc(nst * sizeof(double complex));
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

        // Calculate k1
        ct_rhodot(nst, eenergy, dv, rho, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k1[ist][jst] = edt * rho_dot[ist][jst];
                kfunction[ist][jst] = 0.5 * k1[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        // Calculate k2
        ct_rhodot(nst, eenergy, dv, rho_new, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k2[ist][jst] = edt * rho_dot[ist][jst];
                kfunction[ist][jst] = 0.5 * (- 1.0 + sqrt(2.0)) * k1[ist][jst]
                    + (1.0 - 0.5 * sqrt(2.0)) * k2[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        // Calculate k3
        ct_rhodot(nst, eenergy, dv, rho_new, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k3[ist][jst] = edt * rho_dot[ist][jst];
                kfunction[ist][jst] = - 0.5 * sqrt(2.0) * k2[ist][jst]
                    + (1.0 + 0.5 * sqrt(2.0)) * k3[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        // Calculate k4
        ct_rhodot(nst, eenergy, dv, rho_new, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k4[ist][jst] = edt * rho_dot[ist][jst];
                variation[ist][jst] = (k1[ist][jst] + (2.0 - sqrt(2.0)) * k2[ist][jst]
                    + (2.0 + sqrt(2.0)) * k3[ist][jst] + k4[ist][jst]) / 6.0;
                rho[ist][jst] += variation[ist][jst];
            }
        }

    }

    
    //for(ist = 0; ist < nst; ist++){
    //    na_term[ist] = 0.0;
    //    for(jst = 0; jst < nst; jst++){
    //        if(jst != ist){
    //            na_term[ist] -= 2.0 * dv[ist][jst] * creal(rho[ist][jst]);
    //        }
    //    }
    //    printf("RHODOT_NAC %d %f\n",ist+1, na_term[ist]);
    //}
    //norm = 0.0;
    //for(ist = 0; ist < nst; ist++){
    //    norm += creal(rho[ist][ist]);
    //}
    //printf("RK4_COEF : NORM = %15.8f\n", norm);
    

    for(ist = 0; ist < nst; ist++){
        free(k1[ist]);
        free(k2[ist]);
        free(k3[ist]);
        free(k4[ist]);
        free(kfunction[ist]);
        free(variation[ist]);
        free(rho_dot[ist]);
        free(rho_new[ist]);
        free(dv[ist]);
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(kfunction);
    free(variation);
    free(rho_dot);
    free(rho_new);
    free(eenergy);
//    free(na_term);
    free(dv);

}

*/
