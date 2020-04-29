#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

//double dot(double complex* u, double complex* v, int nst);
//void cdot(int nst, double complex* c, double* e, double **dv, double complex *c_dot);

/*static void rhodot(int nst, double complex **rho, double *e, double **dv, double complex **rho_dot){
    
    int ist, jst, kst;
    for(ist = 0; ist < nst; ist++){
        for(jst = 0; jst < nst; jst++){
            rho_dot[ist][jst] = 0.0 + 0.0 * I;
        }
    }

    for(ist = 0; ist < nst; ist++){
        for(jst = 0; jst < nst; jst++){
            if(ist != jst){
                rho_dot[ist][ist] -= dv[ist][jst] * 2.0 * creal(rho[ist][jst]);
            }
        }
    }

    for(ist = 0; ist < nst; ist++){
        for(jst = ist + 1; jst < nst; jst++){
            rho_dot[ist][jst] -=  1.0 * I * (e[jst] - e[ist]) * rho[ist][jst];

            for(kst = 0; kst < nst; kst++){
                rho_dot[ist][jst] -= dv[ist][kst] * rho[kst][jst] + dv[jst][kst] * rho[ist][kst];
            }
            rho_dot[jst][ist] = conj(rho_dot[ist][jst]);
        }
    }
}*/

static void TD_NAC(int nst, int norb, int nocc, int nvirt, double **nacme, double **ao_overlap, double **mo_coef_old, double **mo_coef_new, double ***ci_coef_old, double ***ci_coef_new){
    double **mo_overlap = malloc(norb * sizeof(double*));
    int ist, iorb, jorb, aorb, borb, mu, nu;
//    double frac, edt, norm;

//    printf("TD_NAC : nst = %8d\n", nst);
//    printf("TD_NAC : norb = %ld\n", sizeof(ao_overlap) / sizeof(ao_overlap[0]));
//    printf("TD_NAC : norb = %ld\n", sizeof(tmp)); // sizeof(ao_overlap[0]));
//    printf("TD_NAC : ao_overlap = %15.8f\n", ao_overlap[0][0]);

//    frac = 1.0 / (double)nesteps;
//    edt = dt * frac;

    for(iorb = 0; iorb < norb; iorb++){
        mo_overlap[iorb] = malloc(norb * sizeof(double));
        // Initialize mo_overlap
        for(jorb = 0; jorb < norb; jorb++){
            mo_overlap[iorb][jorb] = 0.0;
        }
    }

    for(iorb = 0; iorb < norb; iorb++){
        for(jorb = 0; jorb < norb; jorb++){

            for(mu = 0; mu < norb; mu++){
                for(nu = 0; nu < norb; nu++){
                    mo_overlap[iorb][jorb] += ao_overlap[mu][nu] * mo_coef_old[iorb][mu] * mo_coef_new[jorb][nu];
                }
            }

        }
    }
//    // print
//    for(iorb = 0; iorb < norb; iorb++){
//        for(jorb = 0; jorb < norb; jorb++){
//            printf("%15.8f ", mo_overlap[iorb][jorb]);
//        }
//        printf("\n");
//    }
//    printf("\n");
//    printf("\n");

    // TODO : ist = jst should be removed
    for(ist = 0; ist < nst; ist++){
        for(jst = 0; jst < nst; jst++){

            // 1st term in Eq. 15
            for(iorb = 0; iorb < nocc; iorb++){
                for(aorb = 0; aorb < nvirt; aorb++){
                    nacme[ist][jst] += ci_coef_old[ist][iorb][aorb] * ci_coef_new[jst][iorb][aorb];
                }
            }

            // 2nd term in Eq. 15
            for(iorb = 0; iorb < nocc; iorb++){
                for(aorb = 0; aorb < nvirt; aorb++){
                    for(borb = 0; borb < nvirt; aorb++){
                        nacme[ist][jst] += ci_coef_old[ist][iorb][aorb] * ci_coef_new[jst][iorb][borb] * mo_overlap[nocc + aorb][nocc + borb];
                    }
                }
            }

            // 3rd term in Eq. 15
            for(iorb = 0; iorb < nocc; iorb++){
                for(aorb = 0; aorb < nvirt; aorb++){
                    for(jorb = 0; jorb < nocc; jorb++){
                        // TODO : permutation needed
                        nacme[ist][jst] -= ci_coef_old[ist][iorb][aorb] * ci_coef_new[jst][jorb][aorb] * mo_overlap[jorb][iorb];
                    }
                }
            }

            // TODO : dt should be obtained from md.dt
//            nacme[ist][jst] /= dt;

        }
    }
    
/*    for(iestep = 0; iestep < nesteps; iestep++){
        
        for(ist = 0; ist < nst; ist++){
            eenergy[ist] = energy_old[ist] + (energy[ist] - energy_old[ist]) *(double) iestep * frac;
            for(jst = 0; jst < nst; jst++){
                dv[ist][jst] = nacme_old[ist][jst] + (nacme[ist][jst] - nacme_old[ist][jst]) *(double) iestep * frac;
            }
        }

        rhodot(nst, rho, eenergy, dv, rho_dot);

    }
    
    norm = 0.0;
    for(ist = 0; ist < nst; ist++){
        norm += creal(rho[ist][ist]);
    }

    */

    for(iorb = 0; iorb < norb; iorb++){
        free(mo_overlap[iorb]);
    }

    free(mo_overlap);

}


