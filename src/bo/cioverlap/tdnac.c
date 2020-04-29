#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

// Routine to calculate overlap and permutation matrix in MO basis between two time steps
static void calc_MO_over(int norb, double **mo_overlap, double **permut_mat, double **ao_overlap, double **mo_coef_old, double **mo_coef_new){

    int iorb, jorb, mu, nu;

    for(iorb = 0; iorb < norb; iorb++){
        for(jorb = 0; jorb < norb; jorb++){

            for(mu = 0; mu < norb; mu++){
                for(nu = 0; nu < norb; nu++){
                    mo_overlap[iorb][jorb] += ao_overlap[mu][nu] * mo_coef_old[iorb][mu] * mo_coef_new[jorb][nu];
                }
            }
            // Permutation matrix is obtained by rounding off the overlap matrix in MO basis
            permut_mat[iorb][jorb] = round(mo_overlap[iorb][jorb]);

        }
    }

    // Print mo_overlap and permut_mat values
//    for(iorb = 0; iorb < norb; iorb++){
//        for(jorb = 0; jorb < norb; jorb++){
//            printf("%15.8f ", mo_overlap[iorb][jorb]);
//            printf("%15.8f ", permut_mat[iorb][jorb]);
//        }
//        printf("\n");
//    }
//    printf("\n");

}

// Routine to match phase of MO coefficients and orderings between two time steps
static void MO_phase_order(int norb, double **mo_coef_new, double **permut_mat){

    double **tmp_mo = malloc(norb * sizeof(double*));
    int iorb, jorb, mu;

    for(iorb = 0; iorb < norb; iorb++){
        tmp_mo[iorb] = malloc(norb * sizeof(double));
        for(mu = 0; mu < norb; mu++){
            tmp_mo[iorb][mu] = 0.0;
        }
    }

    for(iorb = 0; iorb < norb; iorb++){
        for(mu = 0; mu < norb; mu++){

            // Decide the phase and ordering for MO coefficients
            for(jorb = 0; jorb < norb; jorb++){
                tmp_mo[iorb][mu] += permut_mat[iorb][jorb] * mo_coef_new[jorb][mu];
            }

        }
    }

    // Print mo_coef_new values before matching and ordering
//    for(iorb = 0; iorb < norb; iorb++){
//        for(jorb = 0; jorb < norb; jorb++){
//            printf("%15.8f ", mo_coef_new[iorb][jorb]);
//        }
//        printf("\n");
//    }
//    printf("\n");

    for(iorb = 0; iorb < norb; iorb++){
        for(mu = 0; mu < norb; mu++){
            mo_coef_new[iorb][mu] = tmp_mo[iorb][mu];
        }
    }

    // Print mo_coef_new values after matching and ordering
//    for(iorb = 0; iorb < norb; iorb++){
//        for(jorb = 0; jorb < norb; jorb++){
//            printf("%15.8f ", mo_coef_new[iorb][jorb]);
//        }
//        printf("\n");
//    }
//    printf("\n");

    for(iorb = 0; iorb < norb; iorb++){
        free(tmp_mo[iorb]);
    }

    free(tmp_mo);

}

// Routine to match phase of CI coefficients and orderings between two time steps
// TODO : Is this correct method to match phase (or order) for CI coefficients?
static void CI_phase_order(int nst, int nocc, int nvirt, double ***ci_coef_old, double ***ci_coef_new){

    double val;
    int ist, iorb, aorb;

    for(ist = 0; ist < nst; ist++){
        for(aorb = 0; aorb < nvirt; aorb++){

            val = 0.0;
            for(iorb = 0; iorb < nocc; iorb++){
                val += ci_coef_old[ist][iorb][aorb] * ci_coef_new[ist][iorb][aorb];
            }

            if(val < 0.0){
                for(iorb = 0; iorb < nocc; iorb++){
                    ci_coef_new[ist][iorb][aorb] *= -1.0;
                }
            }

        }
    }

    // Print ci_coef_old and ci_coef_new values
//    for(iorb = 0; iorb < nocc; iorb++){
//        for(aorb = 0; aorb < nvirt; aorb++){
//            printf("%15.8f ", ci_coef_old[0][iorb][aorb]);
//        }
//        printf("\n");
//    }
//    printf("\n");

//    for(iorb = 0; iorb < nocc; iorb++){
//        for(aorb = 0; aorb < nvirt; aorb++){
//            printf("%15.8f ", ci_coef_new[0][iorb][aorb]);
//        }
//        printf("\n");
//    }
//    printf("\n");

}

// Routine to normalize CI coefficients
static void norm_CI_coef(int nst, int nocc, int nvirt, double ***ci_coef){

    double norm;
    int ist, iorb, aorb;

    for(ist = 0; ist < nst; ist++){

        // Calculate norm value for CI coefficients
        norm = 0.0;
        for(iorb = 0; iorb < nocc; iorb++){
            for(aorb = 0; aorb < nvirt; aorb++){
                norm += pow(ci_coef[ist][iorb][aorb], 2);
            }
        }
        norm = sqrt(norm);

        // Normalize the CI coefficients
        for(iorb = 0; iorb < nocc; iorb++){
            for(aorb = 0; aorb < nvirt; aorb++){
                ci_coef[ist][iorb][aorb] /= norm;
            }
        }

        // Print the CI coefficients
//        if(ist == 0){
//            for(iorb = 0; iorb < nocc; iorb++){
//                for(aorb = 0; aorb < nvirt; aorb++){
//                    printf("%15.8f ", ci_coef[ist][iorb][aorb]);
//                }
//                printf("\n");
//            }
//            printf("\n");
//        }

    }

}

// Routine to calculate TDNAC term used in electronic propagation
static void TD_NAC(int nst, int norb, int nocc, int nvirt, double **nacme, double **ao_overlap, double **mo_coef_old, double **mo_coef_new, double ***ci_coef_old, double ***ci_coef_new){
    double **mo_overlap = malloc(norb * sizeof(double*));
    double **permut_mat = malloc(norb * sizeof(double*));
    int ist, jst, iorb, jorb, aorb, borb;
//    double frac, edt, norm;

//    printf("TD_NAC : nst = %8d\n", nst);
//    printf("TD_NAC : norb = %ld\n", sizeof(ao_overlap) / sizeof(ao_overlap[0]));
//    printf("TD_NAC : norb = %ld\n", sizeof(tmp)); // sizeof(ao_overlap[0]));
//    printf("TD_NAC : ao_overlap = %15.8f\n", ao_overlap[0][0]);

//    frac = 1.0 / (double)nesteps;
//    edt = dt * frac;

    for(iorb = 0; iorb < norb; iorb++){
        mo_overlap[iorb] = malloc(norb * sizeof(double));
        permut_mat[iorb] = malloc(norb * sizeof(double));
        // Initialize mo_overlap
        for(jorb = 0; jorb < norb; jorb++){
            mo_overlap[iorb][jorb] = 0.0;
            permut_mat[iorb][jorb] = 0.0;
        }
    }

    calc_MO_over(norb, mo_overlap, permut_mat, ao_overlap, mo_coef_old, mo_coef_new);

    MO_phase_order(norb, mo_coef_new, permut_mat);

    CI_phase_order(nst, nocc, nvirt, ci_coef_old, ci_coef_new);

    norm_CI_coef(nst, nocc, nvirt, ci_coef_old);
    norm_CI_coef(nst, nocc, nvirt, ci_coef_new);

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
                    for(borb = 0; borb < nvirt; borb++){
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

            // print NACME values
            printf("%15.8f ", nacme[ist][jst]);

        }
        printf("\n");
    }
    printf("\n");
    
    for(iorb = 0; iorb < norb; iorb++){
        free(mo_overlap[iorb]);
        free(permut_mat[iorb]);
    }

    free(mo_overlap);
    free(permut_mat);

}


