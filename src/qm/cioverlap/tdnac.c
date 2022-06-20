#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)
// DGEMM prototype
extern void dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
    double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
#endif

// Routine to calculate overlap and permutation matrix in MO basis between two time steps
static void calc_MO_over(int nbasis, int norb, double **mo_overlap, double **permut_mat,
    double **ao_overlap, double **mo_coef_old, double **mo_coef_new);

// Routine to match phase of MO coefficients and orderings between two time steps
static void MO_phase_order(int nbasis, int norb, double **mo_coef_new, double **permut_mat);

// Routine to match phase of CI coefficients and orderings between two time steps
static void CI_phase_order(int nst, int norb, int nocc, int nvirt, int *orb_ini, int *orb_final,
    double ***ci_coef_old, double ***ci_coef_new, double **permut_mat);

// Routine to match phase for the states between two time steps
static void state_phase(int nst, int nocc, int nvirt, int *orb_ini, int *orb_final,
    double ***ci_coef_old, double ***ci_coef_new);

// Routine to normalize CI coefficients
static void norm_CI_coef(int nst, int nocc, int nvirt, int *orb_ini, int *orb_final, double ***ci_coef);

// Routine to calculate TDNAC term used in electronic propagation
static void TD_NAC(int istep, int nst, int nbasis, int norb, int nocc, int nvirt, double dt,
    int *orb_ini, int *orb_final, double **nacme, double **ao_overlap, double **mo_coef_old,
    double **mo_coef_new, double ***ci_coef_old, double ***ci_coef_new){

    double **mo_overlap = malloc(norb * sizeof(double*));
    double **permut_mat = malloc(norb * sizeof(double*));

    int ist, jst, ibasis, iorb, jorb, aorb, borb, exponent;
    double fac;
    int debug;

    // This is temporary option to print several variables
    debug = 0;

    for(iorb = 0; iorb < norb; iorb++){
        mo_overlap[iorb] = malloc(norb * sizeof(double));
        permut_mat[iorb] = malloc(norb * sizeof(double));
        for(jorb = 0; jorb < norb; jorb++){
            mo_overlap[iorb][jorb] = 0.0;
            permut_mat[iorb][jorb] = 0.0;
        }
    }

    calc_MO_over(nbasis, norb, mo_overlap, permut_mat, ao_overlap, mo_coef_old, mo_coef_new);

    if(debug == 1){
        // Print mo_overlap
        printf("mo_overlap \n");
        for(iorb = 0; iorb < norb; iorb++){
            for(jorb = 0; jorb < norb; jorb++){
                printf("%15.8f ", mo_overlap[iorb][jorb]);
            }
            printf("\n");
        }
        printf("\n");

        // Print permut_mat
        printf("permut_mat \n");
        for(iorb = 0; iorb < norb; iorb++){
            for(jorb = 0; jorb < norb; jorb++){
                printf("%15.8f ", permut_mat[iorb][jorb]);
            }
            printf("\n");
        }
        printf("\n");

        // Print mo_coef_old
        printf("mo_coef_old \n");
        for(iorb = 0; iorb < norb; iorb++){
            for(ibasis = 0; ibasis < nbasis; ibasis++){
                printf("%15.8f ", mo_coef_old[iorb][ibasis]);
            }
            printf("\n");
        }
        printf("\n");

        // Print mo_coef_new
        printf("mo_coef_new \n");
        for(iorb = 0; iorb < norb; iorb++){
            for(ibasis = 0; ibasis < nbasis; ibasis++){
                printf("%15.8f ", mo_coef_new[iorb][ibasis]);
            }
            printf("\n");
        }
        printf("\n");
    }

    MO_phase_order(nbasis, norb, mo_coef_new, permut_mat);

    if(debug == 1){
        // Print mo_coef_new
        printf("mo_coef_new after phase correction \n");
        for(iorb = 0; iorb < norb; iorb++){
            for(ibasis = 0; ibasis < nbasis; ibasis++){
                printf("%15.8f ", mo_coef_new[iorb][ibasis]);
            }
            printf("\n");
        }
        printf("\n");

        // Print ci_coef_old
        printf("ci_coef_old \n");
        for(iorb = 0; iorb < nocc; iorb++){
            for(aorb = 0; aorb < nvirt; aorb++){
                printf("%15.8f ", ci_coef_old[1][iorb][aorb]);
            }
            printf("\n");
        }
        printf("\n");

        // Print ci_coef_new
        printf("ci_coef_new \n");
        for(iorb = 0; iorb < nocc; iorb++){
            for(aorb = 0; aorb < nvirt; aorb++){
                printf("%15.8f ", ci_coef_new[1][iorb][aorb]);
            }
            printf("\n");
        }
        printf("\n");
    }

    CI_phase_order(nst, norb, nocc, nvirt, orb_ini, orb_final, ci_coef_old, ci_coef_new, permut_mat);

    if(debug == 1){
        // Print ci_coef_new
        printf("ci_coef_new after phase correction \n");
        for(iorb = 0; iorb < nocc; iorb++){
            for(aorb = 0; aorb < nvirt; aorb++){
                printf("%15.8f ", ci_coef_new[1][iorb][aorb]);
            }
            printf("\n");
        }
        printf("\n");
    }

    state_phase(nst, nocc, nvirt, orb_ini, orb_final, ci_coef_old, ci_coef_new);

    if(debug == 1){
        // Print ci_coef_new
        printf("ci_coef_new after state correction \n");
        for(iorb = 0; iorb < nocc; iorb++){
            for(aorb = 0; aorb < nvirt; aorb++){
                printf("%15.8f ", ci_coef_new[1][iorb][aorb]);
            }
            printf("\n");
        }
        printf("\n");
    }

    if(istep == 0){
        norm_CI_coef(nst, nocc, nvirt, orb_ini, orb_final, ci_coef_old);
    }
    norm_CI_coef(nst, nocc, nvirt, orb_ini, orb_final, ci_coef_new);

    // Re-calculate mo_overlap with phase-corrected MO coefficients
    // Now, mo_overlap is anti-symmetric
    for(iorb = 0; iorb < norb; iorb++){
        for(jorb = 0; jorb < norb; jorb++){
            mo_overlap[iorb][jorb] = 0.0;
        }
    }
    calc_MO_over(nbasis, norb, mo_overlap, permut_mat, ao_overlap, mo_coef_old, mo_coef_new);

    if(debug == 1){
        // Print mo_overlap
        printf("mo_overlap with phase-corrected mo \n");
        for(iorb = 0; iorb < norb; iorb++){
            for(jorb = 0; jorb < norb; jorb++){
                printf("%15.8f ", mo_overlap[iorb][jorb]);
            }
            printf("\n");
        }
        printf("\n");

        // Print permut_mat
        printf("permut_mat with phase-corrected mo \n");
        for(iorb = 0; iorb < norb; iorb++){
            for(jorb = 0; jorb < norb; jorb++){
                printf("%15.8f ", permut_mat[iorb][jorb]);
            }
            printf("\n");
        }
        printf("\n");
    }

    // TODO : ist = jst should be removed
    // TODO : we need to evaluate only ist < jst cases since NACME is anti-symmetric
    for(ist = 0; ist < nst; ist++){
        for(jst = 0; jst < nst; jst++){

            if(ist == 0 || jst == 0){

                if(ist > jst){

                    // TDNAC between S_i and S_0 state
                    for(iorb = orb_ini[0]; iorb < nocc; iorb++){
                        for(aorb = 0; aorb < orb_final[0] - nocc; aorb++){
                            nacme[ist][jst] += 0.5 * ci_coef_new[ist][iorb][aorb] * (mo_overlap[nocc + aorb][iorb] - mo_overlap[iorb][nocc + aorb]);
                        }
                    }

                }
                else{

                    // TDNAC between S_0 and S_j state
                    for(jorb = orb_ini[0]; jorb < nocc; jorb++){
                        for(borb = 0; borb < orb_final[0] - nocc; borb++){
                            nacme[ist][jst] += 0.5 * ci_coef_new[jst][jorb][borb] * (mo_overlap[jorb][nocc + borb] - mo_overlap[nocc + borb][jorb]);
                        }
                    }

                }

            }
            else{

                // TDNAC between S_i and S_j state

                // 1st term in Eq. 15
                for(iorb = orb_ini[0]; iorb < nocc; iorb++){
                    for(aorb = 0; aorb < orb_final[0] - nocc; aorb++){
                        nacme[ist][jst] += 0.5 * (ci_coef_old[ist][iorb][aorb] * ci_coef_new[jst][iorb][aorb] - ci_coef_old[jst][iorb][aorb] * ci_coef_new[ist][iorb][aorb]);
                    }
                }

                // 2nd term in Eq. 15
                for(iorb = orb_ini[0]; iorb < nocc; iorb++){
                    for(aorb = 0; aorb < orb_final[0] - nocc; aorb++){
                        for(borb = 0; borb < orb_final[0] - nocc; borb++){
                            if(aorb != borb){
                                nacme[ist][jst] += 0.5 * ci_coef_new[ist][iorb][aorb] * ci_coef_new[jst][iorb][borb] * (mo_overlap[nocc + aorb][nocc + borb] - mo_overlap[nocc + borb][nocc + aorb]);
                            }
                        }
                    }
                }

                // 3rd term in Eq. 15
                for(iorb = orb_ini[0]; iorb < nocc; iorb++){
                    for(aorb = 0; aorb < orb_final[0] - nocc; aorb++){
                        for(jorb = orb_ini[0]; jorb < nocc; jorb++){
                            if(iorb != jorb){
                                // fac is permutation in 3rd term
                                exponent = abs(jorb - iorb);
                                fac = pow(-1.0, exponent);
                                nacme[ist][jst] -= 0.5 * fac * ci_coef_new[ist][iorb][aorb] * ci_coef_new[jst][jorb][aorb] * (mo_overlap[jorb][iorb] - mo_overlap[iorb][jorb]);
                            }
                        }
                    }
                }

            }

            // NACME is divided by time step (finite numerical differentiation)
            nacme[ist][jst] /= dt;

        }
    }

    if(debug == 1){
        // Print NACME values
        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                printf("%15.8f ", nacme[ist][jst]);
            }
            printf("\n");
        }
        printf("\n");
    }

    for(iorb = 0; iorb < norb; iorb++){
        free(mo_overlap[iorb]);
        free(permut_mat[iorb]);
    }

    free(mo_overlap);
    free(permut_mat);

}

// Routine to calculate overlap and permutation matrix in MO basis between two time steps
static void calc_MO_over(int nbasis, int norb, double **mo_overlap, double **permut_mat,
    double **ao_overlap, double **mo_coef_old, double **mo_coef_new){

    double **tmp_sign = malloc(norb * sizeof(double*));
    int ibasis, jbasis, iorb, jorb;

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)
    double *mo_overlap_1d = malloc((norb * norb) * sizeof(double));
    double *ao_overlap_1d = malloc((nbasis * nbasis) * sizeof(double));
    double *mo_coef_1d = malloc((norb * nbasis) * sizeof(double));
    double *tmp_mat_1d = malloc((norb * nbasis) * sizeof(double));

    double one = 1.0;
    double zero = 0.0;
#endif

    // Initialize temporary array to save sign of overlap in MO basis
    for(iorb = 0; iorb < norb; iorb++){
        tmp_sign[iorb] = malloc(norb * sizeof(double));
        for(jorb = 0; jorb < norb; jorb++){
            tmp_sign[iorb][jorb] = 0.0;
        }
    }

    // Calculate overlap in MO basis using external libraries or direct loop calculation; S' = C * S * C^T
#if defined(HAVE_LAPACK) || defined(HAVE_MKL)
    // Convert mo_coef_old to column-major 1D array
    for(iorb = 0; iorb < norb; iorb++){
        for(ibasis = 0; ibasis < nbasis; ibasis++){
//            mo_coef_1d[iorb * nbasis + ibasis] = mo_coef_old[iorb][ibasis];
            mo_coef_1d[ibasis * norb + iorb] = mo_coef_old[iorb][ibasis];
        }
    }
    // Convert ao_overlap to column-major 1D array
    for(ibasis = 0; ibasis < nbasis; ibasis++){
        for(jbasis = 0; jbasis < nbasis; jbasis++){
//            ao_overlap_1d[ibasis * nbasis + jbasis] = ao_overlap[ibasis][jbasis];
            ao_overlap_1d[jbasis * nbasis + ibasis] = ao_overlap[ibasis][jbasis];
        }
    }

    // Execute following operation; C * S
    dgemm_("N", "N", &norb, &nbasis, &nbasis, &one, mo_coef_1d, &norb, ao_overlap_1d, &nbasis, &zero, tmp_mat_1d, &norb);

    // Convert mo_coef_new to column-major 1D array
    for(iorb = 0; iorb < norb; iorb++){
        for(ibasis = 0; ibasis < nbasis; ibasis++){
//            mo_coef_1d[iorb * nbasis + ibasis] = mo_coef_new[iorb][ibasis];
            mo_coef_1d[ibasis * norb + iorb] = mo_coef_new[iorb][ibasis];
        }
    }

    // Execute following operation; (C * S) * C^T
    dgemm_("N", "T", &norb, &norb, &nbasis, &one, tmp_mat_1d, &norb, mo_coef_1d, &norb, &zero, mo_overlap_1d, &norb);

    for(iorb = 0; iorb < norb; iorb++){
        for(jorb = 0; jorb < norb; jorb++){

            // Convert column-major mo_overlap_1d to 2D array
//            mo_overlap[iorb][jorb] = mo_overlap_1d[iorb * norb + jorb];
            mo_overlap[iorb][jorb] = mo_overlap_1d[jorb * norb + iorb];

            // Save the sign of overlap elements for calculating permutation matrix
            if(mo_overlap[iorb][jorb] < 0.0){
                tmp_sign[iorb][jorb] = -1.0;
            }
            else{
                tmp_sign[iorb][jorb] = 1.0;
            }

        }
    }
#else
    for(iorb = 0; iorb < norb; iorb++){
        for(jorb = 0; jorb < norb; jorb++){

            for(ibasis = 0; ibasis < nbasis; ibasis++){
                for(jbasis = 0; jbasis < nbasis; jbasis++){
                    mo_overlap[iorb][jorb] += mo_coef_old[iorb][ibasis] * ao_overlap[ibasis][jbasis] * mo_coef_new[jorb][jbasis];
                    // Save the sign of overlap elements for calculating permutation matrix
                    if(mo_overlap[iorb][jorb] < 0.0){
                        tmp_sign[iorb][jorb] = -1.0;
                    }
                    else{
                        tmp_sign[iorb][jorb] = 1.0;
                    }
                }
            }

        }
    }
#endif

    for(iorb = 0; iorb < norb; iorb++){
        for(jorb = 0; jorb < norb; jorb++){
            // Permutation matrix is obtained by rounding off the square of overlap matrix in MO basis
            permut_mat[iorb][jorb] = round(pow(mo_overlap[iorb][jorb], 2) * tmp_sign[iorb][jorb]);
        }
    }

    // Deallocate temporary array to save sign of overlap in MO basis
    for(iorb = 0; iorb < norb; iorb++){
        free(tmp_sign[iorb]);
    }

    free(tmp_sign);

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)
    free(mo_overlap_1d);
    free(ao_overlap_1d);
    free(mo_coef_1d);
    free(tmp_mat_1d);
#endif

}

// Routine to match phase of MO coefficients and orderings between two time steps
static void MO_phase_order(int nbasis, int norb, double **mo_coef_new, double **permut_mat){

    double **tmp_mo = malloc(norb * sizeof(double*));
    int ibasis, iorb, jorb;

    // Initialize temporary MO array; C'
    for(iorb = 0; iorb < norb; iorb++){
        tmp_mo[iorb] = malloc(nbasis * sizeof(double));
        for(ibasis = 0; ibasis < nbasis; ibasis++){
            tmp_mo[iorb][ibasis] = 0.0;
        }
    }

    // Decide the phase and ordering for MO coefficients using permutation matrix; C' = O * C
    for(iorb = 0; iorb < norb; iorb++){
        for(ibasis = 0; ibasis < nbasis; ibasis++){

            for(jorb = 0; jorb < norb; jorb++){
                tmp_mo[iorb][ibasis] += permut_mat[iorb][jorb] * mo_coef_new[jorb][ibasis];
            }

        }
    }

    // Apply new phase correction for the MO coefficients; C = C'
    for(iorb = 0; iorb < norb; iorb++){
        for(ibasis = 0; ibasis < nbasis; ibasis++){
            mo_coef_new[iorb][ibasis] = tmp_mo[iorb][ibasis];
        }
    }

    // Deallocate temporary MO array; C'
    for(iorb = 0; iorb < norb; iorb++){
        free(tmp_mo[iorb]);
    }

    free(tmp_mo);

}

// Routine to match phase of CI coefficients and orderings between two time steps
// TODO : Is this correct method to match phase (or order) for CI coefficients?
static void CI_phase_order(int nst, int norb, int nocc, int nvirt, int *orb_ini, int *orb_final,
    double ***ci_coef_old, double ***ci_coef_new, double **permut_mat){

    double **tmp_ci = malloc(norb * sizeof(double*));
    double **tmp_ci_new = malloc(norb * sizeof(double*));
    int ist, iorb, jorb, aorb, borb;

    // Initialize temporary CI arrays; C' and C
    for(iorb = 0; iorb < norb; iorb++){
        tmp_ci[iorb] = malloc(norb * sizeof(double));
        tmp_ci_new[iorb] = malloc(norb * sizeof(double));
    }

    // CI coefficients for S_0 are zero
    for(ist = 1; ist < nst; ist++){

        for(iorb = orb_ini[0]; iorb < orb_final[0]; iorb++){
            for(aorb = orb_ini[0]; aorb < orb_final[0]; aorb++){
                // Assign CI coefficients at time t to new symmetric array
                if(iorb < nocc && aorb >= nocc){
                    tmp_ci[iorb][aorb] = ci_coef_new[ist][iorb][aorb - nocc];
                }
                else if(iorb >= nocc && aorb < nocc){
                    tmp_ci[iorb][aorb] = ci_coef_new[ist][aorb][iorb - nocc];
                }
                else{
                    tmp_ci[iorb][aorb] = 0.0;
                }
                // Initialize new empty array for phase correction
                tmp_ci_new[iorb][aorb] = 0.0;
            }
        }

        // Decide the phase and ordering for CI coefficients using permutation matrix; C' = O * C * O
        // TODO : The phases for occupied and virtual orbitals are matched when permutation is diagonal matrix
        for(jorb = orb_ini[0]; jorb < orb_final[0]; jorb++){
            for(borb = orb_ini[0]; borb < orb_final[0]; borb++){

                for(iorb = orb_ini[0]; iorb < orb_final[0]; iorb++){
                    for(aorb = orb_ini[0]; aorb < orb_final[0]; aorb++){
                        tmp_ci_new[jorb][borb] += permut_mat[jorb][iorb] * tmp_ci[iorb][aorb] * permut_mat[aorb][borb];
                    }
                }

            }
        }

        // Apply new phase correction for the CI coefficients; C = C'
        for(iorb = orb_ini[0]; iorb < nocc; iorb++){
            for(aorb = 0; aorb < orb_final[0] - nocc; aorb++){
                ci_coef_new[ist][iorb][aorb] = tmp_ci_new[iorb][nocc + aorb];
            }
        }

    }

    // Deallocate temporary CI arrays; C' and C
    for(iorb = 0; iorb < norb; iorb++){
        free(tmp_ci[iorb]);
        free(tmp_ci_new[iorb]);
    }

    free(tmp_ci);
    free(tmp_ci_new);

}

// Routine to match phase for the states between two time steps
static void state_phase(int nst, int nocc, int nvirt, int *orb_ini, int *orb_final,
    double ***ci_coef_old, double ***ci_coef_new){

    double val;
    int ist, iorb, aorb;

    // CI coefficients for S_0 are zero
    for(ist = 1; ist < nst; ist++){

        val = 0.0;
        for(iorb = orb_ini[0]; iorb < nocc; iorb++){
            for(aorb = 0; aorb < orb_final[0] - nocc; aorb++){
                val += ci_coef_old[ist][iorb][aorb] * ci_coef_new[ist][iorb][aorb];
            }
        }

        if(val < 0.0){
            for(iorb = orb_ini[0]; iorb < nocc; iorb++){
                for(aorb = 0; aorb < orb_final[0] - nocc; aorb++){
                    ci_coef_new[ist][iorb][aorb] *= -1.0;
                }
            }
        }

    }

}

// Routine to normalize CI coefficients
static void norm_CI_coef(int nst, int nocc, int nvirt, int *orb_ini, int *orb_final, double ***ci_coef){

    double norm;
    int ist, iorb, aorb;

    // CI coefficients for S_0 are zero
    for(ist = 1; ist < nst; ist++){

        // Calculate normalization value for CI coefficients
        norm = 0.0;
        for(iorb = orb_ini[0]; iorb < nocc; iorb++){
            for(aorb = 0; aorb < orb_final[0] - nocc; aorb++){
                norm += pow(ci_coef[ist][iorb][aorb], 2);
            }
        }
        norm = sqrt(norm);

        // Normalize the CI coefficients
        for(iorb = orb_ini[0]; iorb < nocc; iorb++){
            for(aorb = 0; aorb < orb_final[0] - nocc; aorb++){
                ci_coef[ist][iorb][aorb] /= norm;
            }
        }

    }

}

