#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Routine to calculate overlap and permutation matrix in MO basis between two time steps
static void calc_MO_over(int nbasis, int norb, int nspin, double ***mo_overlap, double ***permut_mat,
    double **ao_overlap, double ***mo_coef_old, double ***mo_coef_new);

// Routine to match phase of MO coefficients and orderings between two time steps
static void MO_phase_order(int nbasis, int norb, int nspin, double ***mo_coef_new, double ***permut_mat);

// Routine to match phase of CI coefficients and orderings between two time steps
static void CI_phase_order(int nst, int norb, int nspin, int *nocc, int *nvirt, int *orb_ini, int *orb_final,
    double ****ci_coef_old, double ****ci_coef_new, double ***permut_mat);

// Routine to match phase for the states between two time steps
static void state_phase(int nst, int nspin, int *nocc, int *nvirt, int *orb_ini, int *orb_final,
    double ****ci_coef_old, double ****ci_coef_new);

// Routine to normalize CI coefficients
static void norm_CI_coef(int nst, int nspin, int *nocc, int *nvirt, int *orb_ini, int *orb_final, double ****ci_coef);

// Routine to calculate TDNAC term used in electronic propagation
static void TD_NAC(int istep, int nst, int nbasis, int nspin, int norb, int *nocc, int *nvirt, double dt,
    int *orb_ini, int *orb_final, double **nacme, double **ao_overlap, double ***mo_coef_old,
    double ***mo_coef_new, double ****ci_coef_old, double ****ci_coef_new){

    double ***mo_overlap = malloc(norb * sizeof(double**));
    double ***permut_mat = malloc(norb * sizeof(double**));

    int ist, jst, nsp, nocc_max, nvirt_max, ibasis, iorb, jorb, aorb, borb, exponent;
    double fac;
    int debug;

    // This is temporary option to print several variables
    debug = 0;

    for(iorb = 0; iorb < norb; iorb++){
        mo_overlap[iorb] = malloc(norb * sizeof(double*));
        permut_mat[iorb] = malloc(norb * sizeof(double*));
        for(jorb = 0; jorb < norb; jorb++){
            mo_overlap[iorb][jorb] = malloc(nspin * sizeof(double));
            permut_mat[iorb][jorb] = malloc(nspin * sizeof(double));
	    	for(nsp = 0; nsp < nspin; nsp++){
            		mo_overlap[iorb][jorb][nsp] = 0.0;
            		permut_mat[iorb][jorb][nsp] = 0.0;
		}
        }
    }

    calc_MO_over(nbasis, norb, nspin, mo_overlap, permut_mat, ao_overlap, mo_coef_old, mo_coef_new);

    if(debug == 1){
	    // NOT overworked for spin polarization
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

    MO_phase_order(nbasis, norb, nspin, mo_coef_new, permut_mat);

    if(debug == 1){
	    // NOT overworked for spin polarization
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

    CI_phase_order(nst, norb, nspin, nocc, nvirt, orb_ini, orb_final, ci_coef_old, ci_coef_new, permut_mat);

    if(debug == 1){
	    // NOT overworked for spin polarization
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

    state_phase(nst, nspin, nocc, nvirt, orb_ini, orb_final, ci_coef_old, ci_coef_new);

    if(debug == 1){
	    // NOT overworked for spin polarization
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
        norm_CI_coef(nst, nspin, nocc, nvirt, orb_ini, orb_final, ci_coef_old);
    }
    norm_CI_coef(nst, nspin, nocc, nvirt, orb_ini, orb_final, ci_coef_new);

    // Re-calculate mo_overlap with phase-corrected MO coefficients
    // Now, mo_overlap is anti-symmetric
    for(iorb = 0; iorb < norb; iorb++){
        for(jorb = 0; jorb < norb; jorb++){
        	for(nsp = 0; nsp < nspin; nsp++){
            		mo_overlap[iorb][jorb][nsp] = 0.0;
		}
        }
    }
    calc_MO_over(nbasis, norb, nspin, mo_overlap, permut_mat, ao_overlap, mo_coef_old, mo_coef_new);

    if(debug == 1){
	    // NOT overworked for spin polarization
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
		    for(nsp = 0; nsp < nspin ; nsp++){
                    	for(iorb = orb_ini[nsp]; iorb < nocc[nsp]; iorb++){
                        	for(aorb = 0; aorb < orb_final[nsp] - nocc[nsp]; aorb++){
                            		nacme[ist][jst] += 0.5 * ci_coef_new[ist][iorb][aorb][nsp] * (mo_overlap[nocc[nsp] + aorb][iorb][nsp] - mo_overlap[iorb][nocc[nsp] + aorb][nsp]);
                        }
                    	}
		    }

                }
                else{

                    // TDNAC between S_0 and S_j state
		    for(nsp = 0; nsp < nspin ; nsp++){
                    	for(jorb = orb_ini[nsp]; jorb < nocc[nsp]; jorb++){
                        	for(borb = 0; borb < orb_final[nsp] - nocc[nsp]; borb++){
                            	nacme[ist][jst] += 0.5 * ci_coef_new[jst][jorb][borb][nsp] * (mo_overlap[jorb][nocc[nsp] + borb][nsp] - mo_overlap[nocc[nsp] + borb][jorb][nsp]);
                        	}
                    	}
		    }

                }

            }
            else{

                // TDNAC between S_i and S_j state

                // 1st term in Eq. 15
		
		for(nsp = 0; nsp < nspin ; nsp++){
                for(iorb = orb_ini[nsp]; iorb < nocc[nsp]; iorb++){
                    for(aorb = 0; aorb < orb_final[nsp] - nocc[nsp]; aorb++){
                        nacme[ist][jst] += 0.5 * (ci_coef_old[ist][iorb][aorb][nsp] * ci_coef_new[jst][iorb][aorb][nsp] - ci_coef_old[jst][iorb][aorb][nsp] * ci_coef_new[ist][iorb][aorb][nsp]);
                    }
                }
		}

                // 2nd term in Eq. 15
		for(nsp = 0; nsp < nspin ; nsp++){
                for(iorb = orb_ini[nsp]; iorb < nocc[nsp]; iorb++){
                    for(aorb = 0; aorb < orb_final[nsp] - nocc[nsp]; aorb++){
                        for(borb = 0; borb < orb_final[nsp] - nocc[nsp]; borb++){
                            if(aorb != borb){
                                nacme[ist][jst] += 0.5 * ci_coef_new[ist][iorb][aorb][nsp] * ci_coef_new[jst][iorb][borb][nsp] * (mo_overlap[nocc[nsp] + aorb][nocc[nsp] + borb][nsp] - mo_overlap[nocc[nsp] + borb][nocc[nsp] + aorb][nsp]);
                            }
                        }
                    }
                }
		}

                // 3rd term in Eq. 15
		for(nsp = 0; nsp < nspin ; nsp++){
                for(iorb = orb_ini[nsp]; iorb < nocc[nsp]; iorb++){
                    for(aorb = 0; aorb < orb_final[nsp] - nocc[nsp]; aorb++){
                        for(jorb = orb_ini[nsp]; jorb < nocc[nsp]; jorb++){
                            if(iorb != jorb){
                                // fac is permutation in 3rd term
                                exponent = abs(jorb - iorb);
                                fac = pow(-1.0, exponent);
                                nacme[ist][jst] -= 0.5 * fac * ci_coef_new[ist][iorb][aorb][nsp] * ci_coef_new[jst][jorb][aorb][nsp] * (mo_overlap[jorb][iorb][nsp] - mo_overlap[iorb][jorb][nsp]);
                            }
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
	    // NOT overworked for spin polarization
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
    for(jorb = 0; jorb < norb; jorb++){
        free(mo_overlap[iorb][jorb]);
        free(permut_mat[iorb][jorb]);
    }
    }
    for(jorb = 0; jorb < norb; jorb++){
        free(mo_overlap[jorb]);
        free(permut_mat[jorb]);
    }

    free(mo_overlap);
    free(permut_mat);

}

// Routine to calculate overlap and permutation matrix in MO basis between two time steps
static void calc_MO_over(int nbasis, int norb, int nspin, double ***mo_overlap, double ***permut_mat,
    double **ao_overlap, double ***mo_coef_old, double ***mo_coef_new){

    double ***tmp_sign = malloc(norb * sizeof(double**));
    int ibasis, jbasis, iorb, jorb, nsp;

    // Initialize temporary array to save sign of overlap in MO basis
    for(iorb = 0; iorb < norb; iorb++){
        tmp_sign[iorb] = malloc(norb * sizeof(double*));
        for(jorb = 0; jorb < norb; jorb++){
        tmp_sign[iorb][jorb] = malloc(nspin * sizeof(double*));
        	for(nsp = 0; nsp < nspin; nsp++){
            	tmp_sign[iorb][jorb][nsp] = 0.0;
        }
	}
    }
    
    
    for(nsp = 0; nsp < nspin; nsp++){
    for(iorb = 0; iorb < norb; iorb++){
        for(jorb = 0; jorb < norb; jorb++){

            // Calculate overlap in MO basis; S' = C * S * C^T
            for(ibasis = 0; ibasis < nbasis; ibasis++){
                for(jbasis = 0; jbasis < nbasis; jbasis++){
                    mo_overlap[iorb][jorb][nsp] += mo_coef_old[iorb][ibasis][nsp] * ao_overlap[ibasis][jbasis] * mo_coef_new[jorb][jbasis][nsp];
                    // Save the sign of overlap elements for calculating permutation matrix
                    if(mo_overlap[iorb][jorb][nsp] < 0.0){
                        tmp_sign[iorb][jorb][nsp] = -1.0;
                    }
                    else{
                        tmp_sign[iorb][jorb][nsp] = 1.0;
                    }
                }
            }

        }
    }
    }

    for(nsp = 0; nsp < nspin; nsp++){
    for(iorb = 0; iorb < norb; iorb++){
        for(jorb = 0; jorb < norb; jorb++){
            // Permutation matrix is obtained by rounding off the square of overlap matrix in MO basis
            permut_mat[iorb][jorb][nsp] = round(pow(mo_overlap[iorb][jorb][nsp], 2) * tmp_sign[iorb][jorb][nsp]);
        }
    }
    }

    // Deallocate temporary array to save sign of overlap in MO basis
    for(iorb = 0; iorb < norb; iorb++){
    for(jorb = 0; jorb < norb; jorb++){
        free(tmp_sign[iorb][jorb]);
    }
    }
    for(jorb = 0; jorb < norb; jorb++){
        free(tmp_sign[jorb]);
    }

    free(tmp_sign);

}

// Routine to match phase of MO coefficients and orderings between two time steps
static void MO_phase_order(int nbasis, int norb, int nspin, double ***mo_coef_new, double ***permut_mat){

    double ***tmp_mo = malloc(norb * sizeof(double**));
    int ibasis, iorb, jorb, nsp;

    // Initialize temporary MO array; C'
    for(iorb = 0; iorb < norb; iorb++){
        tmp_mo[iorb] = malloc(nbasis * sizeof(double*));
        for(ibasis = 0; ibasis < nbasis; ibasis++){
        tmp_mo[iorb][ibasis] = malloc(nspin * sizeof(double));
		for(nsp = 0; nsp < nspin; nsp++){
            	tmp_mo[iorb][ibasis][nsp] = 0.0;
		}
        }
    }

    // Decide the phase and ordering for MO coefficients using permutation matrix; C' = O * C
    for(nsp = 0; nsp < nspin; nsp++){
    for(iorb = 0; iorb < norb; iorb++){
        for(ibasis = 0; ibasis < nbasis; ibasis++){

            for(jorb = 0; jorb < norb; jorb++){
                tmp_mo[iorb][ibasis][nsp] += permut_mat[iorb][jorb][nsp] * mo_coef_new[jorb][ibasis][nsp];
            }

        }
    }
    }

    // Apply new phase correction for the MO coefficients; C = C'
    for(nsp = 0; nsp < nspin; nsp++){
    for(iorb = 0; iorb < norb; iorb++){
        for(ibasis = 0; ibasis < nbasis; ibasis++){
            mo_coef_new[iorb][ibasis][nsp] = tmp_mo[iorb][ibasis][nsp];
        }
    }
    }

    // Deallocate temporary MO array; C'
    for(iorb = 0; iorb < norb; iorb++){
    for(ibasis = 0; ibasis < nbasis; ibasis++){
        free(tmp_mo[iorb][ibasis]);
    }
    }

    for(iorb = 0; iorb < norb; iorb++){
        free(tmp_mo[iorb]);
    }

    free(tmp_mo);

}

// Routine to match phase of CI coefficients and orderings between two time steps
// TODO : Is this correct method to match phase (or order) for CI coefficients?
static void CI_phase_order(int nst, int norb, int nspin, int *nocc, int *nvirt, int *orb_ini, int *orb_final,
    double ****ci_coef_old, double ****ci_coef_new, double ***permut_mat){

    double ***tmp_ci = malloc(norb * sizeof(double**));
    double ***tmp_ci_new = malloc(norb * sizeof(double**));
    int ist, iorb, jorb, aorb, borb, nsp;

    // Initialize temporary CI arrays; C' and C
    for(iorb = 0; iorb < norb; iorb++){
        tmp_ci[iorb] = malloc(norb * sizeof(double*));
        tmp_ci_new[iorb] = malloc(norb * sizeof(double*));
    	for(jorb = 0; jorb < norb; jorb++){
        tmp_ci[iorb][jorb] = malloc(nspin * sizeof(double));
        tmp_ci_new[iorb][jorb] = malloc(nspin * sizeof(double));
	}
    }

    // CI coefficients for S_0 are zero
    for(ist = 1; ist < nst; ist++){

	for(nsp = 0 ; nsp < nspin ; nsp++){

        for(iorb = orb_ini[nsp]; iorb < orb_final[nsp]; iorb++){
            for(aorb = orb_ini[nsp]; aorb < orb_final[nsp]; aorb++){
                // Assign CI coefficients at time t to new symmetric array
                if(iorb < nocc[nsp] && aorb >= nocc[nsp]){
                    tmp_ci[iorb][aorb][nsp] = ci_coef_new[ist][iorb][aorb - nocc[nsp]][nsp];
                }
                else if(iorb >= nocc[nsp] && aorb < nocc[nsp]){
                    tmp_ci[iorb][aorb][nsp] = ci_coef_new[ist][aorb][iorb - nocc[nsp]][nsp];
                }
                else{
                    tmp_ci[iorb][aorb][nsp] = 0.0;
                }
                // Initialize new empty array for phase correction
                tmp_ci_new[iorb][aorb][nsp] = 0.0;
            }
        }
	}

        // Decide the phase and ordering for CI coefficients using permutation matrix; C' = O * C * O
        // TODO : The phases for occupied and virtual orbitals are matched when permutation is diagonal matrix
	for(nsp = 0 ; nsp < nspin ; nsp++){
        for(jorb = orb_ini[nsp]; jorb < orb_final[nsp]; jorb++){
            for(borb = orb_ini[nsp]; borb < orb_final[nsp]; borb++){

                for(iorb = orb_ini[nsp]; iorb < orb_final[nsp]; iorb++){
                    for(aorb = orb_ini[nsp]; aorb < orb_final[nsp]; aorb++){
                        tmp_ci_new[jorb][borb][nsp] += permut_mat[jorb][iorb][nsp] * tmp_ci[iorb][aorb][nsp] * permut_mat[aorb][borb][nsp];
                    }
                }

            }
        }
	}

        // Apply new phase correction for the CI coefficients; C = C'
	for(nsp = 0 ; nsp < nspin ; nsp++){
        for(iorb = orb_ini[nsp]; iorb < nocc[nsp]; iorb++){
            for(aorb = 0; aorb < orb_final[nsp] - nocc[nsp]; aorb++){
                ci_coef_new[ist][iorb][aorb][nsp] = tmp_ci_new[iorb][nocc[nsp] + aorb][nsp];
            }
        }
	}

    }

    // Deallocate temporary CI arrays; C' and C
    for(iorb = 0; iorb < norb; iorb++){
    for(jorb = 0; jorb < norb; jorb++){
        free(tmp_ci[iorb][jorb]);
        free(tmp_ci_new[iorb][jorb]);
    }
    }
    for(jorb = 0; jorb < norb; jorb++){
        free(tmp_ci[jorb]);
        free(tmp_ci_new[jorb]);
    }

    free(tmp_ci);
    free(tmp_ci_new);

}

// Routine to match phase for the states between two time steps
static void state_phase(int nst, int nspin, int *nocc, int *nvirt, int *orb_ini, int *orb_final,
    double ****ci_coef_old, double ****ci_coef_new){

    double val;
    int ist, iorb, aorb, nsp;

    // CI coefficients for S_0 are zero
    for(ist = 1; ist < nst; ist++){

        val = 0.0;
	for(nsp = 0 ; nsp < nspin; nsp++){
        for(iorb = orb_ini[nsp]; iorb < nocc[nsp]; iorb++){
            for(aorb = 0; aorb < orb_final[nsp] - nocc[nsp]; aorb++){
                val += ci_coef_old[ist][iorb][aorb][nsp] * ci_coef_new[ist][iorb][aorb][nsp];
            }
        }

        if(val < 0.0){
            for(iorb = orb_ini[nsp]; iorb < nocc[nsp]; iorb++){
                for(aorb = 0; aorb < orb_final[nsp] - nocc[nsp]; aorb++){
                    ci_coef_new[ist][iorb][aorb][nsp] *= -1.0;
                }
            }
        }
	}

    }

}

// Routine to normalize CI coefficients
static void norm_CI_coef(int nst, int nspin, int *nocc, int *nvirt, int *orb_ini, int *orb_final, double ****ci_coef){

    double norm;
    int ist, iorb, aorb, nsp;

    // CI coefficients for S_0 are zero
    for(ist = 1; ist < nst; ist++){

        // Calculate normalization value for CI coefficients
        norm = 0.0;
	for(nsp = 0; nsp < nspin ; nsp++){
        for(iorb = orb_ini[nsp]; iorb < nocc[nsp]; iorb++){
            for(aorb = 0; aorb < orb_final[nsp] - nocc[nsp]; aorb++){
                norm += pow(ci_coef[ist][iorb][aorb][nsp], 2);
            }
        }
	}
        norm = sqrt(norm);

        // Normalize the CI coefficients
	for(nsp = 0; nsp < nspin ; nsp++){
        for(iorb = orb_ini[nsp]; iorb < nocc[nsp]; iorb++){
            for(aorb = 0; aorb < orb_final[nsp] - nocc[nsp]; aorb++){
                ci_coef[ist][iorb][aorb][nsp] /= norm;
            }
        }
	}

    }

}

