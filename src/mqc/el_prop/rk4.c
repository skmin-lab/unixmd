#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

//double dot(double complex* u, double complex* v, int nst);
//void cdot(int nst, double complex* c, double* e, double **dv, double complex *c_dot);

static double dot(double complex *u, double complex *v, int nst){
    double complex sum = 0.0 + 0.0 * I;
    double norm;
    int ist;
    for(ist = 0; ist < nst; ist++){
        sum += conj(u[ist]) * v[ist];
    }
    norm = creal(sum);
    return norm;
}

static void cdot(int nst, double complex *c, double *e, double **dv, double complex *c_dot){

    double complex *na_term = malloc(nst * sizeof(double complex));
    int ist, jst;
    double egs;

    for(ist = 0; ist < nst; ist++){
        na_term[ist] = 0.0 + 0.0 * I;
        for(jst = 0; jst < nst; jst++){
//            if(ist != jst){
                na_term[ist] -= dv[ist][jst] * c[jst];
//            }
        }
    }
    egs = e[0];
    for(ist = 0; ist < nst; ist++){
        c_dot[ist] = -1.0 * I * c[ist] * (e[ist] - egs) + na_term[ist];
    }

    free(na_term);
}

static void rhodot(int nst, double complex **rho, double *e, double **dv, double complex **rho_dot){
    
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
}
static void xf_cdot(int nst, int nat, int nsp, int *l_coh, double wsigma, 
                        double *mass, double complex *c, double **pos, double ***aux_pos, double ***phase, 
                        double complex *xfcdot){
    
    int ist, jst, iat, isp;
    double **qmom = malloc(nat * sizeof(double*));
    double **dec = malloc(nst * sizeof(double*));
    double *rho = malloc(nst * sizeof(double)); 
    
    for(iat = 0; iat < nat; iat++){
        qmom[iat] = malloc(nsp * sizeof(double));
        for(isp = 0; isp < nsp; isp++){
            qmom[iat][isp] = 0.0;
        }
    }
    for(ist = 0; ist < nst; ist++){
        dec[ist] = malloc(nst * sizeof(double));
        rho[ist] = creal(conj(c[ist]) * c[ist]);
    }
    
    // get qmom
    for(ist = 0; ist < nst; ist++){
        if(l_coh[ist] == 1){
            for(iat = 0; iat < nat; iat++){
                for(isp = 0; isp < nsp; isp++){
                    qmom[iat][isp] = qmom[iat][isp] + rho[ist] * 
                    0.5 * (pos[iat][isp] - aux_pos[ist][iat][isp]) / pow(wsigma,2.0) / mass[iat];
                }
            }
        }
    }
    for(ist = 0; ist < nst; ist++){
        dec[ist][ist] = 0.0;
        for(jst = ist + 1; jst < nst; jst++){
            dec[ist][jst] = 0.0;
            if(l_coh[ist] == 1 && l_coh[jst] == 1){
                for(iat = 0; iat < nat; iat++){
                    for(isp = 0; isp < nsp; isp++){
                        dec[ist][jst] = dec[ist][jst] + qmom[iat][isp]*(phase[ist][iat][isp]-phase[jst][iat][isp]);
                    }
                }
            }
            dec[jst][ist] = -1.0 * dec[ist][jst];
        }
    }
    for(ist = 0; ist < nst; ist++){
        xfcdot[ist] = 0.0 + 0.0 * I;
        for(jst = 0; jst < nst; jst++){
            xfcdot[ist] = xfcdot[ist] - rho[jst] * dec[jst][ist] * c[ist];
        }
    }
    free(rho);
    for(ist = 0; ist < nst; ist++){
        free(dec[ist]);
    }
    for(iat = 0; iat < nat; iat++){
        free(qmom[iat]);
    }
    free(dec);
    free(qmom);
}

static void xf_rhodot(int nst, int nat, int nsp, int *l_coh, double wsigma, 
                      double *mass, double complex **rho, double **pos, double ***aux_pos, double ***phase, 
                      double complex **xfrhodot){
    int ist, jst, kst, iat, isp;
    double **qmom = malloc(nat * sizeof(double*));
    double **dec = malloc(nst * sizeof(double*));
    
    for(iat = 0; iat < nat; iat++){
        qmom[iat] = malloc(nsp * sizeof(double));
        for(isp = 0; isp < nsp; isp++){
            qmom[iat][isp] = 0.0;
        }
    }
    for(ist = 0; ist < nst; ist++){
        dec[ist] = malloc(nst * sizeof(double));
    }
    
    // get qmom
    for(ist = 0; ist < nst; ist++){
        if(l_coh[ist] == 1){
            for(iat = 0; iat < nat; iat++){
                for(isp = 0; isp < nsp; isp++){
                    qmom[iat][isp] = qmom[iat][isp] + creal(rho[ist][ist]) * 
                    0.5 * (pos[iat][isp] - aux_pos[ist][iat][isp]) / pow(wsigma,2.0) / mass[iat];
                }
            }
        }
    }
    for(ist = 0; ist < nst; ist++){
        dec[ist][ist] = 0.0;
        for(jst = ist + 1; jst < nst; jst++){
            dec[ist][jst] = 0.0;
            if(l_coh[ist] == 1 && l_coh[jst] == 1){
                for(iat = 0; iat < nat; iat++){
                    for(isp = 0; isp < nsp; isp++){
                        dec[ist][jst] = dec[ist][jst] + qmom[iat][isp]*(phase[ist][iat][isp]-phase[jst][iat][isp]);
                    }
                }
            }
            dec[jst][ist] = -1.0 * dec[ist][jst];
        }
    }
    for(ist = 0; ist < nst; ist++){
        xfrhodot[ist][ist] = 0.0 + 0.0 * I; 
        for(kst = 0; kst < nst; kst++){
            xfrhodot[ist][ist] -= 2.0 * dec[kst][ist] * rho[ist][kst] * rho[kst][ist];
        }
        for(jst = ist + 1; jst < nst; jst++){
            xfrhodot[ist][jst] = 0.0 + 0.0 * I;
            for(kst = 0; kst < nst; kst++){
                xfrhodot[ist][jst] -= (dec[kst][ist]+dec[kst][jst]) * rho[ist][kst] * rho[kst][jst];
            }
            xfrhodot[jst][ist] = conj(xfrhodot[ist][jst]);
        }
    }
    for(ist = 0; ist < nst; ist++){
        free(dec[ist]);
    }
    for(iat = 0; iat < nat; iat++){
        free(qmom[iat]);
    }
    free(dec);
    free(qmom);
}

static void RK4_coef(int nst, int nesteps, double dt, double complex *coef, double *energy, double *energy_old, double **nacme, double **nacme_old){
    double complex *c_dot = malloc(nst * sizeof(double complex));
    double complex *coef_new = malloc(nst * sizeof(double complex));
    double complex *k1 = malloc(nst * sizeof(double complex));
    double complex *k2 = malloc(nst * sizeof(double complex));
    double complex *k3 = malloc(nst * sizeof(double complex));
    double complex *k4 = malloc(nst * sizeof(double complex));
    double complex *kfunction = malloc(nst * sizeof(double complex));
    double complex *variation = malloc(nst * sizeof(double complex));
    double *eenergy = malloc(nst * sizeof(double));
    double **dv = malloc(nst * sizeof(double*));
    double *na_term = malloc(nst * sizeof(double));
    int iestep, ist, jst;
    double frac, edt, norm;
    frac = 1.0 / (double)nesteps;
    edt = dt * frac;

    for(ist = 0; ist < nst; ist++){
        dv[ist] = malloc(nst * sizeof(double));
    }
    
    for(iestep = 0; iestep < nesteps; iestep++){
        
        for(ist = 0; ist < nst; ist++){
            eenergy[ist] = energy_old[ist] + (energy[ist] - energy_old[ist]) *(double) iestep * frac;
            for(jst = 0; jst < nst; jst++){
                dv[ist][jst] = nacme_old[ist][jst] + (nacme[ist][jst] - nacme_old[ist][jst]) *(double) iestep * frac;
            }
        }

        cdot(nst, coef, eenergy, dv, c_dot);

        for(ist = 0; ist < nst; ist++){
            k1[ist] = edt * c_dot[ist];
            kfunction[ist] = 0.5 * k1[ist];
            coef_new[ist] = coef[ist] + kfunction[ist];
        }

        cdot(nst, coef_new, eenergy, dv, c_dot);

        for(ist = 0; ist < nst; ist++){
            k2[ist] = edt * c_dot[ist];
            kfunction[ist] = 0.5 * (-1.0 + sqrt(2.0)) *  k1[ist] + (1.0 - 0.5 * sqrt(2.0)) * k2[ist];
            coef_new[ist] = coef[ist] + kfunction[ist];
        }

        cdot(nst, coef_new, eenergy, dv, c_dot);

        for(ist = 0; ist < nst; ist++){
            k3[ist] = edt * c_dot[ist];
            kfunction[ist] = - 0.5 * sqrt(2.0) * k2[ist] + (1.0 + 0.5 * sqrt(2.0)) * k3[ist];
            coef_new[ist] = coef[ist] + kfunction[ist];
        }

        cdot(nst, coef_new, eenergy, dv, c_dot);

        for(ist = 0; ist < nst; ist++){
            k4[ist] = edt * c_dot[ist];
            variation[ist] = (k1[ist] + (2.0 - sqrt(2.0)) * k2[ist] + (2.0 + sqrt(2.0)) * k3[ist] + k4[ist]) / 6.0;
            coef_new[ist] = coef[ist] + variation[ist];
        }
        norm = dot(coef_new, coef_new, nst);
        
        for(ist = 0; ist < nst; ist++){
            coef_new[ist] /= sqrt(norm);
            coef[ist] = coef_new[ist];
        }
    }
    
    for(ist = 0; ist < nst; ist++){
        na_term[ist] = 0.0;
        for(jst = 0; jst < nst; jst++){
            if(jst != ist){
                na_term[ist] -= dv[ist][jst] * coef[jst];
            }
        }
        //printf("RHODOT_NAC %d %f\n",ist+1,2.0*creal(conj(coef[ist])) * na_term[ist]);
    }
    /*printf("RK4_COEF : NORM = %15.8f\n", creal(norm));*/

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(kfunction);
    free(variation);
    free(c_dot);
    free(coef_new);
    free(eenergy);
    free(na_term);
    for(ist = 0; ist < nst; ist++){
        free(dv[ist]);
    }
    free(dv);
}

static void RK4_rho(int nst, int nesteps, double dt, double complex **rho, double *energy, double *energy_old, double **nacme, double **nacme_old){
    double complex **rho_dot = malloc(nst * sizeof(double complex *));
    double complex **rho_new = malloc(nst * sizeof(double complex *));
    double complex **k1 = malloc(nst * sizeof(double complex *));
    double complex **k2 = malloc(nst * sizeof(double complex *));
    double complex **k3 = malloc(nst * sizeof(double complex *));
    double complex **k4 = malloc(nst * sizeof(double complex *));
    double complex **kfunction = malloc(nst * sizeof(double complex *));
    double complex **variation = malloc(nst * sizeof(double complex *));
    double *eenergy = malloc(nst * sizeof(double));
    double **dv = malloc(nst * sizeof(double*));
    double *na_term = malloc(nst * sizeof(double));
    int iestep, ist, jst;
    double frac, edt, norm;
    frac = 1.0 / (double)nesteps;
    edt = dt * frac;

    for(ist = 0; ist < nst; ist++){
        dv[ist] = malloc(nst * sizeof(double));
        rho_dot[ist] = malloc(nst * sizeof(double complex));
        rho_new[ist] = malloc(nst * sizeof(double complex));
        k1[ist] = malloc(nst * sizeof(double complex));
        k2[ist] = malloc(nst * sizeof(double complex));
        k3[ist] = malloc(nst * sizeof(double complex));
        k4[ist] = malloc(nst * sizeof(double complex));
        kfunction[ist] = malloc(nst * sizeof(double complex));
        variation[ist] = malloc(nst * sizeof(double complex));
    }
    
    for(iestep = 0; iestep < nesteps; iestep++){
        
        for(ist = 0; ist < nst; ist++){
            eenergy[ist] = energy_old[ist] + (energy[ist] - energy_old[ist]) *(double) iestep * frac;
            for(jst = 0; jst < nst; jst++){
                dv[ist][jst] = nacme_old[ist][jst] + (nacme[ist][jst] - nacme_old[ist][jst]) *(double) iestep * frac;
            }
        }

        rhodot(nst, rho, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k1[ist][jst] = edt * rho_dot[ist][jst];
                kfunction[ist][jst] = 0.5 * k1[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        rhodot(nst, rho_new, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k2[ist][jst] = edt * rho_dot[ist][jst];
                kfunction[ist][jst] = 0.5 * (-1.0 + sqrt(2.0)) *  k1[ist][jst] + (1.0 - 0.5 * sqrt(2.0)) * k2[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        rhodot(nst, rho_new, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k3[ist][jst] = edt * rho_dot[ist][jst];
                kfunction[ist][jst] = - 0.5 * sqrt(2.0) * k2[ist][jst] + (1.0 + 0.5 * sqrt(2.0)) * k3[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        rhodot(nst, rho_new, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k4[ist][jst] = edt * rho_dot[ist][jst];
                variation[ist][jst] = (k1[ist][jst] + (2.0 - sqrt(2.0)) * k2[ist][jst] + (2.0 + sqrt(2.0)) * k3[ist][jst] + k4[ist][jst]) / 6.0;
                rho[ist][jst] += variation[ist][jst];
            }
        }
    }
    
    for(ist = 0; ist < nst; ist++){
        na_term[ist] = 0.0;
        for(jst = 0; jst < nst; jst++){
            if(jst != ist){
                na_term[ist] -= 2.0 * dv[ist][jst] * creal(rho[ist][jst]);
            }
        }
        /*printf("RHODOT_NAC %d %f\n",ist+1, na_term[ist]);*/
    }
    norm = 0.0;
    for(ist = 0; ist < nst; ist++){
        norm += creal(rho[ist][ist]);
    }
    /*printf("RK4_COEF : NORM = %15.8f\n", norm);*/

    free(eenergy);
    free(na_term);
    for(ist = 0; ist < nst; ist++){
        free(dv[ist]);
        free(k1[ist]);
        free(k2[ist]);
        free(k3[ist]);
        free(k4[ist]);
        free(kfunction[ist]);
        free(variation[ist]);
        free(rho_dot[ist]);
        free(rho_new[ist]);
    }
    free(dv);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(kfunction);
    free(variation);
    free(rho_dot);
    free(rho_new);
}
static void RK4_coef_xf(int nst, int nesteps, double dt, double complex *coef, 
                     double *energy, double *energy_old, double **nacme, double **nacme_old,  // normal rk4_coef up to this arg
                     int nat, int nsp, int *l_coh, double wsigma, double *mass, double **pos, double ***aux_pos, double ***phase){
    double complex *c_dot = malloc(nst * sizeof(double complex));
    double complex *coef_new = malloc(nst * sizeof(double complex));
    double complex *k1 = malloc(nst * sizeof(double complex));
    double complex *k2 = malloc(nst * sizeof(double complex));
    double complex *k3 = malloc(nst * sizeof(double complex));
    double complex *k4 = malloc(nst * sizeof(double complex));
    double complex *kfunction = malloc(nst * sizeof(double complex));
    double complex *variation = malloc(nst * sizeof(double complex));
    double *eenergy = malloc(nst * sizeof(double));
    double **dv = malloc(nst * sizeof(double*));
    double *na_term = malloc(nst * sizeof(double));
    int iestep, ist, jst;
    double frac, edt, norm;
    
    double complex *xf_c_dot = malloc(nst * sizeof(double complex));
    
    frac = 1.0 / (double)nesteps;
    edt = dt * frac;

    for(ist = 0; ist < nst; ist++){
        dv[ist] = malloc(nst * sizeof(double));
    }
    
    for(iestep = 0; iestep < nesteps; iestep++){
        
        xf_cdot(nst, nat, nsp, l_coh, wsigma, mass, coef, pos, aux_pos, phase, xf_c_dot);
        for(ist = 0; ist < nst; ist++){
            eenergy[ist] = energy_old[ist] + (energy[ist] - energy_old[ist]) *(double) iestep * frac;
            for(jst = 0; jst < nst; jst++){
                dv[ist][jst] = nacme_old[ist][jst] + (nacme[ist][jst] - nacme_old[ist][jst]) *(double) iestep * frac;
            }
        }

        cdot(nst, coef, eenergy, dv, c_dot);
        //xf_cdot(nst, nat, nsp, l_coh, wsigma, mass, coef, pos, aux_pos, phase, xf_c_dot);
        
        for(ist = 0; ist < nst; ist++){
            //k1[ist] = edt * c_dot[ist];
            k1[ist] = edt * (c_dot[ist] + xf_c_dot[ist]);
            kfunction[ist] = 0.5 * k1[ist];
            coef_new[ist] = coef[ist] + kfunction[ist];
        }

        cdot(nst, coef_new, eenergy, dv, c_dot);
        //xf_cdot(nst, nat, nsp, l_coh, wsigma, mass, coef_new, pos, aux_pos, phase, xf_c_dot);

        for(ist = 0; ist < nst; ist++){
            //k2[ist] = edt * c_dot[ist];
            k2[ist] = edt * (c_dot[ist] + xf_c_dot[ist]);
            kfunction[ist] = 0.5 * (-1.0 + sqrt(2.0)) *  k1[ist] + (1.0 - 0.5 * sqrt(2.0)) * k2[ist];
            coef_new[ist] = coef[ist] + kfunction[ist];
        }

        cdot(nst, coef_new, eenergy, dv, c_dot);
        //xf_cdot(nst, nat, nsp, l_coh, wsigma, mass, coef_new, pos, aux_pos, phase, xf_c_dot);

        for(ist = 0; ist < nst; ist++){
            //k3[ist] = edt * c_dot[ist];
            k3[ist] = edt * (c_dot[ist] + xf_c_dot[ist]);
            kfunction[ist] = - 0.5 * sqrt(2.0) * k2[ist] + (1.0 + 0.5 * sqrt(2.0)) * k3[ist];
            coef_new[ist] = coef[ist] + kfunction[ist];
        }

        cdot(nst, coef_new, eenergy, dv, c_dot);
        //xf_cdot(nst, nat, nsp, l_coh, wsigma, mass, coef_new, pos, aux_pos, phase, xf_c_dot);

        for(ist = 0; ist < nst; ist++){
            //k4[ist] = edt * c_dot[ist];
            k4[ist] = edt * (c_dot[ist] + xf_c_dot[ist]);
            variation[ist] = (k1[ist] + (2.0 - sqrt(2.0)) * k2[ist] + (2.0 + sqrt(2.0)) * k3[ist] + k4[ist]) / 6.0;
            coef_new[ist] = coef[ist] + variation[ist];
        }
        norm = dot(coef_new, coef_new, nst);
        
        for(ist = 0; ist < nst; ist++){
            coef_new[ist] /= sqrt(norm);
            coef[ist] = coef_new[ist];
        }
    }
    
    for(ist = 0; ist < nst; ist++){
        //printf(" TSHXF RHODOT           %d %22.15E\n",ist+1,2.0*creal(xf_c_dot[ist] * conj(coef[ist])));
        na_term[ist] = 0.0;
        c_dot[ist] = variation[ist]/edt;
        for(jst = 0; jst < nst; jst++){
            if(jst != ist){
                na_term[ist] -= dv[ist][jst] * coef[jst];
            }
        }
        //printf("RHODOT_NAC %d %f\n",ist+1,2.0*creal(conj(coef[ist])) * na_term[ist]);
        //printf("RHODOT_TOT %d %f\n",ist+1,2.0*creal(conj(coef[ist])) * c_dot[ist]);
    }
    /*printf("RK4_COEF : NORM = %15.8f\n", creal(norm));*/

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(kfunction);
    free(variation);
    free(c_dot);
    free(coef_new);
    free(eenergy);
    free(na_term);
    for(ist = 0; ist < nst; ist++){
        free(dv[ist]);
    }
    free(dv);
    
    free(xf_c_dot);
}

static void RK4_rho_xf(int nst, int nesteps, double dt, double complex **rho, 
                       double *energy, double *energy_old, double **nacme, double **nacme_old,
                       int nat, int nsp, double *mass, double **pos,
                       int *l_coh, double wsigma, double ***aux_pos, double ***phase){
    double complex **rho_dot = malloc(nst * sizeof(double complex *));
    double complex **rho_new = malloc(nst * sizeof(double complex *));
    double complex **k1 = malloc(nst * sizeof(double complex *));
    double complex **k2 = malloc(nst * sizeof(double complex *));
    double complex **k3 = malloc(nst * sizeof(double complex *));
    double complex **k4 = malloc(nst * sizeof(double complex *));
    double complex **kfunction = malloc(nst * sizeof(double complex *));
    double complex **variation = malloc(nst * sizeof(double complex *));
    double *eenergy = malloc(nst * sizeof(double));
    double **dv = malloc(nst * sizeof(double*));
    double *na_term = malloc(nst * sizeof(double));
    int iestep, ist, jst;
    double frac, edt, norm;

    double complex **xf_rho_dot = malloc(nst * sizeof(double complex *));

    frac = 1.0 / (double)nesteps;
    edt = dt * frac;

    for(ist = 0; ist < nst; ist++){
        dv[ist] = malloc(nst * sizeof(double));
        rho_dot[ist] = malloc(nst * sizeof(double complex));
        rho_new[ist] = malloc(nst * sizeof(double complex));
        k1[ist] = malloc(nst * sizeof(double complex));
        k2[ist] = malloc(nst * sizeof(double complex));
        k3[ist] = malloc(nst * sizeof(double complex));
        k4[ist] = malloc(nst * sizeof(double complex));
        kfunction[ist] = malloc(nst * sizeof(double complex));
        variation[ist] = malloc(nst * sizeof(double complex));
        
        xf_rho_dot[ist] = malloc(nst * sizeof(double complex));
    }
     
    for(ist = 0; ist < nst; ist++){
        for(jst = 0; jst < nst; jst++){
            xf_rho_dot[ist][jst] = 0.0 + 0.0 * I;
        }
    }
    for(iestep = 0; iestep < nesteps; iestep++){
        
        xf_rhodot(nst, nat, nsp, l_coh, wsigma, mass, rho, pos, aux_pos, phase, xf_rho_dot);
        for(ist = 0; ist < nst; ist++){
            eenergy[ist] = energy_old[ist] + (energy[ist] - energy_old[ist]) *(double) iestep * frac;
            for(jst = 0; jst < nst; jst++){
                dv[ist][jst] = nacme_old[ist][jst] + (nacme[ist][jst] - nacme_old[ist][jst]) *(double) iestep * frac;
            }
        }

        rhodot(nst, rho, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k1[ist][jst] = edt * (rho_dot[ist][jst] + xf_rho_dot[ist][jst]);
                kfunction[ist][jst] = 0.5 * k1[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        rhodot(nst, rho_new, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k2[ist][jst] = edt * (rho_dot[ist][jst] + xf_rho_dot[ist][jst]);
                kfunction[ist][jst] = 0.5 * (-1.0 + sqrt(2.0)) *  k1[ist][jst] + (1.0 - 0.5 * sqrt(2.0)) * k2[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        rhodot(nst, rho_new, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k3[ist][jst] = edt * (rho_dot[ist][jst] + xf_rho_dot[ist][jst]);
                kfunction[ist][jst] = - 0.5 * sqrt(2.0) * k2[ist][jst] + (1.0 + 0.5 * sqrt(2.0)) * k3[ist][jst];
                rho_new[ist][jst] = rho[ist][jst] + kfunction[ist][jst];
            }
        }

        rhodot(nst, rho_new, eenergy, dv, rho_dot);

        for(ist = 0; ist < nst; ist++){
            for(jst = 0; jst < nst; jst++){
                k4[ist][jst] = edt * (rho_dot[ist][jst] + xf_rho_dot[ist][jst]);
                variation[ist][jst] = (k1[ist][jst] + (2.0 - sqrt(2.0)) * k2[ist][jst] + (2.0 + sqrt(2.0)) * k3[ist][jst] + k4[ist][jst]) / 6.0;
                rho[ist][jst] += variation[ist][jst];
            }
        }
    }
    
    for(ist = 0; ist < nst; ist++){
        na_term[ist] = 0.0;
        for(jst = 0; jst < nst; jst++){
            if(jst != ist){
                na_term[ist] -= 2.0 * dv[ist][jst] * creal(rho[ist][jst]);
            }
        }
        //printf("RHODOT_NAC %d %f\n",ist+1, na_term[ist]);
    }
    norm = 0.0;
    for(ist = 0; ist < nst; ist++){
        norm += creal(rho[ist][ist]);
    }
    //printf("RK4_COEF : NORM = %15.8f\n", norm);

    free(eenergy);
    free(na_term);
    for(ist = 0; ist < nst; ist++){
        free(dv[ist]);
        free(k1[ist]);
        free(k2[ist]);
        free(k3[ist]);
        free(k4[ist]);
        free(kfunction[ist]);
        free(variation[ist]);
        free(rho_dot[ist]);
        free(rho_new[ist]);

        free(xf_rho_dot[ist]);
    }
    free(dv);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(kfunction);
    free(variation);
    free(rho_dot);

    free(xf_rho_dot);
    free(rho_new);
}
