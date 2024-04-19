#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

// Routine to transform the coefficient (diabatic to adiabatic)
static void transform_d2a(int pst, double **unitary, double complex *coef_d, double complex *coef_a){

    double tmp_real, tmp_imag;
    int ist, jst;

    // C = U^T * D; U^-1 = U^T
    // Index for adiabatic states
    for(ist = 0; ist < pst; ist++){
        tmp_real = 0.0;
        tmp_imag = 0.0;
        // Index for diabatic states
        for(jst = 0; jst < pst; jst++){
            tmp_real += unitary[jst][ist] * creal(coef_d[jst]);
            tmp_imag += unitary[jst][ist] * cimag(coef_d[jst]);
        }
        coef_a[ist] = tmp_real + tmp_imag * I;
    }

}
#endif


