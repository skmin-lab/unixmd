# cython: language_level=3
import numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.complex cimport complex
cimport numpy as np

cdef extern from "rk4.c":
    void RK4_coef(int nst, int nesteps, double dt, double complex *coef, double *energy, double *energy_old, double **nacme, double **nacme_old)
    void RK4_rho(int nst, int nesteps, double dt, double complex **rho, double *energy, double *energy_old, double **nacme, double **nacme_old)
    void RK4_coef_xf(int nst, int nesteps, double dt, double complex *coef, \
                  double *energy, double *energy_old, double **nacme, double **nacme_old, \
                  int nat, int nsp, bint *l_coh, double *wsigma, double *mass, double **pos,\
                  double ***aux_pos, double ***phase)
    void RK4_rho_xf(int nst, int nesteps, double dt, double complex **rho, 
                  double *energy, double *energy_old, double **nacme, double **nacme_old,
                  int nat, int nsp, double *mass, double **pos,
                  bint *l_coh, double *wsigma, double ***aux_pos, double ***phase)

def el_coef_xf(xfvars, molecule):
    cdef:
        double *energy
        double *energy_old
        double **nacme
        double **nacme_old
        double complex *coef
        double dt
        int ist, jst, nst, nesteps
        
        double *wsigma
        double *mass
        double **pos
        double ***aux_pos
        double ***phase
        int iat, aux_nat, aux_nsp
        bint *l_coh

    #
    nst = molecule.nst
    nesteps, dt = xfvars.nesteps, xfvars.dt
    energy = <double*> PyMem_Malloc(nst * sizeof(double))
    energy_old = <double*> PyMem_Malloc(nst * sizeof(double))
    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))
    coef = <double complex*> PyMem_Malloc(nst * sizeof(double complex))

    aux_nat, aux_nsp = xfvars.aux.nat, xfvars.aux.nsp
    wsigma = <double*> PyMem_Malloc(aux_nat * sizeof(double))
    mass = <double*> PyMem_Malloc(aux_nat * sizeof(double))
    pos = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
    aux_pos = <double***> PyMem_Malloc(nst * sizeof(double**))
    phase = <double***> PyMem_Malloc(nst * sizeof(double**))
    l_coh = <bint*> PyMem_Malloc(nst * sizeof(bint))
    
    #
    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        energy[ist] = molecule.states[ist].energy
        energy_old[ist] = molecule.states[ist].energy_old
        coef[ist] = molecule.states[ist].coef
        for jst in range(nst):
            nacme[ist][jst] = molecule.nacme[ist, jst]
            nacme_old[ist][jst] = molecule.nacme_old[ist, jst]

    for iat in range(aux_nat):
        wsigma[iat] = xfvars.wsigma[iat]
    
    for iat in range(aux_nat):
        pos[iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
        mass[iat] = xfvars.aux.mass[iat]
        for isp in range(aux_nsp):
            pos[iat][isp] = xfvars.pos_0[iat, isp]
    for ist in range(nst):
        aux_pos[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        phase[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        l_coh[ist] = xfvars.l_coh[ist]
        for iat in range(aux_nat):
            aux_pos[ist][iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
            phase[ist][iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
            for isp in range(aux_nsp):
                aux_pos[ist][iat][isp] = xfvars.aux.pos[ist, iat, isp]
                phase[ist][iat][isp] = xfvars.phase[ist, iat, isp]


    #
    #RK4_coef(nst, nesteps, dt, coef, energy, energy_old, nacme, nacme_old)
    RK4_coef_xf(nst, nesteps, dt, coef, energy, energy_old, nacme, nacme_old, \
        aux_nat, aux_nsp, l_coh, wsigma, mass, pos, aux_pos, phase)

    for ist in range(nst):
        molecule.states[ist].coef = coef[ist]

    for ist in range(nst):
        for jst in range(nst):
            molecule.rho[ist, jst] = np.conj(molecule.states[ist].coef) * molecule.states[jst].coef

    PyMem_Free(energy)
    PyMem_Free(energy_old)
    PyMem_Free(coef)
    for ist in range(nst):
        PyMem_Free(nacme[ist])
        PyMem_Free(nacme_old[ist])
    PyMem_Free(nacme)
    PyMem_Free(nacme_old)
    
    PyMem_Free(wsigma)

    PyMem_Free(mass)
    PyMem_Free(l_coh)
    for ist in range(nst):
        for iat in range(aux_nat):
            PyMem_Free(aux_pos[ist][iat])
            PyMem_Free(phase[ist][iat])
        PyMem_Free(aux_pos[ist])
        PyMem_Free(phase[ist])
    for iat in range(aux_nat):
        PyMem_Free(pos[iat])
    PyMem_Free(aux_pos)
    PyMem_Free(phase)
    PyMem_Free(pos)

def el_coef(nesteps_py, dt_py, molecule):
    cdef:
        double *energy
        double *energy_old
        double **nacme
        double **nacme_old
        double complex *coef
        int ist, jst, nst, nesteps
        double dt

    nst = molecule.nst
    nesteps, dt = nesteps_py, dt_py

    energy = <double*> PyMem_Malloc(nst * sizeof(double))
    energy_old = <double*> PyMem_Malloc(nst * sizeof(double))
    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))
    coef = <double complex*> PyMem_Malloc(nst * sizeof(double complex))

    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))

    for ist in range(nst):
        energy[ist] = molecule.states[ist].energy
        energy_old[ist] = molecule.states[ist].energy_old
        coef[ist] = molecule.states[ist].coef
        for jst in range(nst):
            nacme[ist][jst] = molecule.nacme[ist, jst]
            nacme_old[ist][jst] = molecule.nacme_old[ist, jst]

    RK4_coef(nst, nesteps, dt, coef, energy, energy_old, nacme, nacme_old)

    for ist in range(nst):
        molecule.states[ist].coef = coef[ist]

    for ist in range(nst):
        for jst in range(nst):
            molecule.rho[ist, jst] = np.conj(molecule.states[ist].coef) * molecule.states[jst].coef

    PyMem_Free(energy)
    PyMem_Free(energy_old)
    PyMem_Free(coef)
    for ist in range(nst):
        PyMem_Free(nacme[ist])
        PyMem_Free(nacme_old[ist])
    PyMem_Free(nacme)
    PyMem_Free(nacme_old)


def el_rho_xf(xfvars, molecule):
    cdef:
        double *energy
        double *energy_old
        double **nacme
        double **nacme_old
        double complex **rho
        double dt
        int ist, jst, nst, nesteps

        double *wsigma
        double *mass
        double **pos
        double ***aux_pos
        double ***phase
        int iat, aux_nat, aux_nsp
        bint *l_coh

    nst = molecule.nst
    nesteps, dt = xfvars.nesteps, xfvars.dt
    energy = <double*> PyMem_Malloc(nst * sizeof(double))
    energy_old = <double*> PyMem_Malloc(nst * sizeof(double))
    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))
    rho = <double complex**> PyMem_Malloc(nst * sizeof(double complex*))

    aux_nat, aux_nsp = xfvars.aux.nat, xfvars.aux.nsp
    wsigma = <double*> PyMem_Malloc(aux_nat * sizeof(double))
    mass = <double*> PyMem_Malloc(aux_nat * sizeof(double))
    pos = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
    aux_pos = <double***> PyMem_Malloc(nst * sizeof(double**))
    phase = <double***> PyMem_Malloc(nst * sizeof(double**))
    l_coh = <bint*> PyMem_Malloc(nst * sizeof(bint))


    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        energy[ist] = molecule.states[ist].energy
        energy_old[ist] = molecule.states[ist].energy_old
        rho[ist] = <double complex*> PyMem_Malloc(nst * sizeof(double complex))
        for jst in range(nst):
            nacme[ist][jst] = molecule.nacme[ist, jst]
            nacme_old[ist][jst] = molecule.nacme_old[ist, jst]
            rho[ist][jst] = molecule.rho[ist, jst]


    for iat in range(aux_nat):
        wsigma[iat] = xfvars.wsigma[iat]

    
    # xf related variables allocataion
    for iat in range(aux_nat):
        pos[iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
        mass[iat] = xfvars.aux.mass[iat]
        for isp in range(aux_nsp):
            pos[iat][isp] = xfvars.pos_0[iat, isp]
    for ist in range(nst):
        aux_pos[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        phase[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        l_coh[ist] = xfvars.l_coh[ist]
        for iat in range(aux_nat):
            aux_pos[ist][iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
            phase[ist][iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
            for isp in range(aux_nsp):
                aux_pos[ist][iat][isp] = xfvars.aux.pos[ist, iat, isp]
                phase[ist][iat][isp] = xfvars.phase[ist, iat, isp]


    RK4_rho_xf(nst, nesteps, dt, rho, energy, energy_old, nacme, nacme_old, \
        aux_nat, aux_nsp, mass, pos, l_coh, wsigma, aux_pos, phase)

    for ist in range(nst):
        for jst in range(nst):
            molecule.rho[ist, jst] = rho[ist][jst]

    PyMem_Free(energy)
    PyMem_Free(energy_old)
    for ist in range(nst):
        PyMem_Free(nacme[ist])
        PyMem_Free(nacme_old[ist])
        PyMem_Free(rho[ist])
    PyMem_Free(nacme)
    PyMem_Free(nacme_old)
    PyMem_Free(rho)

    PyMem_Free(wsigma)
    
    PyMem_Free(mass)
    PyMem_Free(l_coh)
    for ist in range(nst):
        for iat in range(aux_nat):
            PyMem_Free(aux_pos[ist][iat])
            PyMem_Free(phase[ist][iat])
        PyMem_Free(aux_pos[ist])
        PyMem_Free(phase[ist])
    for iat in range(aux_nat):
        PyMem_Free(pos[iat])
    PyMem_Free(aux_pos)
    PyMem_Free(phase)
    PyMem_Free(pos)

def el_rho(nesteps_py, dt_py, molecule):
    cdef:
        double *energy
        double *energy_old
        double **nacme
        double **nacme_old
        double complex **rho
        int ist, jst, nst, nesteps
        double dt

    nst = molecule.nst
    nesteps, dt = nesteps_py, dt_py

    energy = <double*> PyMem_Malloc(nst * sizeof(double))
    energy_old = <double*> PyMem_Malloc(nst * sizeof(double))
    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))
    rho = <double complex**> PyMem_Malloc(nst * sizeof(double complex*))

    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        rho[ist] = <double complex*> PyMem_Malloc(nst * sizeof(double complex))

    for ist in range(nst):
        energy[ist] = molecule.states[ist].energy
        energy_old[ist] = molecule.states[ist].energy_old
        for jst in range(nst):
            nacme[ist][jst] = molecule.nacme[ist, jst]
            nacme_old[ist][jst] = molecule.nacme_old[ist, jst]
            rho[ist][jst] = molecule.rho[ist, jst]

    RK4_rho(nst, nesteps, dt, rho, energy, energy_old, nacme, nacme_old)

    for ist in range(nst):
        for jst in range(nst):
            molecule.rho[ist, jst] = rho[ist][jst]

    PyMem_Free(energy)
    PyMem_Free(energy_old)
    for ist in range(nst):
        PyMem_Free(nacme[ist])
        PyMem_Free(nacme_old[ist])
        PyMem_Free(rho[ist])
    PyMem_Free(nacme)
    PyMem_Free(nacme_old)
    PyMem_Free(rho)


