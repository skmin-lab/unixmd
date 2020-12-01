# cython: language_level=3
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.complex cimport complex
import numpy as np
cimport numpy as np

cdef extern from "rk4_xf.c":
    void rk4(int nat, int nsp, int nst, int nesteps, double dt, char *propagation, \
        bint *l_coh, double *mass, double *energy, double *energy_old, double *wsigma, \
        double **nacme, double **nacme_old, double **pos, double **qmom, double ***aux_pos, \
        double ***phase, double complex *coef, double complex **rho, int verbosity, \
        double *dotpopd)

def el_run(md):
    cdef:
        char *propagation_c
        bint *l_coh
        double *mass
        double *energy
        double *energy_old
        double *wsigma
        double **nacme
        double **nacme_old
        double **pos
        double **qmom
        double ***aux_pos
        double ***phase
        double complex *coef
        double complex **rho
        double *dotpopd

        bytes py_bytes
        int ist, jst, nst, nesteps, iat, aux_nat, aux_nsp, verbosity
        double dt

    # Assign size variables
    nst = md.mol.nst
    nesteps, dt = md.nesteps, md.dt
    aux_nat, aux_nsp = md.aux.nat, md.aux.nsp

    # Allocate variables
    l_coh = <bint*> PyMem_Malloc(nst * sizeof(bint))
    mass = <double*> PyMem_Malloc(aux_nat * sizeof(double))
    energy = <double*> PyMem_Malloc(nst * sizeof(double))
    energy_old = <double*> PyMem_Malloc(nst * sizeof(double))
    wsigma = <double*> PyMem_Malloc(aux_nat * sizeof(double))

    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))
    pos = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
    qmom = <double**> PyMem_Malloc(aux_nat * sizeof(double*))

    aux_pos = <double***> PyMem_Malloc(nst * sizeof(double**))
    phase = <double***> PyMem_Malloc(nst * sizeof(double**))

    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))

    for iat in range(aux_nat):
        pos[iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
        qmom[iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))

    for ist in range(nst):
        aux_pos[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        phase[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        for iat in range(aux_nat):
            aux_pos[ist][iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
            phase[ist][iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
    
    # Debug related
    verbosity = md.verbosity
    dotpopd = <double*> PyMem_Malloc(nst * sizeof(double))

    # Assign variables from python to C
    for ist in range(nst):
        energy[ist] = md.mol.states[ist].energy
        energy_old[ist] = md.mol.states[ist].energy_old

    for iat in range(aux_nat):
        mass[iat] = md.aux.mass[iat]
        wsigma[iat] = md.wsigma[iat]

    for ist in range(nst):
        for jst in range(nst):
            nacme[ist][jst] = md.mol.nacme[ist, jst]
            nacme_old[ist][jst] = md.mol.nacme_old[ist, jst]

    for iat in range(aux_nat):
        for isp in range(aux_nsp):
            pos[iat][isp] = md.pos_0[iat, isp]

    for ist in range(nst):
        l_coh[ist] = md.l_coh[ist]
        for iat in range(aux_nat):
            for isp in range(aux_nsp):
                aux_pos[ist][iat][isp] = md.aux.pos[ist, iat, isp]
                phase[ist][iat][isp] = md.phase[ist, iat, isp]

    # Assign coef or rho with respect to propagation scheme
    if (md.propagation == "coefficient"):

        coef = <double complex*> PyMem_Malloc(nst * sizeof(double complex))

        for ist in range(nst):
            coef[ist] = md.mol.states[ist].coef

    elif (md.propagation == "density"):

        rho = <double complex**> PyMem_Malloc(nst * sizeof(double complex*))
        for ist in range(nst):
            rho[ist] = <double complex*> PyMem_Malloc(nst * sizeof(double complex))

        for ist in range(nst):
            for jst in range(nst):
                rho[ist][jst] = md.mol.rho[ist, jst]

    py_bytes = md.propagation.encode()
    propagation_c = py_bytes

    # Propagate electrons depending on the solver
    if (md.solver == "rk4"):
        rk4(aux_nat, aux_nsp, nst, nesteps, dt, propagation_c, l_coh, mass, energy, \
            energy_old, wsigma, nacme, nacme_old, pos, qmom, aux_pos, phase, coef, rho, verbosity, dotpopd)

    # Assign variables from C to python
    if (md.propagation == "coefficient"):

        for ist in range(nst):
            md.mol.states[ist].coef = coef[ist]

        for ist in range(nst):
            for jst in range(nst):
                md.mol.rho[ist, jst] = np.conj(md.mol.states[ist].coef) * md.mol.states[jst].coef

        PyMem_Free(coef)

    elif (md.propagation == "density"):

        for ist in range(nst):
            for jst in range(nst):
                md.mol.rho[ist, jst] = rho[ist][jst]

        for ist in range(nst):
            PyMem_Free(rho[ist])
        PyMem_Free(rho)

    # Debug
    for ist in range(nst):
        md.dotpopd[ist] = dotpopd[ist]

    if (verbosity >= 2):
        for iat in range(aux_nat):
            for isp in range(aux_nsp):
                md.qmom[iat, isp] = qmom[iat][isp]

    # Deallocate variables
    for ist in range(nst):
        for iat in range(aux_nat):
            PyMem_Free(aux_pos[ist][iat])
            PyMem_Free(phase[ist][iat])

    for ist in range(nst):
        PyMem_Free(nacme[ist])
        PyMem_Free(nacme_old[ist])

    for iat in range(aux_nat):
        PyMem_Free(pos[iat])
        PyMem_Free(qmom[iat])

    for ist in range(nst):
        PyMem_Free(aux_pos[ist])
        PyMem_Free(phase[ist])

    PyMem_Free(l_coh)
    PyMem_Free(mass)
    PyMem_Free(energy)
    PyMem_Free(energy_old)
    PyMem_Free(wsigma)

    PyMem_Free(nacme)
    PyMem_Free(nacme_old)
    PyMem_Free(pos)
    PyMem_Free(qmom)

    PyMem_Free(aux_pos)
    PyMem_Free(phase)

    PyMem_Free(dotpopd)

