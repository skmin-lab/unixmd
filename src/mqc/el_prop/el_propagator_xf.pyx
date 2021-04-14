# cython: language_level=3
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.complex cimport complex
import numpy as np
cimport numpy as np

cdef extern from "rk4_xf.c":
    void rk4(int nat, int ndim, int nst, int nesteps, double dt, char *elec_object, \
        bint *l_coh, double *mass, double *energy, double *energy_old, double *sigma, \
        double **nacme, double **nacme_old, double **pos, double **qmom, double ***aux_pos, \
        double ***phase, double complex *coef, double complex **rho, int verbosity, \
        double *dotpopdec)

def el_run(md):
    cdef:
        char *elec_object_c
        bint *l_coh
        double *mass
        double *energy
        double *energy_old
        double *sigma
        double **nacme
        double **nacme_old
        double **pos
        double **qmom
        double ***aux_pos
        double ***phase
        double complex *coef
        double complex **rho
        double *dotpopdec

        bytes py_bytes
        int ist, jst, nst, nesteps, iat, aux_nat, aux_ndim, verbosity
        double dt

    # Assign size variables
    nst = md.mol.nst
    nesteps, dt = md.nesteps, md.dt
    aux_nat, aux_ndim = md.aux.nat, md.aux.ndim

    # Allocate variables
    l_coh = <bint*> PyMem_Malloc(nst * sizeof(bint))
    mass = <double*> PyMem_Malloc(aux_nat * sizeof(double))
    energy = <double*> PyMem_Malloc(nst * sizeof(double))
    energy_old = <double*> PyMem_Malloc(nst * sizeof(double))
    sigma = <double*> PyMem_Malloc(aux_nat * sizeof(double))

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
        pos[iat] = <double*> PyMem_Malloc(aux_ndim * sizeof(double))
        qmom[iat] = <double*> PyMem_Malloc(aux_ndim * sizeof(double))

    for ist in range(nst):
        aux_pos[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        phase[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        for iat in range(aux_nat):
            aux_pos[ist][iat] = <double*> PyMem_Malloc(aux_ndim * sizeof(double))
            phase[ist][iat] = <double*> PyMem_Malloc(aux_ndim * sizeof(double))
    
    # Debug related
    verbosity = md.verbosity
    if (verbosity >= 1):
        dotpopdec = <double*> PyMem_Malloc(nst * sizeof(double))

    # Assign variables from python to C
    for ist in range(nst):
        energy[ist] = md.mol.states[ist].energy
        energy_old[ist] = md.mol.states[ist].energy_old

    for iat in range(aux_nat):
        mass[iat] = md.aux.mass[iat]
        sigma[iat] = md.sigma[iat]

    for ist in range(nst):
        for jst in range(nst):
            nacme[ist][jst] = md.mol.nacme[ist, jst]
            nacme_old[ist][jst] = md.mol.nacme_old[ist, jst]

    for iat in range(aux_nat):
        for isp in range(aux_ndim):
            pos[iat][isp] = md.pos_0[iat, isp]

    for ist in range(nst):
        l_coh[ist] = md.l_coh[ist]
        for iat in range(aux_nat):
            for isp in range(aux_ndim):
                aux_pos[ist][iat][isp] = md.aux.pos[ist, iat, isp]
                phase[ist][iat][isp] = md.phase[ist, iat, isp]

    # Assign coef or rho with respect to propagation scheme
    if (md.elec_object == "coefficient"):

        coef = <double complex*> PyMem_Malloc(nst * sizeof(double complex))

        for ist in range(nst):
            coef[ist] = md.mol.states[ist].coef

    elif (md.elec_object == "density"):

        rho = <double complex**> PyMem_Malloc(nst * sizeof(double complex*))
        for ist in range(nst):
            rho[ist] = <double complex*> PyMem_Malloc(nst * sizeof(double complex))

        for ist in range(nst):
            for jst in range(nst):
                rho[ist][jst] = md.mol.rho[ist, jst]

    py_bytes = md.elec_object.encode()
    elec_object_c = py_bytes

    # Propagate electrons depending on the propagator
    if (md.propagator == "rk4"):
        rk4(aux_nat, aux_ndim, nst, nesteps, dt, elec_object_c, l_coh, mass, energy, \
            energy_old, sigma, nacme, nacme_old, pos, qmom, aux_pos, phase, coef, rho, verbosity, dotpopdec)

    # Assign variables from C to python
    if (md.elec_object == "coefficient"):

        for ist in range(nst):
            md.mol.states[ist].coef = coef[ist]

        for ist in range(nst):
            for jst in range(nst):
                md.mol.rho[ist, jst] = np.conj(md.mol.states[ist].coef) * md.mol.states[jst].coef

        PyMem_Free(coef)

    elif (md.elec_object == "density"):

        for ist in range(nst):
            for jst in range(nst):
                md.mol.rho[ist, jst] = rho[ist][jst]

        for ist in range(nst):
            PyMem_Free(rho[ist])
        PyMem_Free(rho)

    # Debug
    if (verbosity >= 1):
        for ist in range(nst):
            md.dotpopdec[ist] = dotpopdec[ist]
            md.dotpopnac[ist] = 0.
            for jst in range(nst):
                if (jst != ist):
                    md.dotpopnac[ist] -= 2. * nacme[ist][jst] * md.mol.rho.real[jst, ist]

    if (verbosity >= 2):
        for iat in range(aux_nat):
            for isp in range(aux_ndim):
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
    PyMem_Free(sigma)

    PyMem_Free(nacme)
    PyMem_Free(nacme_old)
    PyMem_Free(pos)
    PyMem_Free(qmom)

    PyMem_Free(aux_pos)
    PyMem_Free(phase)

    if (verbosity >= 1):
        PyMem_Free(dotpopdec)

