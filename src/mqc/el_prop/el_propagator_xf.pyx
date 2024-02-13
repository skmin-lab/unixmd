# cython: language_level=3
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.complex cimport complex
import numpy as np
cimport numpy as np

cdef extern from "rk4_xf.c":
    void rk4(int nst, int nesteps, double dt, char *elec_object, double *energy, double *energy_old, 
        double **nacme, double **nacme_old, double **k_lk, double complex *coef, double complex **rho, int verbosity, \
        double *dotpopdec)

def el_run(md):
    cdef:
        char *elec_object_c
        double *energy
        double *energy_old
        double **nacme
        double **nacme_old
        double **k_lk
        double complex *coef
        double complex **rho
        double *dotpopdec

        bytes py_bytes
        int ist, jst, nst, nesteps, verbosity
        double dt

    # Assign size variables
    nst = md.mol.nst
    nesteps, dt = md.nesteps, md.dt

    # Allocate variables
    energy = <double*> PyMem_Malloc(nst * sizeof(double))
    energy_old = <double*> PyMem_Malloc(nst * sizeof(double))

    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))
    k_lk = <double**> PyMem_Malloc(nst * sizeof(double*))

    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        k_lk[ist] = <double*> PyMem_Malloc(nst * sizeof(double))

    # Debug related
    verbosity = md.verbosity
    if (verbosity >= 1):
        dotpopdec = <double*> PyMem_Malloc(nst * sizeof(double))

    # Assign variables from python to C
    for ist in range(nst):
        energy[ist] = md.mol.states[ist].energy
        energy_old[ist] = md.mol.states[ist].energy_old

    for ist in range(nst):
        for jst in range(nst):
            nacme[ist][jst] = md.mol.nacme[ist, jst]
            nacme_old[ist][jst] = md.mol.nacme_old[ist, jst]
            k_lk[ist][jst] = md.k_lk[ist, jst]

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
        rk4(nst, nesteps, dt, elec_object_c, energy, energy_old, nacme, nacme_old, \
            k_lk, coef, rho, verbosity, dotpopdec)

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

    # Deallocate variables
    for ist in range(nst):
        PyMem_Free(nacme[ist])
        PyMem_Free(nacme_old[ist])
        PyMem_Free(k_lk[ist])

    PyMem_Free(energy)
    PyMem_Free(energy_old)

    PyMem_Free(nacme)
    PyMem_Free(nacme_old)
    PyMem_Free(k_lk)

    if (verbosity >= 1):
        PyMem_Free(dotpopdec)

