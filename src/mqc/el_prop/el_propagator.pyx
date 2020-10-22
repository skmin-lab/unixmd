# cython: language_level=3
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.complex cimport complex
import numpy as np
cimport numpy as np

cdef extern from "rk4.c":
    void rk4(int nst, int nesteps, double dt, char *propagation, double *energy, \
        double *energy_old, double **nacme, double **nacme_old, double complex *coef, \
        double complex **rho)

def el_run(md, molecule):
    cdef:
        char *propagation_c
        double *energy
        double *energy_old
        double **nacme
        double **nacme_old
        double complex *coef
        double complex **rho

        bytes py_bytes
        int ist, jst, nst, nesteps
        double dt

    # Assign size variables
    nst = molecule.nst
    nesteps, dt = md.nesteps, md.dt

    # Allocate variables
    energy = <double*> PyMem_Malloc(nst * sizeof(double))
    energy_old = <double*> PyMem_Malloc(nst * sizeof(double))

    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))

    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))

    # Assign variables from python to C
    for ist in range(nst):
        energy[ist] = molecule.states[ist].energy
        energy_old[ist] = molecule.states[ist].energy_old

    for ist in range(nst):
        for jst in range(nst):
            nacme[ist][jst] = molecule.nacme[ist, jst]
            nacme_old[ist][jst] = molecule.nacme_old[ist, jst]

    # Assign coef or rho with respect to propagation scheme
    if (md.propagation == "coefficient"):

        coef = <double complex*> PyMem_Malloc(nst * sizeof(double complex))

        for ist in range(nst):
            coef[ist] = molecule.states[ist].coef

    elif (md.propagation == "density"):

        rho = <double complex**> PyMem_Malloc(nst * sizeof(double complex*))
        for ist in range(nst):
            rho[ist] = <double complex*> PyMem_Malloc(nst * sizeof(double complex))

        for ist in range(nst):
            for jst in range(nst):
                rho[ist][jst] = molecule.rho[ist, jst]

    py_bytes = md.propagation.encode()
    propagation_c = py_bytes

    # Propagate electrons depending on the solver
    if (md.solver == "rk4"):
        rk4(nst, nesteps, dt, propagation_c, energy, energy_old, nacme, nacme_old, coef, rho)

    # Assign variables from C to python
    if (md.propagation == "coefficient"):

        for ist in range(nst):
            molecule.states[ist].coef = coef[ist]

        for ist in range(nst):
            for jst in range(nst):
                molecule.rho[ist, jst] = np.conj(molecule.states[ist].coef) * molecule.states[jst].coef

        PyMem_Free(coef)

    elif (md.propagation == "density"):

        for ist in range(nst):
            for jst in range(nst):
                molecule.rho[ist, jst] = rho[ist][jst]

        for ist in range(nst):
            PyMem_Free(rho[ist])
        PyMem_Free(rho)

    # Deallocate variables
    for ist in range(nst):
        PyMem_Free(nacme[ist])
        PyMem_Free(nacme_old[ist])

    PyMem_Free(energy)
    PyMem_Free(energy_old)

    PyMem_Free(nacme)
    PyMem_Free(nacme_old)


