# cython: language_level=3
import numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.complex cimport complex
cimport numpy as np

cdef extern from "rk4.c":
    void rk4(char *propagation, int nst, int nesteps, double dt, double complex *coef, double complex **rho, double *energy, double *energy_old, double **nacme, double **nacme_old)

def el_run(md, molecule):
    cdef:
        double *energy
        double *energy_old
        double **nacme
        double **nacme_old
        double complex *coef
        double complex **rho
        int ist, jst, nst, nesteps
        bytes py_bytes
        char *propagation_c
        double dt

    nst = molecule.nst
    nesteps, dt = md.nesteps, md.dt

    energy = <double*> PyMem_Malloc(nst * sizeof(double))
    energy_old = <double*> PyMem_Malloc(nst * sizeof(double))
    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))

    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))

    for ist in range(nst):
        energy[ist] = molecule.states[ist].energy
        energy_old[ist] = molecule.states[ist].energy_old
        for jst in range(nst):
            nacme[ist][jst] = molecule.nacme[ist, jst]
            nacme_old[ist][jst] = molecule.nacme_old[ist, jst]

    if md.propagation == "coefficient":
        coef = <double complex*> PyMem_Malloc(nst * sizeof(double complex))
        for ist in range(nst):
            coef[ist] = molecule.states[ist].coef
    elif md.propagation == "density":
        rho = <double complex**> PyMem_Malloc(nst * sizeof(double complex*))
        for ist in range(nst):
            rho[ist] = <double complex*> PyMem_Malloc(nst * sizeof(double complex))
            for jst in range(nst):
                rho[ist][jst] = molecule.rho[ist, jst]
    #else:
    #    raise ValueError (f"( {md.md_type}.{call_name()} ) Other propagator not implemented! {self.propagation}")

    py_bytes = md.propagation.encode() 
    propagation_c = py_bytes
    if md.solver == "rk4":
        rk4(propagation_c, nst, nesteps, dt, coef, rho, energy, energy_old, nacme, nacme_old)
    #else:  
    #    raise ValueError (f"( {md.md_type}.{call_name()} ) Other solver not implemented! {self.propagation}")
      
    if md.propagation == "coefficient":
        for ist in range(nst):
            molecule.states[ist].coef = coef[ist]

        for ist in range(nst):
            for jst in range(nst):
                molecule.rho[ist, jst] = np.conj(molecule.states[ist].coef) * molecule.states[jst].coef
        PyMem_Free(coef)
    elif md.propagation == "density":
        for ist in range(nst):
            for jst in range(nst):
                molecule.rho[ist, jst] = rho[ist][jst]
    
        for ist in range(nst):
            PyMem_Free(rho[ist])
        PyMem_Free(rho)

    PyMem_Free(energy)
    PyMem_Free(energy_old)
    for ist in range(nst):
        PyMem_Free(nacme[ist])
        PyMem_Free(nacme_old[ist])
    PyMem_Free(nacme)
    PyMem_Free(nacme_old)

