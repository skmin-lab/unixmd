# cython: language_level=3
import numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.complex cimport complex
cimport numpy as np

cdef extern from "tdnac.c":
    void TD_NAC(int nst)
#    void TD_NAC(int nst, int nesteps, double dt, double complex **rho, double *energy, double *energy_old, double **nacme, double **nacme_old)

def wf_overlap(theory, molecule):
    cdef:
#        double *energy
#        double *energy_old
#        double **nacme
#        double **nacme_old
#        double complex **rho
#        int ist, jst, nst, nesteps
#        double dt
        int nst

    nst = molecule.nst

    TD_NAC(nst)
#    print (theory.scc_tol)

#    nesteps, dt = nesteps_py, dt_py
#
#    energy = <double*> PyMem_Malloc(nst * sizeof(double))
#    energy_old = <double*> PyMem_Malloc(nst * sizeof(double))
#    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
#    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))
#    rho = <double complex**> PyMem_Malloc(nst * sizeof(double complex*))
#
#    for ist in range(nst):
#        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
#        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
#        rho[ist] = <double complex*> PyMem_Malloc(nst * sizeof(double complex))
#
#    for ist in range(nst):
#        energy[ist] = molecule.states[ist].energy
#        energy_old[ist] = molecule.states[ist].energy_old
#        for jst in range(nst):
#            nacme[ist][jst] = molecule.nacme[ist, jst]
#            nacme_old[ist][jst] = molecule.nacme_old[ist, jst]
#            rho[ist][jst] = molecule.rho[ist, jst]
#
#    RK4_rho(nst, nesteps, dt, rho, energy, energy_old, nacme, nacme_old)
#
#    for ist in range(nst):
#        for jst in range(nst):
#            molecule.rho[ist, jst] = rho[ist][jst]
#
#    PyMem_Free(energy)
#    PyMem_Free(energy_old)
#    for ist in range(nst):
#        PyMem_Free(nacme[ist])
#        PyMem_Free(nacme_old[ist])
#        PyMem_Free(rho[ist])
#    PyMem_Free(nacme)
#    PyMem_Free(nacme_old)
#    PyMem_Free(rho)


