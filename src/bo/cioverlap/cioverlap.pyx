# cython: language_level=3
#import numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Free
#from cpython.complex cimport complex
#cimport numpy as np

cdef extern from "tdnac.c":
    void TD_NAC(int nst, double **ao_overlap, double **mo_coef_old, double **mo_coef_new, double ***ci_coef_old, double ***ci_coef_new)

def wf_overlap(theory, molecule):
    cdef:
        double **ao_overlap
        double **mo_coef_old
        double **mo_coef_new
        double ***ci_coef_old
        double ***ci_coef_new
        int ist, nst, iorb, jorb, norb, nocc, nvirt

    # Assign size variables
    nst = molecule.nst
    norb = theory.norb
    nocc = theory.nocc
    nvirt = theory.nvirt

    # Allocate NACME variables
    ao_overlap = <double**> PyMem_Malloc(norb * sizeof(double*))
    mo_coef_old = <double**> PyMem_Malloc(norb * sizeof(double*))
    mo_coef_new = <double**> PyMem_Malloc(norb * sizeof(double*))

    for iorb in range(norb):
        ao_overlap[iorb] = <double*> PyMem_Malloc(norb * sizeof(double))
        mo_coef_old[iorb] = <double*> PyMem_Malloc(norb * sizeof(double))
        mo_coef_new[iorb] = <double*> PyMem_Malloc(norb * sizeof(double))

    ci_coef_old = <double***> PyMem_Malloc(nst * sizeof(double**))
    ci_coef_new = <double***> PyMem_Malloc(nst * sizeof(double**))

    for ist in range(nst):
        ci_coef_old[ist] = <double**> PyMem_Malloc(nvirt * sizeof(double*))
        ci_coef_new[ist] = <double**> PyMem_Malloc(nvirt * sizeof(double*))

    for ist in range(nst):
        for iorb in range(nvirt):
            ci_coef_old[ist][iorb] = <double*> PyMem_Malloc(nocc * sizeof(double))
            ci_coef_new[ist][iorb] = <double*> PyMem_Malloc(nocc * sizeof(double))

    # Assign NACME variables from python to C
    for iorb in range(norb):
        for jorb in range(norb):
            ao_overlap[iorb][jorb] = theory.ao_overlap[iorb, jorb]
            mo_coef_old[iorb][jorb] = theory.mo_coef_old[iorb, jorb]
            mo_coef_new[iorb][jorb] = theory.mo_coef_new[iorb, jorb]

    for ist in range(nst):
        for iorb in range(nvirt):
            for jorb in range(nocc):
                ci_coef_old[ist][iorb][jorb] = theory.ci_coef_old[ist, iorb, jorb]
                ci_coef_new[ist][iorb][jorb] = theory.ci_coef_new[ist, iorb, jorb]

    # Calculate TDNAC term for CIoverlap
    TD_NAC(nst, ao_overlap, mo_coef_old, mo_coef_new, ci_coef_old, ci_coef_new)

    # Assign NACME variables from C to python
#    for iorb in range(norb):
#        for jorb in range(norb):
#            theory.ao_overlap[iorb, jorb] = ao_overlap[iorb][jorb]
#            theory.mo_coef_old[iorb, jorb] = mo_coef_old[iorb][jorb]
#            theory.mo_coef_new[iorb, jorb] = mo_coef_new[iorb][jorb]
#
#    for ist in range(nst):
#        for iorb in range(nvirt):
#            for jorb in range(nocc):
#                theory.ci_coef_old[ist, iorb, jorb] = ci_coef_old[ist][iorb][jorb]
#                theory.ci_coef_new[ist, iorb, jorb] = ci_coef_new[ist][iorb][jorb]

    # Deallocate NACME variables
    for iorb in range(norb):
        PyMem_Free(ao_overlap[iorb])
        PyMem_Free(mo_coef_old[iorb])
        PyMem_Free(mo_coef_new[iorb])

    PyMem_Free(ao_overlap)
    PyMem_Free(mo_coef_old)
    PyMem_Free(mo_coef_new)

    for ist in range(nst):
        for iorb in range(nvirt):
            PyMem_Free(ci_coef_old[ist][iorb])
            PyMem_Free(ci_coef_new[ist][iorb])

    for ist in range(nst):
        PyMem_Free(ci_coef_old[ist])
        PyMem_Free(ci_coef_new[ist])

    PyMem_Free(ci_coef_old)
    PyMem_Free(ci_coef_new)


