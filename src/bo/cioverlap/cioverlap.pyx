# cython: language_level=3
from cpython.mem cimport PyMem_Malloc, PyMem_Free

cdef extern from "tdnac.c":
    void TD_NAC(int istep, int nst, int norb, int nocc, int nvirt, double dt, double **nacme, double **ao_overlap, double **mo_coef_old, double **mo_coef_new, double ***ci_coef_old, double ***ci_coef_new)

def wf_overlap(theory, molecule, istep_py, dt_py):
    cdef:
        double **nacme
        double **ao_overlap
        double **mo_coef_old
        double **mo_coef_new
        double ***ci_coef_old
        double ***ci_coef_new
        double dt
        int istep, ist, nst, iorb, jorb, norb, nocc, nvirt

    # Assign size variables
    dt = dt_py
    istep = istep_py
    nst = molecule.nst
    norb = theory.norb
    nocc = theory.nocc
    nvirt = theory.nvirt

    # Allocate NACME variables
    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))

    ao_overlap = <double**> PyMem_Malloc(norb * sizeof(double*))
    mo_coef_old = <double**> PyMem_Malloc(norb * sizeof(double*))
    mo_coef_new = <double**> PyMem_Malloc(norb * sizeof(double*))

    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))

    for iorb in range(norb):
        ao_overlap[iorb] = <double*> PyMem_Malloc(norb * sizeof(double))
        mo_coef_old[iorb] = <double*> PyMem_Malloc(norb * sizeof(double))
        mo_coef_new[iorb] = <double*> PyMem_Malloc(norb * sizeof(double))

    ci_coef_old = <double***> PyMem_Malloc(nst * sizeof(double**))
    ci_coef_new = <double***> PyMem_Malloc(nst * sizeof(double**))

    for ist in range(nst):
        ci_coef_old[ist] = <double**> PyMem_Malloc(nocc * sizeof(double*))
        ci_coef_new[ist] = <double**> PyMem_Malloc(nocc * sizeof(double*))

    for ist in range(nst):
        for iorb in range(nocc):
            ci_coef_old[ist][iorb] = <double*> PyMem_Malloc(nvirt * sizeof(double))
            ci_coef_new[ist][iorb] = <double*> PyMem_Malloc(nvirt * sizeof(double))

    # Assign NACME variables from python to C
    for ist in range(nst):
        for jst in range(nst):
            nacme[ist][jst] = 0.

    for iorb in range(norb):
        for jorb in range(norb):
            ao_overlap[iorb][jorb] = theory.ao_overlap[iorb, jorb]
            mo_coef_old[iorb][jorb] = theory.mo_coef_old[iorb, jorb]
            mo_coef_new[iorb][jorb] = theory.mo_coef_new[iorb, jorb]

    for ist in range(nst):
        for iorb in range(nocc):
            for jorb in range(nvirt):
                ci_coef_old[ist][iorb][jorb] = theory.ci_coef_old[ist, iorb, jorb]
                ci_coef_new[ist][iorb][jorb] = theory.ci_coef_new[ist, iorb, jorb]

    # Calculate TDNAC term for CIoverlap
    TD_NAC(istep, nst, norb, nocc, nvirt, dt, nacme, ao_overlap, mo_coef_old, mo_coef_new, ci_coef_old, ci_coef_new)

    # Assign NACME variables from C to python
    for ist in range(nst):
        for jst in range(nst):
             molecule.nacme[ist, jst] = nacme[ist][jst]

    for iorb in range(norb):
        for jorb in range(norb):
            theory.mo_coef_old[iorb, jorb] = mo_coef_new[iorb][jorb]

    for ist in range(nst):
        for iorb in range(nocc):
            for jorb in range(nvirt):
                theory.ci_coef_old[ist, iorb, jorb] = ci_coef_new[ist][iorb][jorb]

    # Deallocate NACME variables
    for ist in range(nst):
        PyMem_Free(nacme[ist])

    PyMem_Free(nacme)

    for iorb in range(norb):
        PyMem_Free(ao_overlap[iorb])
        PyMem_Free(mo_coef_old[iorb])
        PyMem_Free(mo_coef_new[iorb])

    PyMem_Free(ao_overlap)
    PyMem_Free(mo_coef_old)
    PyMem_Free(mo_coef_new)

    for ist in range(nst):
        for iorb in range(nocc):
            PyMem_Free(ci_coef_old[ist][iorb])
            PyMem_Free(ci_coef_new[ist][iorb])

    for ist in range(nst):
        PyMem_Free(ci_coef_old[ist])
        PyMem_Free(ci_coef_new[ist])

    PyMem_Free(ci_coef_old)
    PyMem_Free(ci_coef_new)


