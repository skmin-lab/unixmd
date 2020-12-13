# cython: language_level=3
from cpython.mem cimport PyMem_Malloc, PyMem_Free

cdef extern from "tdnac.c":
    void TD_NAC(int istep, int nst, int nbasis, int nspin, int norb, int *nocc, int *nvirt, double dt, \
        int *orb_ini, int *orb_final, double **nacme, double **ao_overlap, \
        double ***mo_coef_old, double ***mo_coef_new, double ****ci_coef_old, double ****ci_coef_new)

def wf_overlap(qm, molecule, istep_py, dt_py):
    cdef:
        int *orb_ini
        int *orb_final
        int *nocc
        int *nvirt
        double **nacme
        double **ao_overlap
        double ***mo_coef_old
        double ***mo_coef_new
        double ****ci_coef_old
        double ****ci_coef_new

        int istep, ist, nst, nsp, nocc_max, nvirt_max, ibasis, jbasis, iorb, jorb, nbasis, norb, nspin
        double dt

    # Assign size variables
    dt = dt_py
    istep = istep_py
    nst = molecule.nst
    nbasis = qm.nbasis
    norb = qm.norb
    nspin = qm.nspin

    # Allocate NACME variables
    orb_ini = <int*> PyMem_Malloc(nspin * sizeof(int))
    orb_final = <int*> PyMem_Malloc(nspin * sizeof(int))
    nocc = <int*> PyMem_Malloc(nspin * sizeof(int))
    nvirt = <int*> PyMem_Malloc(nspin * sizeof(int))

    for nsp in range(nspin):
        nocc[nsp] = qm.nocc[nsp]
        nvirt[nsp] = qm.nvirt[nsp]
    
    if nspin > 1:
        if nocc[0] > nocc[1]:
            nocc_max = nocc[0]
            nvirt_max = nvirt[1]
        else:
            nocc_max = nocc[1]
            nvirt_max = nvirt[0]
    else:
        nocc_max = nocc[0]
        nvirt_max = nvirt[0]


    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))

    ao_overlap = <double**> PyMem_Malloc(nbasis * sizeof(double*))

    mo_coef_old = <double***> PyMem_Malloc(norb * sizeof(double**))
    mo_coef_new = <double***> PyMem_Malloc(norb * sizeof(double**))

    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))

    for ibasis in range(nbasis):
        ao_overlap[ibasis] = <double*> PyMem_Malloc(nbasis * sizeof(double))

    for iorb in range(norb):
        mo_coef_old[iorb] = <double**> PyMem_Malloc(nbasis * sizeof(double*))
        mo_coef_new[iorb] = <double**> PyMem_Malloc(nbasis * sizeof(double*))

    for iorb in range(norb):
        for ibasis in range(nbasis):
            mo_coef_old[iorb][ibasis] = <double*> PyMem_Malloc(nspin * sizeof(double))
            mo_coef_new[iorb][ibasis] = <double*> PyMem_Malloc(nspin * sizeof(double))


    ci_coef_old = <double****> PyMem_Malloc(nst * sizeof(double***))
    ci_coef_new = <double****> PyMem_Malloc(nst * sizeof(double***))

    for ist in range(nst):
        ci_coef_old[ist] = <double***> PyMem_Malloc(nocc_max * sizeof(double**))
        ci_coef_new[ist] = <double***> PyMem_Malloc(nocc_max * sizeof(double**))

    for ist in range(nst):
        for iorb in range(nocc_max):
            ci_coef_old[ist][iorb] = <double**> PyMem_Malloc(nvirt_max * sizeof(double*))
            ci_coef_new[ist][iorb] = <double**> PyMem_Malloc(nvirt_max * sizeof(double*))

    for ist in range(nst):
        for iorb in range(nocc_max):
            for ibasis in range(nvirt_max):
                ci_coef_old[ist][iorb][ibasis] = <double*> PyMem_Malloc(nspin * sizeof(double))
                ci_coef_new[ist][iorb][ibasis] = <double*> PyMem_Malloc(nspin * sizeof(double))


    # Assign NACME variables from python to C
    orb_ini[0] = qm.orb_ini[0]
    orb_final[0] = qm.orb_final[0]
    if nspin > 1:
        orb_ini[1] = qm.orb_ini[1]
        orb_final[1] = qm.orb_final[1]

    for ist in range(nst):
        for jst in range(nst):
            nacme[ist][jst] = 0.

    for ibasis in range(nbasis):
        for jbasis in range(nbasis):
            ao_overlap[ibasis][jbasis] = qm.ao_overlap[ibasis, jbasis]
    
    for iorb in range(norb):
        for ibasis in range(nbasis):
            for nsp in range(nspin):
                mo_coef_old[iorb][ibasis][nsp] = qm.mo_coef_old[iorb, ibasis, nsp]
                mo_coef_new[iorb][ibasis][nsp] = qm.mo_coef_new[iorb, ibasis, nsp]

    for ist in range(nst):
        for iorb in range(nocc_max):
            for jorb in range(nvirt_max):
                for nsp in range(nspin):
                    ci_coef_old[ist][iorb][jorb][nsp] = qm.ci_coef_old[ist, iorb, jorb, nsp]
                    ci_coef_new[ist][iorb][jorb][nsp] = qm.ci_coef_new[ist, iorb, jorb, nsp]

    # Calculate TDNAC term for CIoverlap
    TD_NAC(istep, nst, nbasis, nspin, norb, nocc, nvirt, dt, orb_ini, orb_final, nacme, \
        ao_overlap, mo_coef_old, mo_coef_new, ci_coef_old, ci_coef_new)

    # Assign NACME variables from C to python
    for ist in range(nst):
        for jst in range(nst):
             molecule.nacme[ist, jst] = nacme[ist][jst]

    for iorb in range(norb):
        for ibasis in range(nbasis):
            for nsp in range(nspin):
                qm.mo_coef_old[iorb, ibasis, nsp] = mo_coef_new[iorb][ibasis][nsp]

    for ist in range(nst):
        for iorb in range(nocc_max):
            for jorb in range(nvirt_max):
                for nsp in range(nspin):
                    qm.ci_coef_old[ist, iorb, jorb, nsp] = ci_coef_new[ist][iorb][jorb][nsp]

    # Deallocate NACME variables
    PyMem_Free(orb_ini)
    PyMem_Free(orb_final)
    PyMem_Free(nocc)
    PyMem_Free(nvirt)

    for ist in range(nst):
        PyMem_Free(nacme[ist])

    PyMem_Free(nacme)

    for ibasis in range(nbasis):
        PyMem_Free(ao_overlap[ibasis])

    for iorb in range(norb):
        for ibasis in range(nbasis):
            PyMem_Free(mo_coef_old[iorb][ibasis])
            PyMem_Free(mo_coef_new[iorb][ibasis])

    for iorb in range(norb):
        PyMem_Free(mo_coef_old[iorb])
        PyMem_Free(mo_coef_new[iorb])

    PyMem_Free(ao_overlap)
    PyMem_Free(mo_coef_old)
    PyMem_Free(mo_coef_new)

    for ist in range(nst):
        for iorb in range(nocc_max):
            for ibasis range(nvirt_max):
                PyMem_Free(ci_coef_old[ist][iorb][ibasis])
                PyMem_Free(ci_coef_new[ist][iorb][ibasis])

    for ist in range(nst):
        for iorb in range(nocc_max):
            PyMem_Free(ci_coef_old[ist][iorb])
            PyMem_Free(ci_coef_new[ist][iorb])

    for ist in range(nst):
        PyMem_Free(ci_coef_old[ist])
        PyMem_Free(ci_coef_new[ist])

    PyMem_Free(ci_coef_old)
    PyMem_Free(ci_coef_new)


