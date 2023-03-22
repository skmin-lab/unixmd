# cython: language_level=3
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.complex cimport complex
import numpy as np
cimport numpy as np

cdef extern from "rk4.c":
    void rk4(int pst, int nesteps, double dt, char *elec_object, int **get_d_ind, double **ham_d, \
        double **ham_d_old, double **nacme, double **nacme_old, double complex *coef)

def el_run(md, qed):
    cdef:
        char *elec_object_c
        int **get_d_ind
        double **ham_d
        double **ham_d_old
        double **nacme
        double **nacme_old
        double complex *coef_a, *coef_d
#        double complex **rho

        bytes py_bytes
        int ist, jst, nst, pst, nesteps, verbosity
        double dt

    # Assign size variables
    nst, pst = md.pol.nst, md.pol.pst
    nesteps, dt = md.nesteps, md.dt

    # Allocate variables
    get_d_ind = <int**> PyMem_Malloc(pst * sizeof(int*))

    ham_d = <double**> PyMem_Malloc(pst * sizeof(double*))
    ham_d_old = <double**> PyMem_Malloc(pst * sizeof(double*))

    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))

    for ist in range(pst):
        get_d_ind[ist] = <int*> PyMem_Malloc(2 * sizeof(int))

    for ist in range(pst):
        ham_d[ist] = <double*> PyMem_Malloc(pst * sizeof(double))
        ham_d_old[ist] = <double*> PyMem_Malloc(pst * sizeof(double))

    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))

    # Assign variables from python to C
    for ist in range(pst):
        for jst in range(2):
            get_d_ind[ist][jst] = qed.get_d_ind[ist, jst]

    for ist in range(pst):
        for jst in range(pst):
            ham_d[ist][jst] = qed.ham_d[ist, jst]
            ham_d_old[ist][jst] = qed.ham_d_old[ist, jst]

    for ist in range(nst):
        for jst in range(nst):
            nacme[ist][jst] = md.pol.nacme[ist, jst]
            nacme_old[ist][jst] = md.pol.nacme_old[ist, jst]

    # Debug related
    verbosity = md.verbosity

    # Assign coef or rho with respect to propagation scheme
    if (md.elec_object == "coefficient"):

        coef_d = <double complex*> PyMem_Malloc(pst * sizeof(double complex))

        for ist in range(pst):
            coef_d[ist] = md.pol.pol_states[ist].coef_d

#    elif (md.elec_object == "density"):
#
#        rho = <double complex**> PyMem_Malloc(nst * sizeof(double complex*))
#        for ist in range(nst):
#            rho[ist] = <double complex*> PyMem_Malloc(nst * sizeof(double complex))
#
#        for ist in range(nst):
#            for jst in range(nst):
#                rho[ist][jst] = md.mol.rho[ist, jst]

    py_bytes = md.elec_object.encode()
    elec_object_c = py_bytes

    # Propagate electrons depending on the propagator
    if (md.propagator == "rk4"):
        rk4(pst, nesteps, dt, elec_object_c, get_d_ind, ham_d, ham_d_old, nacme, nacme_old, coef_d)

    # Assign variables from C to python
    if (md.elec_object == "coefficient"):

        for ist in range(pst):
            md.pol.pol_states[ist].coef_d = coef_d[ist]

        for ist in range(pst):
            for jst in range(pst):
                md.pol.rho_d[ist, jst] = np.conj(md.pol.pol_states[ist].coef_d) * md.pol.pol_states[jst].coef_d

        PyMem_Free(coef_d)

#    elif (md.elec_object == "density"):
#
#        for ist in range(nst):
#            for jst in range(nst):
#                md.mol.rho[ist, jst] = rho[ist][jst]
#
#        for ist in range(nst):
#            PyMem_Free(rho[ist])
#        PyMem_Free(rho)

    # Debug
    if (verbosity >= 1):
        for ist in range(pst):
            ind_mol1 = get_d_ind[ist][0]
            md.dotpopnac_d[ist] = 0.
            for jst in range(pst):
                ind_mol2 = get_d_ind[jst][0]
                if (jst != ist):
                    md.dotpopnac_d[ist] -= 2. * (ham_d[ist][jst] * md.pol.rho_d.imag[jst, ist] \
                        + nacme[ind_mol1][ind_mol2] * md.pol.rho_d.real[jst, ist])

    # Deallocate variables
    for ist in range(2):
        PyMem_Free(get_d_ind[ist])

    for ist in range(pst):
        PyMem_Free(ham_d[ist])
        PyMem_Free(ham_d_old[ist])

    for ist in range(nst):
        PyMem_Free(nacme[ist])
        PyMem_Free(nacme_old[ist])

    PyMem_Free(get_d_ind)

    PyMem_Free(ham_d)
    PyMem_Free(ham_d_old)

    PyMem_Free(nacme)
    PyMem_Free(nacme_old)


