# cython: language_level=3
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.complex cimport complex
import numpy as np
cimport numpy as np

cdef extern from "rk4_xf.c":
    void rk4(int nat, int ndim, int pst, int nesteps, int verbosity, double dt, \
        char *elec_object, bint *l_coh, double *mass, double *sigma, int **get_d_ind, \
        double **unitary, double **ham_d, double **ham_d_old, double **nacme, \
        double **nacme_old, double **pos, double ***aux_pos, double ***phase, \
        double *dotpopdec, double complex *coef_d, double **qmom)

cdef extern from "exponential_xf.c":
    void exponential(int nat, int ndim, int pst, int nesteps, int verbosity, double dt, \
        char *elec_object, bint *l_coh, double *mass, double *sigma, int **get_d_ind, \
        double **unitary, double **ham_d, double **ham_d_old, double **nacme, \
        double **nacme_old, double **pos, double ***aux_pos, double ***phase, \
        double *dotpopdec, double complex *coef_d, double **qmom)

def el_run(md, qed):
    cdef:
        char *elec_object_c
        bint *l_coh
        double *mass
        double *sigma
        int **get_d_ind
        double **unitary
        double **ham_d
        double **ham_d_old
        double **nacme
        double **nacme_old
        double **pos
        double **qmom
        double ***aux_pos
        double ***phase
        double *dotpopdec_d
        double complex *coef_d
#        double complex **rho

        bytes py_bytes
        int ist, jst, nst, pst, nesteps, iat, aux_nat, isp, aux_ndim, verbosity
        int ind_mol1, ind_mol2
        double dt

    # Assign size variables
    nst, pst = md.pol.nst, md.pol.pst
    nesteps, dt = md.nesteps, md.dt
    aux_nat, aux_ndim = md.aux.nat, md.aux.ndim

    # Allocate variables
    l_coh = <bint*> PyMem_Malloc(pst * sizeof(bint))
    mass = <double*> PyMem_Malloc(aux_nat * sizeof(double))
    sigma = <double*> PyMem_Malloc(aux_nat * sizeof(double))

    get_d_ind = <int**> PyMem_Malloc(pst * sizeof(int*))

    unitary = <double**> PyMem_Malloc(pst * sizeof(double*))
    ham_d = <double**> PyMem_Malloc(pst * sizeof(double*))
    ham_d_old = <double**> PyMem_Malloc(pst * sizeof(double*))

    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))

    pos = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
    qmom = <double**> PyMem_Malloc(aux_nat * sizeof(double*))

    aux_pos = <double***> PyMem_Malloc(pst * sizeof(double**))
    phase = <double***> PyMem_Malloc(pst * sizeof(double**))

    for ist in range(pst):
        get_d_ind[ist] = <int*> PyMem_Malloc(2 * sizeof(int))

    for ist in range(pst):
        unitary[ist] = <double*> PyMem_Malloc(pst * sizeof(double))
        ham_d[ist] = <double*> PyMem_Malloc(pst * sizeof(double))
        ham_d_old[ist] = <double*> PyMem_Malloc(pst * sizeof(double))

    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))

    for iat in range(aux_nat):
        pos[iat] = <double*> PyMem_Malloc(aux_ndim * sizeof(double))
        qmom[iat] = <double*> PyMem_Malloc(aux_ndim * sizeof(double))

    for ist in range(pst):
        aux_pos[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        phase[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        for iat in range(aux_nat):
            aux_pos[ist][iat] = <double*> PyMem_Malloc(aux_ndim * sizeof(double))
            phase[ist][iat] = <double*> PyMem_Malloc(aux_ndim * sizeof(double))
    
    # Assign variables from python to C
    for iat in range(aux_nat):
        mass[iat] = md.aux.mass[iat]
        sigma[iat] = md.sigma[iat]

    for ist in range(pst):
        for jst in range(2):
            get_d_ind[ist][jst] = qed.get_d_ind[ist, jst]

    for ist in range(pst):
        for jst in range(pst):
            unitary[ist][jst] = qed.unitary[ist, jst]
            ham_d[ist][jst] = qed.ham_d[ist, jst]
            ham_d_old[ist][jst] = qed.ham_d_old[ist, jst]

    for ist in range(nst):
        for jst in range(nst):
            nacme[ist][jst] = md.pol.nacme[ist, jst]
            nacme_old[ist][jst] = md.pol.nacme_old[ist, jst]

    for iat in range(aux_nat):
        for isp in range(aux_ndim):
            pos[iat][isp] = md.pos_0[iat, isp]

    for ist in range(pst):
        l_coh[ist] = md.l_coh[ist]
        for iat in range(aux_nat):
            for isp in range(aux_ndim):
                aux_pos[ist][iat][isp] = md.aux.pos[ist, iat, isp]
                phase[ist][iat][isp] = md.phase[ist, iat, isp]

    # Debug related
    verbosity = md.verbosity
    if (verbosity >= 1):
        dotpopdec_d = <double*> PyMem_Malloc(pst * sizeof(double))

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
        rk4(aux_nat, aux_ndim, pst, nesteps, verbosity, dt, elec_object_c, l_coh, mass, \
            sigma, get_d_ind, unitary, ham_d, ham_d_old, nacme, nacme_old, pos, aux_pos, \
            phase, dotpopdec_d, coef_d, qmom)

    elif (md.propagator == "expon"):
        exponential(aux_nat, aux_ndim, pst, nesteps, verbosity, dt, elec_object_c, l_coh, \
            mass, sigma, get_d_ind, unitary, ham_d, ham_d_old, nacme, nacme_old, pos, aux_pos, \
            phase, dotpopdec_d, coef_d, qmom)

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

        for ist in range(pst):
            md.dotpopdec_d[ist] = dotpopdec_d[ist]

    if (verbosity >= 2):
        for iat in range(aux_nat):
            for isp in range(aux_ndim):
                md.qmom[iat, isp] = qmom[iat][isp]

    # Deallocate variables
    for ist in range(pst):
        for iat in range(aux_nat):
            PyMem_Free(aux_pos[ist][iat])
            PyMem_Free(phase[ist][iat])

    for ist in range(2):
        PyMem_Free(get_d_ind[ist])

    for ist in range(pst):
        PyMem_Free(unitary[ist])
        PyMem_Free(ham_d[ist])
        PyMem_Free(ham_d_old[ist])

    for ist in range(nst):
        PyMem_Free(nacme[ist])
        PyMem_Free(nacme_old[ist])

    for iat in range(aux_nat):
        PyMem_Free(pos[iat])
        PyMem_Free(qmom[iat])

    for ist in range(pst):
        PyMem_Free(aux_pos[ist])
        PyMem_Free(phase[ist])

    PyMem_Free(l_coh)
    PyMem_Free(mass)
    PyMem_Free(sigma)

    PyMem_Free(get_d_ind)

    PyMem_Free(unitary)
    PyMem_Free(ham_d)
    PyMem_Free(ham_d_old)

    PyMem_Free(nacme)
    PyMem_Free(nacme_old)

    PyMem_Free(pos)
    PyMem_Free(qmom)

    PyMem_Free(aux_pos)
    PyMem_Free(phase)


