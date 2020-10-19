# cython: language_level=3
import numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.complex cimport complex
cimport numpy as np

cdef extern from "rk4_xf.c":
    void rk4(char *propagation, int nst, int nesteps, double dt, double complex *coef, \
             double complex **rho, double *energy, double *energy_old, double **nacme, \
             double **nacme_old, int nat, int nsp, bint *l_coh, double *wsigma, \
             double *mass, double **pos, double ***aux_pos, double ***phase)

def el_run(md, molecule):
    cdef:
        double *energy
        double *energy_old
        double **nacme
        double **nacme_old
        double complex *coef
        double complex **rho
        bytes py_bytes
        char *propagation_c
        double dt
        int ist, jst, nst, nesteps
        
        double *wsigma
        double *mass
        double **pos
        double ***aux_pos
        double ***phase
        int iat, aux_nat, aux_nsp
        bint *l_coh

    #
    nst = molecule.nst
    nesteps, dt = md.nesteps, md.dt
    energy = <double*> PyMem_Malloc(nst * sizeof(double))
    energy_old = <double*> PyMem_Malloc(nst * sizeof(double))
    nacme = <double**> PyMem_Malloc(nst * sizeof(double*))
    nacme_old = <double**> PyMem_Malloc(nst * sizeof(double*))

    aux_nat, aux_nsp = md.aux.nat, md.aux.nsp
    wsigma = <double*> PyMem_Malloc(aux_nat * sizeof(double))
    mass = <double*> PyMem_Malloc(aux_nat * sizeof(double))
    pos = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
    aux_pos = <double***> PyMem_Malloc(nst * sizeof(double**))
    phase = <double***> PyMem_Malloc(nst * sizeof(double**))
    l_coh = <bint*> PyMem_Malloc(nst * sizeof(bint))
    
    #
    for ist in range(nst):
        nacme[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        nacme_old[ist] = <double*> PyMem_Malloc(nst * sizeof(double))
        energy[ist] = molecule.states[ist].energy
        energy_old[ist] = molecule.states[ist].energy_old
        for jst in range(nst):
            nacme[ist][jst] = molecule.nacme[ist, jst]
            nacme_old[ist][jst] = molecule.nacme_old[ist, jst]

    for iat in range(aux_nat):
        wsigma[iat] = md.wsigma[iat]
    
    for iat in range(aux_nat):
        pos[iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
        mass[iat] = md.aux.mass[iat]
        for isp in range(aux_nsp):
            pos[iat][isp] = md.pos_0[iat, isp]
    for ist in range(nst):
        aux_pos[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        phase[ist] = <double**> PyMem_Malloc(aux_nat * sizeof(double*))
        l_coh[ist] = md.l_coh[ist]
        for iat in range(aux_nat):
            aux_pos[ist][iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
            phase[ist][iat] = <double*> PyMem_Malloc(aux_nsp * sizeof(double))
            for isp in range(aux_nsp):
                aux_pos[ist][iat][isp] = md.aux.pos[ist, iat, isp]
                phase[ist][iat][isp] = md.phase[ist, iat, isp]

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
        rk4(propagation_c, nst, nesteps, dt, coef, rho, energy, energy_old, nacme, nacme_old, \
            aux_nat, aux_nsp, l_coh, wsigma, mass, pos, aux_pos, phase)
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
    
    PyMem_Free(wsigma)

    PyMem_Free(mass)
    PyMem_Free(l_coh)
    for ist in range(nst):
        for iat in range(aux_nat):
            PyMem_Free(aux_pos[ist][iat])
            PyMem_Free(phase[ist][iat])
        PyMem_Free(aux_pos[ist])
        PyMem_Free(phase[ist])
    for iat in range(aux_nat):
        PyMem_Free(pos[iat])
    PyMem_Free(aux_pos)
    PyMem_Free(phase)
    PyMem_Free(pos)

