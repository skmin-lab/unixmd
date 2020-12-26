
TeraChem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TeraChem :cite:`Ufimtsev2008_1,Ufimtsev2009_1,Ufimtsev2009_2,Ufimtsev2008_2,Titov2013,Song2016` is general
purpose quantum chemistry software designed to run on NVIDIA GPU
architectures under a 64-bit Linux operating system. It includes many functionalities
such as DFT or wave function based methods. Among them, SSR method is interfaced with
current version of UNI-xMD.

- In general, spin-restricted ensemble-referenced Kohn-Sham (REKS) method can be classified
  as single-state REKS, state-averaged REKS (SA-REKS) and state-interaction SA-REKS (SSR).
  In single-state REKS, only ground state is calculated and it can treat the multireference
  character. SA-REKS and SSR can calculate the excited state as well as ground state. The
  difference is that the state-interaction term is considered in SSR so that more accurate
  states can be generated. SSR can provide nonadiabatic couplings so it can be used for
  surface hopping or Ehrenfest dynamics.

+-------------------+------+----+----+-----+
|                   | BOMD | SH | Eh | nac |
+===================+======+====+====+=====+
| single-state REKS | o    | x  | x  | x   |
+-------------------+------+----+----+-----+
| SA-REKS           | o    | x  | x  | x   |
+-------------------+------+----+----+-----+
| SI-SA-REKS (SSR)  | o    | o  | o  | o   |
+-------------------+------+----+----+-----+

SSR
"""""""""""""""""""""""""""""""""""""

UNI-xMD automatically determines single-state REKS as BO interfaces for ground state BOMD.
When we include the excited states, SA-REKS or SSR methods can be exploited and these are
determined from the **use_ssr_state** argument.

.. note:: Our interface script is generated with 1.93, 1.99 version of TeraChem program.
   Here, you should refer to manual of TeraChem program if you want to see detailed
   lists for **basis_set**, **functional**, **precision** variables.

.. note:: Currently, **guess** variable reads the following two strings.
   One is *'dft'*, which uses DFT orbitlas as initial guess of REKS method for every time step.
   The other is **read**, which reads c0 file generated from previous step.
   If **guess_file** exists, then c0 file is used as initial guess at t = 0.0 s.


+-------------------------+---------------------------------------------+-------------+
| Keywords                | Work                                        | Default     |
+=========================+=============================================+=============+
| **molecule**            | molecular object                            |             |  
| (:class:`Molecule`)     |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **basis_set**           | basis set information                       | *'sto-3g'*  |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **functional**          | functional in the calculations              | *'hf'*      |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **precision**           | precision in the calculations               | *'dynamic'* |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **scf_rho_tol**         | wavefunction convergence for SCF iterations | *1E-2**     |
| *(double)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **scf_max_iter**        | maximum number of SCF iterations            | *300*       |
| *(integer)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **ssr22**               | use REKS(2,2) calculation?                  | *True*      |
| *(boolean)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **guess**               | initial guess for REKS SCF iterations       | *dft*       |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **guess_file**          | initial guess file                          | *'./c0'*    |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **reks_rho_tol**        | DIIS error for REKS SCF iterations          | *1E-6*      |
| *(double)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **reks_max_iter**       | maximum number of REKS SCF iterations       | *1000*      |
| *(integer)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **shift**               | level shifting value in REKS SCF iterations | *0.3*       |
| *(double)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **use_ssr_state**       | calculate SSR state, if not, treat SA-REKS  | *True*      |
| *(boolean)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **cpreks_grad_tol**     | gradient tolerance for CP-REKS equations    | *1E-6*      |
| *(double)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **cpreks_max_iter**     | maximum number of CP-REKS iterations        | *1000*      |
| *(integer)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **qm_path**             | path for QM binary                          | *'./'*      |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **ngpus**               | number of GPUs                              | *1*         |
| *(integer)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **gpu_id**              | ID of used GPUs                             | *1*         |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **version**             | version of TeraChem program                 | *1.93*      |
| *(double)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+

