
TeraChem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TeraChem :cite:`Ufimtsev2008_1,Ufimtsev2009_1,Ufimtsev2009_2,Ufimtsev2008_2,Titov2013,Song2016` is general
purpose quantum chemistry software designed to run on NVIDIA GPU
architectures under a 64-bit Linux operating system. It includes many functionalities
such as DFT or wave function based methods. Among them, SSR method is interfaced with
current version of PyUNIxMD.

- In general, spin-restricted ensemble-referenced Kohn-Sham (REKS) method can be classified
  as single-state REKS, state-averaged REKS (SA-REKS) and state-interaction SA-REKS (SSR).
  In single-state REKS, only ground state is calculated and it can treat the multireference
  character. SA-REKS and SSR can calculate the excited state as well as ground state. The
  difference is that the state-interaction term is considered in SSR so that more accurate
  states can be generated. SSR can provide nonadiabatic couplings so it can be used for
  surface hopping or Ehrenfest dynamics.

+-------------------+------+--------+----+-----+
|                   | BOMD | SH(XF) | Eh | nac |
+===================+======+========+====+=====+
| single-state REKS | o    | x      | x  | x   |
+-------------------+------+--------+----+-----+
| SA-REKS           | o    | x      | x  | x   |
+-------------------+------+--------+----+-----+
| SI-SA-REKS (SSR)  | o    | o      | o  | o   |
+-------------------+------+--------+----+-----+

SSR
"""""""""""""""""""""""""""""""""""""

PyUNIxMD automatically determines single-state REKS as BO interfaces for ground state BOMD.
When we include the excited states, SA-REKS or SSR methods can be exploited and these are
determined from the **state_interactions** argument.

+-------------------------+---------------------------------------------+-------------+
| Keywords                | Work                                        | Default     |
+=========================+=============================================+=============+
| **molecule**            | Molecular object                            |             |  
| (:class:`Molecule`)     |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **basis_set**           | Basis set information                       | *'sto-3g'*  |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **functional**          | Functional in the calculations              | *'hf'*      |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **precision**           | Precision in the calculations               | *'dynamic'* |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **scf_rho_tol**         | Wavefunction convergence for SCF iterations | *1E-2**     |
| *(double)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **scf_max_iter**        | Maximum number of SCF iterations            | *300*       |
| *(integer)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **ssr22**               | Use SSR(2,2) calculation?                   | *False*     |
| *(boolean)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **guess**               | Initial guess for REKS SCF iterations       | *'dft'*     |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **guess_file**          | Initial guess file                          | *'./c0'*    |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **reks_rho_tol**        | DIIS error for REKS SCF iterations          | *1E-6*      |
| *(double)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **reks_max_iter**       | Maximum number of REKS SCF iterations       | *1000*      |
| *(integer)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **shift**               | Level shifting value in REKS SCF iterations | *0.3*       |
| *(double)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **state_interactions**  | Include state-interaction terms to SA-REKS  | *False*     |
| *(boolean)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **cpreks_grad_tol**     | Gradient tolerance for CP-REKS equations    | *1E-6*      |
| *(double)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **cpreks_max_iter**     | Maximum number of CP-REKS iterations        | *1000*      |
| *(integer)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **qm_path**             | Path for QM binary                          | *'./'*      |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **ngpus**               | Number of GPUs                              | *1*         |
| *(integer)*             |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **gpu_id**              | ID of used GPUs                             | *'1'*       |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **version**             | Version of TeraChem program                 | *'1.93'*    |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **basis_set** *(string)* - Default: *'sto-3g'*

  Basis set information used in the calculations.
  These arguments are same with the original arguments in used in TeraChem.
  If you want to know the detailed list for basis sets, see the manual of the TeraChem program.

\

- **functional** *(string)* - Default: *'hf'*

  Exchange-correlation functional used in the calculations.
  These arguments are same with the original arguments in used in TeraChem.
  If you want to know the detailed list for functionals, see the manual of the TeraChem program.

\

- **precision** *(string)* - Default: *'dynamic'*

  Method to determine the accuracy of the evaluation of the integrals.
  These arguments are same with the original arguments in used in TeraChem.
  If you want to know the detailed list for functionals, see the manual of the TeraChem program.

\

- **scf_rho_tol** *(double)* - Default: *1E-2*

  Wavefunction convergence for the SCF iterations.

\

- **scf_max_iter** *(integer)* - Default: *300*

  Maximum number of the SCF iteractions.

\

- **ssr22** *(boolean)* - Default: *False*

  Uses SSR(2,2) calculation. When this sets to *True*, detailed type of the REKS calculation is
  automatically determined from the number of states and **state_interactions** argument. If the number of states is one,
  the single-state REKS calculation is carried out. When the number of states is larger than one,
  the SA-REKS or SI-SA-REKS calculation is executed according to the **state_interactions** argument.

\

- **guess** *(string)* - Default: *'dft'*

  Initial guess method for the REKS SCF iteration.

  + *'dft'*: Initial orbitals are generated from the DFT calculation with **scf_rho_tol** tolerance.
  + *'read'*: Reads "c0" file generated from previous step.
    If **guess_file** exists, then "c0" file is used as initial guess at t = 0.0 s.

\

- **guess_file** *(string)* - Default: *'./c0'*

  Initial guess file for eigenvectors. It is vaild when **guess** is *'read'* option.

\

- **reks_rho_tol** *(double)* - Default: *1E-6*

  DIIS error for the REKS SCF iterations.

\

- **reks_max_iter** *(integer)* - Default: *1000*

  Maximum number of the REKS SCF iteractions.

\

- **shift** *(double)* - Default: *0.3*

  Level shifting value used in REKS SCF iterations. It can be helpful to increase **shift** when
  it is hard to converge the REKS SCF iterations.

\

- **state_interactions** *(boolean)* - Default: *False*

  Includes state-interaction terms to SA-REKS calculation. If this sets to *True*, the SI-SA-REKS states are calculated.
  Otherwise, the SA-REKS states are obtained. It is valid when the number of states is larger
  than one. In general, it generates more reliable adiabatic states.

\

- **cpreks_grad_tol** *(double)* - Default: *1E-6*

  Tolerance used in the conjugate-gradient based algorithm for solving the CP-REKS equations.
  Sometimes, it can be helpful to use loose tolerance for the stable molecular dynamics.
  In this case, *4E-6* is recommended for the tolerance.

\

- **cpreks_max_iter** *(integer)* - Default: *1000*

  Maximum number of the CP-REKS iterations.

\

- **qm_path** *(string)* - Default: *'./'*

  Path for TeraChem binary.

\

- **ngpus** *(integer)* - Default: *1*

  Number of GPUs used in the calculations.

\

- **gpu_id** *(string)* - Default: *'1'*

  ID of used GPUs. If you want to use 2 GPUs with ID of 0 and 1, then put *'0 1'* into **gpu_id**.

\

- **version** *(string)* - Default: *'1.93'*

  Version of TeraChem program. Currently, 1.93 and 1.99 versions are interfaced with PyUNIxMD.

