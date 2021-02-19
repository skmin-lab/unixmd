
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

.. note:: In the case of SSR method, the calculation is possible only when the number
   of states (``molecule.nst``) is smaller than 4 due to the limited active space.
   If you want to treat more excited states, then increase the active space.

+-------------------------+---------------------------------------------+-------------+
| Keywords                | Work                                        | Default     |
+=========================+=============================================+=============+
| **molecule**            | Molecule object                             |             |  
| (:class:`Molecule`)     |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **basis_set**           | Basis set information                       | *'sto-3g'*  |
| *(string)*              |                                             |             |
+-------------------------+---------------------------------------------+-------------+
| **functional**          | Exchange-correlation functional information | *'hf'*      |
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
| **reks_rho_tol**        | wavefunction error for REKS SCF iterations  | *1E-6*      |
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

  This argument specifies a basis set used in TeraChem.
  The available options of this argument are same as the original arguments of TeraChem.
  It is recommended to check a TeraChem manual for the detailed list of **basis_set**.

\

- **functional** *(string)* - Default: *'hf'*

  This argument specifies exchange-correlation functional used in TeraChem.
  The available options of this argument are same as the original arguments of TeraChem.
  It is recommended to check a TeraChem manual for the detailed list of **functional**.

\

- **precision** *(string)* - Default: *'dynamic'*

  This argument specifies a method to determine the accuracy of the evaluation of the integrals.
  The available options of this argument are same as the original arguments of TeraChem.
  It is recommended to check a TeraChem manual for the detailed list of **precision**.

\

- **scf_rho_tol** *(double)* - Default: *1E-2*

  SCF cycles are considered converged when the wavefunction error is less than **scf_rho_tol**.

\

- **scf_max_iter** *(integer)* - Default: *300*

  This argument determines maximum number of SCF iterations.

\

- **ssr22** *(boolean)* - Default: *False*

  When **ssr22** is set to *True*, SSR(2,2) calculation is carried out, and detailed type of the REKS calculation is
  automatically determined from ``molecule.nst`` and **state_interactions** arguments. If ``molecule.nst`` is one,
  the single-state REKS calculation is carried out. When ``molecule.nst`` is larger than one,
  the SA-REKS or SI-SA-REKS calculation is executed according to the **state_interactions** argument.

\

- **guess** *(string)* - Default: *'dft'*

  This argument determines initial guess method for SSR method.

  + *'dft'*: Initial orbitals for SSR method are generated from the DFT calculation with **scf_rho_tol** tolerance.
  + *'read'*: Use orbitals calculated at the previous time step as the initial guess for the SSR calculation.

\

- **guess_file** *(string)* - Default: *'./c0'*

  The **guess_file** determines the name of file containing orbitals for
  the initial guess of orbitals for the SSR calculation at the first MD step.
  This argument is effective only if **guess** = *'read'*.
  If the file does not exist, *'dft'* option is requested for the initial guess for the SSR calculation.

\

- **reks_rho_tol** *(double)* - Default: *1E-6*

  REKS SCF cycles are considered converged when the wavefunction error is less than **reks_rho_tol**.

\

- **reks_max_iter** *(integer)* - Default: *1000*

  This argument determines maximum number of REKS SCF iterations.

\

- **shift** *(double)* - Default: *0.3*

  This argument specifies level shifting value used in REKS SCF iterations. It can be helpful to increase **shift** when
  it is hard to converge the SCC iterations.

\

- **state_interactions** *(boolean)* - Default: *False*

  When **state_interactions** is set to *True*, state-interaction terms are included so that SI-SA-REKS states are generated.
  Otherwise, the SA-REKS states are obtained. It is valid when ``molecule.nst`` is larger
  than one. In general, it generates more reliable adiabatic states.

\

- **cpreks_grad_tol** *(double)* - Default: *1E-6*

  This argument determines tolerance used in the conjugate-gradient based algorithm for solving the CP-REKS equations.
  Sometimes, it can be helpful to use slightly loose tolerance for the stable molecular dynamics.
  In this case, *4E-6* is recommended for **cpreks_grad_tol**.

\

- **cpreks_max_iter** *(integer)* - Default: *1000*

  This argument determines maximum number of CP-REKS iterations.

\

- **qm_path** *(string)* - Default: *'./'*

  This argument determines path for QM binary file for TeraChem. The `$TeraChem` environment
  variable determines the directory where the licensing file can be found, i.e. '`$TeraChem`/license.dat'
  (For example, `$TeraChem` is '/my_disk/my_name/TeraChem/').
  Thus, **qm_path** must be *'`$TeraChem`/bin/'*, not *'`$TeraChem`/'*.

\

- **ngpus** *(integer)* - Default: *1*

  This argument determines number of GPUs used in TeraChem.

\

- **gpu_id** *(string)* - Default: *'1'*

  This argument specifies the ID of used GPUs. If you want to use 2 GPUs with ID of 0 and 1,
  then put *'0 1'* into **gpu_id**.

\

- **version** *(string)* - Default: *'1.93'*

  This argument determines version of TeraChem program.
  PyUNIxMD is currently based on 1.93 and 1.99 versions of TeraChem program.

