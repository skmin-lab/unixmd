
TeraChem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TeraChem :cite:`Ufimtsev2008_1,Ufimtsev2009_1,Ufimtsev2009_2,Ufimtsev2008_2,Titov2013,Song2016` is general
purpose quantum chemistry software designed to run on NVIDIA GPU
architectures under a 64-bit Linux operating system. It includes many functionalities
such as DFT or wave function based methods. Among them, the SSR method is interfaced with
the current version of PyUNIxMD.

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
When we include the excited states, the SA-REKS or SSR methods can be exploited and these are
determined from the **l_state_interactions** parameter.

.. note:: In the case of the SSR method, the calculation is possible only when the number
   of states (``molecule.nst``) is smaller than 4 due to the limited active space.
   If you want to treat more excited states, then increase the active space.

+--------------------------+---------------------------------------------+-------------+
| Parameters               | Work                                        | Default     |
+==========================+=============================================+=============+
| **molecule**             | Molecule object                             |             |  
| (:class:`Molecule`)      |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **basis_set**            | Basis set information                       | *'sto-3g'*  |
| *(string)*               |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **functional**           | Exchange-correlation functional information | *'hf'*      |
| *(string)*               |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **precision**            | Precision in the calculations               | *'dynamic'* |
| *(string)*               |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **scf_wf_tol**           | Wavefunction convergence for SCF iterations | *1E-2*      |
| *(double)*               |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **scf_max_iter**         | Maximum number of SCF iterations            | *300*       |
| *(integer)*              |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **active_space**         | Active space for SSR calculation            | *2*         |
| *(integer)*              |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **guess**                | Initial guess for REKS SCF iterations       | *'dft'*     |
| *(string)*               |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **guess_file**           | Initial guess file                          | *'./c0'*    |
| *(string)*               |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **reks_diis_tol**        | DIIS error for REKS SCF iterations          | *1E-6*      |
| *(double)*               |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **reks_max_iter**        | Maximum number of REKS SCF iterations       | *1000*      |
| *(integer)*              |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **shift**                | Level shifting value in REKS SCF iterations | *0.3*       |
| *(double)*               |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **l_state_interactions** | Include state-interaction terms to SA-REKS  | *False*     |
| *(boolean)*              |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **cpreks_grad_tol**      | Gradient tolerance for CP-REKS equations    | *1E-6*      |
| *(double)*               |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **cpreks_max_iter**      | Maximum number of CP-REKS iterations        | *1000*      |
| *(integer)*              |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **root_path**            | Path for TeraChem root directory            | *'./'*      |
| *(string)*               |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **ngpus**                | Number of GPUs                              | *1*         |
| *(integer)*              |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **gpu_id**               | ID of used GPUs                             | *None*      |
| *(integer, list)*        |                                             |             |
+--------------------------+---------------------------------------------+-------------+
| **version**              | Version of TeraChem                         | *'1.93'*    |
| *(string)*               |                                             |             |
+--------------------------+---------------------------------------------+-------------+

Example input for SSR
''''''''''''''''''''''''''''''''''''

.. code-block:: python
   :linenos:

   from molecule import Molecule
   import qm, mqc

   geom = '''
   3
   example
   O  1.14  3.77  0.00  0.00  0.00  0.00
   H  2.11  3.77  0.00  0.00  0.00  0.00
   H  0.81  4.45  0.60  0.00  0.00  0.00
   '''

   mol = Molecule(geometry=geom, ndim=3, nstates=2, unit_pos='angs')

   qm = qm.terachem.SSR(molecule=mol, active_space=2, guess='dft', basis_set='sto-3g'\
       l_state_interactions=True, shift=0.3, root_path='/opt/terachem1.93/TeraChem/')

   md = mqc.SHXF(molecule=mol, nsteps=100, nesteps=20, dt=0.5, unit_dt='au', \
       sigma=0.1, istate=1, hop_rescale='energy', elec_object='density')

   md.run(qm=qm)

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **basis_set** *(string)* - Default: *'sto-3g'*

  This parameter specifies the basis set used in TeraChem.
  The available options of this parameter are the same as the original ones of TeraChem.
  It is recommended to check a TeraChem manual for the detailed list of **basis_set**.

\

- **functional** *(string)* - Default: *'hf'*

  This parameter specifies the exchange-correlation functional used in TeraChem.
  The available options of this parameter are same as the original ones of TeraChem.
  It is recommended to check a TeraChem manual for the detailed list of **functional**.

\

- **precision** *(string)* - Default: *'dynamic'*

  This parameter specifies a method to determine the accuracy of the evaluation of the integrals.
  The available options of this parameter are same as the original ones of TeraChem.
  It is recommended to check a TeraChem manual for the detailed list of **precision**.

\

- **scf_wf_tol** *(double)* - Default: *1E-2*

  SCF cycles are considered converged when the wavefunction error is less than **scf_wf_tol**.

\

- **scf_max_iter** *(integer)* - Default: *300*

  This parameter determines the maximum number of SCF iterations.

\

- **active_space** *(integer)* - Default: *2*

  This parameter specifies the active space for SSR calculation. Detailed types of the REKS calculation are
  automatically determined by ``molecule.nst`` and **l_state_interactions** parameters. If ``molecule.nst`` is *1*,
  the single-state REKS calculation is carried out. When ``molecule.nst`` is larger than *1*,
  the SA-REKS or the SI-SA-REKS calculation is executed according to the **l_state_interactions** parameter.
  Currently, only (2,2) space is available for SSR calculation.

  + *2*: The numbers of electrons and orbitals are 2 and 2, respectively.

\

- **guess** *(string)* - Default: *'dft'*

  This parameter determines the initial guess method for the SSR calculations.

  + *'dft'*: Initial guess orbitals for the SSR calculations are generated from the DFT calculations.
  + *'read'*: Initial guess orbitals are read from the 'c0' file which contains the orbitals calculated at the previous time step.

\

- **guess_file** *(string)* - Default: *'./c0'*

  The **guess_file** determines the name of the file containing orbitals for the initial guess of orbitals for the SSR calculation at the first MD step.
  This parameter is effective only if **guess** = *'read'*.
  If the file does not exist, *'dft'* option is requested for the initial guess for the SSR calculation at the first MD step.

\

- **reks_diis_tol** *(double)* - Default: *1E-6*

  The REKS SCF cycles are considered converged when the DIIS error is less than **reks_diis_tol**.

\

- **reks_max_iter** *(integer)* - Default: *1000*

  This parameter determines the maximum number of the REKS SCF iterations.

\

- **shift** *(double)* - Default: *0.3*

  This parameter specifies the level shifting value used in the REKS SCF iterations. It can be helpful to increase **shift** when
  it is hard to converge the SCC iterations.

\

- **l_state_interactions** *(boolean)* - Default: *False*

  When **l_state_interactions** is set to *True*, state-interaction terms are included so that the SI-SA-REKS states are generated.
  Otherwise, the SA-REKS states are obtained. It is valid when ``molecule.nst`` is larger
  than *1*. In general, it generates more reliable adiabatic states.

\

- **cpreks_grad_tol** *(double)* - Default: *1E-6*

  This parameter determines the tolerance used in the conjugate-gradient based algorithm for solving the CP-REKS equations.
  Sometimes, it can be helpful to use slightly loose tolerance for the stable molecular dynamics.
  In this case, *4E-6* is recommended for **cpreks_grad_tol**.

\

- **cpreks_max_iter** *(integer)* - Default: *1000*

  This parameter determines the maximum number of the CP-REKS iterations.

\

- **root_path** *(string)* - Default: *'./'*

  This parameter determines the path for the TeraChem root directory. The `$TeraChem` environment
  variable determines the directory where the licensing file can be found, i.e., '`$TeraChem`/license.dat'
  (For example, `$TeraChem` is '/my_disk/my_name/TeraChem/').
  Thus, **root_path** must be *'`$TeraChem`/'*, not *'`$TeraChem`/bin/'*.

\

- **ngpus** *(integer)* - Default: *1*

  This parameter determines the number of GPUs used in TeraChem.

\

- **gpu_id** *(integer, list)* - Default: *None*

  This parameter specifies the ID of used GPUs. If you want to use 2 GPUs with ID of 0 and 1,
  then put *[0, 1]* into **gpu_id**.
  The length of **gpu_id** should be same to **ngpus**

\

- **version** *(string)* - Default: *'1.93'*

  This parameter determines the version of TeraChem.
  PyUNIxMD is currently based on version 1.93 and 1.99 of TeraChem.

