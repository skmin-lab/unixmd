
Columbus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Columbus :cite:`Lischka2011` is one of open-source softwares for high-level ab initio
quantum calculation. It is designed primarily for extended multi-reference (MR) calculations
on electronic ground and excited states of atoms and molecules.
In the current version of PyUNIxMD, (SA-)CASSCF and MRCI methods are available.

- (SA-)CASSCF is the state-averaged complete active space self-consistent field method. It provides analytical gradients as
  well as nonadiabatic couplings, thus it can be used for nonadiabatic molecular dynamics.

- MRCI is the multireference configuration interaction method. Simillar with (SA-)CASSCF, it supports analytical gradient and nonadiabatic couplings calculation,
  thus excited states dynamics simulations are possible.  

+-------------+------+--------+----+-----+
|             | BOMD | SH(XF) | Eh | nac |
+=============+======+========+====+=====+
| (SA-)CASSCF | o    | o      | o  | o   |
+-------------+------+--------+----+-----+
| MRCI        | o    | o      | o  | o   |
+-------------+------+--------+----+-----+

(SA-)CASSCF
"""""""""""""""""""""""""""""""""""""

+------------------------+-----------------------------------------------------+----------------+
| Parameters             | Work                                                | Default        |
+========================+=====================================================+================+
| **molecule**           | Molecule object                                     |                |
| (:class:`Molecule`)    |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **basis_set**          | Basis set information                               | *'6-31g\*'*    |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **memory**             | Allocatable memory in the calculations              | *500*          |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **guess**              | Initial guess method for (SA-)CASSCF method         | *'hf'*         |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **guess_file**         | Initial guess file                                  | *'./mocoef'*   |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **scf_en_tol**         | Energy convergence for SCF iterations               | *9*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **scf_max_iter**       | Maximum number of SCF iterations                    | *40*           |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **mcscf_en_tol**       | Energy convergence for (SA-)CASSCF iterations       | *8*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **mcscf_max_iter**     | Maximum number of (SA-)CASSCF iterations            | *100*          |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **state_avg**          | Number of states to be averaged for (SA-)CASSCF     | *None*         |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **cpscf_grad_tol**     | Gradient tolerance for CP-CASSCF equations          | *6*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **cpscf_max_iter**     | Maximum number of iterations for CP-CASSCF equations| *100*          |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **active_elec**        | Number of electrons in active space                 | *2*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **active_orb**         | Number of orbitals in active space                  | *2*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **qm_path**            | Path for QM binary                                  | *'./'*         |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **version**            | Version of Columbus                                 | *'7.0'*        |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+

Example input for CASSCF
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

   qm = qm.columbus.CASSCF(molecule=mol, basis_set='6-31g*', guess='hf', \
       state_avg=None, active_elec=2, active_orb=2, qm_path='/opt/Columbus7.0/Columbus/')

   md = mqc.SHXF(molecule=mol, nsteps=100, nesteps=20, dt=0.5, unit_dt='au', \
       sigma=0.1, istate=1, hop_rescale='energy', elec_object='density')

   md.run(qm=qm)

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **basis_set** *(string)* - Default: *'6-31g\*'*

  This parameter specifies a basis set used in the Columbus calculation.
  Not all basis sets are supported, so it is recommended to check the Columbus manual for the compatibility with PyUNIxMD.
  In PyUNIxMD, currently 10 basis sets are supported; {*'cc-pvdz'*, *'cc-pvtz'*, *'cc-pvqz'*, *'3-21g\*'*, *'3-21+g\*'*, *'6-31g'*, *'6-31g\*'*, *'6-31+g\*'*, *'6-311g\*'*, *'6-311+g\*'*}.

\

- **memory** *(integer)* - Default: *500*

  This parameter determines how much memory will be allocated in the QM calculation. The unit is MB.

\

- **guess** *(string)* - Default: *'hf'*

  This parameter determines the initial guess method for the (SA-)CASSCF calculations. 

  + *'hf'*: Initial guess orbitals for the (SA-)CASSCF calculations are generated from the HF calculations.
  + *'read'*: Initial guess orbitals are read from the 'mocoef' file which contains the orbitals calculated at the previous time step.

\

- **guess_file** *(string)* - Default: *'./mocoef'*

  The **guess_file** determines the name of the file containing orbitals for the initial guess of orbitals for the (SA-)CASSCF calculation at the first MD step.
  This parameter is effective only if **guess** = *'read'*.
  If the file does not exist, *'hf'* option is applied for the initial guess for the (SA-)CASSCF calculation at the first MD step.

\

- **scf_en_tol** *(integer)* - Default: *9*

  SCF is considered converged when the energy error is less than :math:`10^{-\textbf{scf_en_tol}}`.

\

- **scf_max_iter** *(integer)* - Default: *40*

  This parameter determines the maximum number of SCF iterations.

\

- **mcscf_en_tol** *(integer)* - Default: *8*

  (SA-)CASSCF is considered converged when the energy error is less than :math:`10^{-\textbf{mcscf_en_tol}}`.

\

- **mcscf_max_iter** *(integer)* - Default: *100*

  This parameter determines the maximum number of The (SA-)CASSCF iterations.

\

- **state_avg** *(integer)* - Default: *None*

  This parameter determines the number of states to be averaged for (SA-)CASSCF.
  The actual calculation is performed based on **state_avg**, not ``molecule.nst``.
  If it is not determined by a user, **state_avg** = ``molecule.nst``.

\

- **cpscf_grad_tol** *(integer)* - Default: *6*

  CP-CASSCF is considered converged when the gradient error is less than :math:`10^{-\textbf{cpscf_grad_tol}}`.

\

- **cpscf_max_iter** *(integer)* - Default: *100*

  This parameter determines the maximum number of iterations for the CP-CASSCF equations.

\

- **active_elec** *(integer)* - Default: *2*

  This parameter determines the number of electrons included in the active space. Currently, only closed shell system is supported. 

\

- **active_orb** *(integer)* - Default: *2*

  This parameter determines the number of orbitals in the active space.

\

- **qm_path** *(string)* - Default: *'./'*

  This parameter designates the path for QM binary files for Columbus.
  The `$COLUMBUS` environment variable determines the directory where Columbus is installed, not the binary files themselves (For example, `$COLUMBUS` is '/my_disk/my_name/Columbus7.0/Columbus/').
  Thus, **qm_path** must be *'`$COLUMBUS`'*, not *'`$COLUMBUS`/runc'*.

\

- **version** *(string)* - Default: *'7.0'*

  This parameter determines the version of Columbus. PyUNIxMD is currently based on version 7.0.

MRCI
"""""""""""""""""""""""""""""""""""""

+------------------------+-----------------------------------------------------+----------------+
| Parameters             | Work                                                | Default        |
+========================+=====================================================+================+
| **molecule**           | Molecule object                                     |                |
| (:class:`Molecule`)    |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **basis_set**          | Basis set information                               | *'6-31g\*'*    |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **memory**             | Allocatable memory in the calculations              | *500*          |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **guess**              | Initial guess method for MRCI method                | *'hf'*         |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **guess_file**         | Initial guess file                                  | *'./mocoef'*   |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **scf_en_tol**         | Energy convergence for SCF iterations               | *9*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **scf_max_iter**       | Maximum number of SCF iterations                    | *40*           |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **mcscf_en_tol**       | Energy convergence for (SA-)CASSCF iterations       | *8*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **mcscf_max_iter**     | Maximum number of (SA-)CASSCF iterations            | *100*          |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **mrci_en_tol**        | Energy convergence for MRCI iterations              | *4*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **mrci_max_iter**      | Maximum number of MRCI iterations                   | *None*         |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **state_avg**          | Number of states to be averaged                     | *None*         |
| *(integer)*            | for (SA-)CASSCF and MRCI                            |                |
+------------------------+-----------------------------------------------------+----------------+
| **active_elec**        | Number of electrons in active space                 | *2*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **active_orb**         | Number of orbitals in active space                  | *2*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **frozen_core_orb**    | Number of frozen core orbitals in                   | *0*            |
| *(integer)*            | doubly occupied space                               |                |
+------------------------+-----------------------------------------------------+----------------+
| **frozen_virt_orb**    | Number of frozen virtual orbitals from the          | *0*            |
| *(integer)*            | highest unoccupied space                            |                |
+------------------------+-----------------------------------------------------+----------------+
| **cpscf_grad_tol**     | Gradient tolerance for CP-MRCI equations            | *6*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **cpscf_max_iter**     | Maximum number of iterations for CP-MRCI equations  | *100*          |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **qm_path**            | Path for QM binary                                  | *'./'*         |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **version**            | Version of Columbus                                 | *'7.0'*        |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+

Example input for MRCI
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

   qm = qm.columbus.MRCI(molecule=mol, basis_set='6-31g*', guess='hf', \
       state_avg=None, active_elec=2, active_orb=2, frozen_core_orb=1, \
       frozen_virt_orb=0, qm_path='/opt/Columbus7.0/Columbus/')

   md = mqc.SHXF(molecule=mol, nsteps=100, nesteps=20, dt=0.5, unit_dt='au', \
       sigma=0.1, istate=1, hop_rescale='energy', elec_object='density')

   md.run(qm=qm)

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **basis_set** *(string)* - Default: *'6-31g\*'*

  This parameter specifies a basis set used in the Columbus calculation.
  Not all basis sets are supported, so it is recommended to check the Columbus manual for the compatibility with PyUNIxMD.
  In PyUNIxMD, currently 10 basis sets are supported; {*'cc-pvdz'*, *'cc-pvtz'*, *'cc-pvqz'*, *'3-21g\*'*, *'3-21+g\*'*, *'6-31g'*, *'6-31g\*'*, *'6-31+g\*'*, *'6-311g\*'*, *'6-311+g\*'*}.

\

- **memory** *(integer)* - Default: *500*

  This parameter determines how much memory will be allocated in the QM calculation. The unit is MB.

\

- **guess** *(string)* - Default: *'hf'*

  This parameter determines the initial guess method for the (SA-)CASSCF calculations. 
  (SA-)CASSCF must be performed before MRCI, and the optimized orbitals are used in MRCI calculations.

  + *'hf'*: Initial guess orbitals for the (SA-)CASSCF calculations are generated from the HF calculations.
  + *'read'*: Initial guess orbitals for the (SA-)CASSCF calculations are read from the 'mocoef' file
    which contains the orbitals calculated at the previous time step.

\

- **guess_file** *(string)* - Default: *'./mocoef'*

  The **guess_file** determines the name of the file containing orbitals for the initial guess of orbitals for the MCSCF calculation at the first MD step.
  This parameter is effective only if **guess** = *'read'*.
  If the file does not exist, *'hf'* option is applied for the initial guess for the (SA-)CASSCF calculation at the first MD step.

\

- **scf_en_tol** *(integer)* - Default: *9*

  SCF is considered converged when the energy error is less than :math:`10^{-\textbf{scf_en_tol}}`.

\

- **scf_max_iter** *(integer)* - Default: *40*

  This parameter determines the maximum number of SCF iterations.

\

- **mcscf_en_tol** *(integer)* - Default: *8*

  (SA-)CASSCF is considered converged when the energy error is less than :math:`10^{-\textbf{mcscf_en_tol}}`.

\

- **mcscf_max_iter** *(integer)* - Default: *100*

  This parameter determines the maximum number of The (SA-)CASSCF iterations.

\

- **mrci_en_tol** *(integer)* - Default: *4*

  MRCI is considered converged when the energy error is less than :math:`10^{-\textbf{mrci_en_tol}}`.

\

- **mrci_max_iter** *(integer)* - Default: *None*

  This parameter determines the maximum number of The MRCI iterations.
  If it is not determined by a user, **mrci_max_iter** :math:`= 30\,\times` **state_avg**.

\

- **state_avg** *(integer)* - Default: *None*

  This parameter determines the number of states to be averaged for (SA-)CASSCF and MRCI.
  The actual calculation is performed based on **state_avg**, not ``molecule.nst``.
  If it is not determined by a user, **state_avg** = ``molecule.nst``.

\

- **active_elec** *(integer)* - Default: *2*

  This parameter determines the number of electrons included in the active space. Currently, only closed shell system is supported. 

\

- **active_orb** *(integer)* - Default: *2*

  This parameter determines the number of orbitals in the active space.

\

- **forzen_core_elec** *(integer)* - Default: *0*

  This parameter determines the number of frozen core electrons included in doubly occupied space.

\

- **frozen_virt_orb** *(integer)* - Default: *0*

  This parameter determines the number of frozen virtual orbitals from the highest unoccupied space.

\

- **cpscf_grad_tol** *(integer)* - Default: *6*

  CP-MRCI is considered converged when the gradient error is less than :math:`10^{-\textbf{cpscf_grad_tol}}`.

\

- **cpscf_max_iter** *(integer)* - Default: *100*

  This parameter determines the maximum number of iterations for the CP-MRCI equations.

\

- **qm_path** *(string)* - Default: *'./'*

  This parameter designates the path for QM binary files for Columbus.
  The `$COLUMBUS` environment variable determines the directory where Columbus is installed, not the binary files themselves (For example, `$COLUMBUS` is '/my_disk/my_name/Columbus7.0/Columbus/').
  Thus, **qm_path** must be *'`$COLUMBUS`'*, not *'`$COLUMBUS`/runc'*.

\

- **version** *(string)* - Default: *'7.0'*

  This parameter determines the version of Columbus. PyUNIxMD is currently based on version 7.0.

