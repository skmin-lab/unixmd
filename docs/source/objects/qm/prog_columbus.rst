
Columbus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Columbus :cite:`Lischka2011` is one of open-source softwares for high-level ab initio
quantum calculation. It is designed primarily for extended multi-reference (MR) calculations
on electronic ground and excited states of atoms and molecules.
In the current version of PyUNIxMD, only (SA-)CASSCF method is available.

- (SA-)CASSCF is the state-averaged complete active space self-consistent field method. It provides analytical gradients as
  well as nonadiabatic couplings, thus it can be used for nonadiabatic molecular dynamics.

+-------------+------+--------+----+-----+
|             | BOMD | SH(XF) | Eh | nac |
+=============+======+========+====+=====+
| (SA-)CASSCF | o    | o      | o  | o   |
+-------------+------+--------+----+-----+

(SA-)CASSCF
"""""""""""""""""""""""""""""""""""""

+------------------------+-----------------------------------------------------+----------------+
| Keywords               | Work                                                | Default        |
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

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **basis_set** *(string)* - Default: *'6-31g\*'*

  This argument specifies a basis set used in the Columbus calculation.
  Not all basis sets are supported, so it is recommended to check the Columbus manual for the compatibility with PyUNIxMD.
  In PyUNIxMD, currently 10 basis sets are supported; {*'cc-pvdz'*, *'cc-pvtz'*, *'cc-pvqz'*, *'3-21g\*'*, *'3-21+g\*'*, *'6-31g'*, *'6-31g\*'*, *'6-31+g\*'*, *'6-311g\*'*, *'6-311+g\*'*}.

\

- **memory** *(integer)* - Default: *500*

  This argument determines how much memory will be allocated in the QM calculation. The unit is MB.

\

- **guess** *(string)* - Default: *'hf'*

  This argument determines the initial guess method for (SA-)CASSCF calculations. 

  + *'hf'*: Initial guess orbitals for (SA-)CASSCF calculations are generated from the HF calculations.
  + *'read'*: Initial guess orbitals are read from the 'mocoef' file which contains the orbitals calculated at the previous time step.

\

- **guess_file** *(string)* - Default: *'./mocoef'*

  The **guess_file** determines the name of the file containing orbitals for the initial guess of orbitals for the (SA-)CASSCF calculation at the first MD step.
  This argument is effective only if **guess** = *'read'*.
  If the file does not exist, *'hf'* option is applied for the initial guess for the (SA-)CASSCF calculation at the first MD step.

\

- **scf_en_tol** *(integer)* - Default: *9*

  SCF is considered converged when the energy error is less than :math:`10^{-\textbf{scf_en_tol}}`.

\

- **scf_max_iter** *(integer)* - Default: *40*

  This argument determines the maximum number of SCF iterations.

\

- **mcscf_en_tol** *(integer)* - Default: *8*

  (SA-)CASSCF is considered converged when the energy error is less than :math:`10^{-\textbf{mcscf_en_tol}}`.

\

- **mcscf_max_iter** *(integer)* - Default: *100*

  This argument determines the maximum number of (SA-)CASSCF iterations.

\

- **cpscf_grad_tol** *(integer)* - Default: *6*

  CP-CASSCF is considered converged when the gradient error is less than :math:`10^{-\textbf{cpscf_grad_tol}}`.

\

- **cpscf_max_iter** *(integer)* - Default: *100*

  This argument determines the maximum number of iterations for CP-CASSCF equations.

\

- **active_elec** *(integer)* - Default: *2*

  This argument determines the number of electrons included in the active space. Currently, only closed shell system is supported. 

\

- **active_orb** *(integer)* - Default: *2*

  This argument determines the number of orbitals in the active space.

\

- **qm_path** *(string)* - Default: *'./'*

  This argument designates a path for QM binary files for Columbus.
  The `$COLUMBUS` environment variable determines the directory where Columbus is installed, not the binary files themselves (For example, `$COLUMBUS` is '/my_disk/my_name/Columbus7.0/Columbus/').
  Thus, **qm_path** must be *'`$COLUMBUS`'*, not *'`$COLUMBUS`/runc'*.

\

- **version** *(string)* - Default: *'7.0'*

  This argument determines the version of Columbus. PyUNIxMD is currently based on version 7.0.

\

