
Columbus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Columbus :cite:`Lischka2011` is one of open-source software for high-level ab initio
quantum calculation. It is designed primarily for extended multi-reference (MR) calculations
on electronic ground and excited states of atoms and molecules.
In the current version of PyUNIxMD, only (SA-)CASSCF method is available.

- (SA-)CASSCF is state-averaged complete active space self-consistent field method. It provides analytical gradients as
  well as nonadiabatic couplings, thus it can be used for excited state molecular dynamics.

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
| **version**            | Version of Columbus program                         | *'7.0'*        |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **basis_set** *(string)* - Default: *'6-31g\*'*

  This argument specifies a basis set used in Columbus calculation.
  Not all basis sets are supported, so it is recommended to check a Columbus manual for the compatibility with PyUNIxMD.
  In PyUNIxMD, currently 10 basis sets are supported; {*'cc-pvdz'*, *'cc-pvtz'*, *'cc-pvqz'*, *'3-21g\*'*, *'3-21+g\*'*, *'6-31g'*, *'6-31g\*'*, *'6-31+g\*'*, *'6-311g\*'*, *'6-311+g\*'*}.

\

- **memory** *(integer)* - Default: *500*

  This argument determines how much memory will be allocated in a QM calculation. The unit is MB.

\

- **guess** *(string)* - Default: *'hf'*

  This argument determines initial guess method for (SA-)CASSCF method. 

  + *'hf'*: Initial orbitals for (SA-)CASSCF method are generated from the HF calculation.
  + *'read'*: Reads 'mocoef' file generated from previous step as initial guess.
    At t = 0.0 s, **guess_file** will be used as initial guess.

\

- **guess_file** *(string)* - Default: *'./mocoef'*

  This argument designates initial molecular orbital file for (SA-)CASSCF method. It is valid when **guess** = *'read'*.
  It will be used as initial guess for (SA-)CASSCF calculation in first MD step.

\

- **scf_en_tol** *(integer)* - Default: *9*

  This argument determines energy convergence threshold for SCF iterations. Convergence threshold is :math:`10^{-\textbf{scf_en_tol}}`.

\

- **scf_max_iter** *(integer)* - Default: *40*

  This argument determines maximum number of SCF iterations.

\

- **mcscf_en_tol** *(integer)* - Default: *8*

  This argument determines energy convergence threshold for (SA-)CASSCF iterations. Convergence threshold is :math:`10^{-\textbf{mcscf_en_tol}}`.

\

- **mcscf_max_iter** *(integer)* - Default: *100*

  This argument determines maximum number of (SA-)CASSCF iterations.

\

- **cpscf_grad_tol** *(integer)* - Default: *6*

  This arugment determines gradient convergence threshold for CP-CASSCF equations. Convergence threshold is :math:`10^{-\textbf{cpscf_grad_tol}}`.

\

- **cpscf_max_iter** *(integer)* - Default: *100*

  This argument determines maximum number of iterations for CP-CASSCF equations.

\

- **active_elec** *(integer)* - Default: *2*

  This argument determines number of electrons included in active space. Currently, only closed shell system is supported. 

\

- **active_orb** *(integer)* - Default: *2*

  This argument determines number of orbitals in active space.

\

- **qm_path** *(string)* - Default: *'./'*

  This argument designates a path for QM binary files for the Columbus.
  The `$COLUMBUS` environment variable determines the directory where Columbus is installed, not the binary files themselves (ex. `$COLUMBUS` = /opt/Columbus7.0/Columbus/).
  Thus, **qm_path** must be a *'`$COLUMBUS`'*, not a *'`$COLUMBUS`/runc'*.

\

- **version** *(string)* - Default: *'7.0'*

  This argument determines version of Columbus program. PyUNIxMD is currently based on 7.0 version.

\

