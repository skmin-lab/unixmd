
Columbus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Columbus :cite:`Lischka2011` is one of famous open-source software for high-level *ab initio*
quantum calculation. Similar with other softwares, it can do various types of fundamental quantum
calculation. However, the major competitiveness of Columbus compared to other softwares is that
it is mainly designed to compute multireference calculations on electonic ground and excited states.
This feature is indeed well suited for dynamics in UNI-xMD, it is implemented for various types of dynamics.
In the current version of UNI-xMD, only CASSCF method is available.

- CASSCF is complete active space self-consistent field method. It provides analytical gradients as
  well as nonadiabatic couplings, thus it can be used for excited state molecular dynamics.

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| CASSCF | o    | o  | o  | o   |
+--------+------+----+----+-----+

CASSCF
"""""""""""""""""""""""""""""""""""""

.. note:: Our interface script is generated with 7.0 version of Columbus program.
   Here, you should refer to manual of Columbus program if you want to see detailed
   lists for **basis_set** variable.

.. note:: Currently, **guess** variable reads the following two strings.
   One is *'hf'*, which uses HF orbitlas as initial guess of CASSCF method for every time step.
   The other is *'read'*, which reads mocoef file generated from previous step.
   If **guess_file** exists, then mocoef file is used as initial guess at t = 0.0 s.

+------------------------+-----------------------------------------------------+----------------+
| Keywords               | Work                                                | Default        |
+========================+=====================================================+================+
| **molecule**           | molecular object                                    |                |  
| (:class:`Molecule`)    |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **basis_set**          | basis set information                               | *'6-31g\*'*    |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **memory**             | allocatable memory in the calculations              | *'500'*        |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **guess**              | initial guess for MCSCF method                      | *'hf'*         |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **guess_file**         | initial guess file                                  | *'./mocoef'*   |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **scf_en_tol**         | energy convergence for SCF iterations               | *9*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **scf_max_iter**       | maximum number of SCF iterations                    | *40*           |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **mcscf_en_tol**       | energy convergence for MCSCF iterations             | *8*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **mcscf_max_iter**     | maximum number of MCSCF iterations                  | *100*          |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **cpscf_grad_tol**     | gradient tolerance for CP-MCSCF equations           | *6*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **cpscf_max_iter**     | maximum number of iterations for CP-MCSCF equations | *100*          |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **active_elec**        | number of electrons in active space                 | *2*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **active_orb**         | number of orbitals in active space                  | *2*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **qm_path**            | path for QM binary                                  | *'./'*         |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **nthreads**           | number of threads in the calculations               | *1*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **version**            | version of Molpro program                           | *'7.0'*        |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+

