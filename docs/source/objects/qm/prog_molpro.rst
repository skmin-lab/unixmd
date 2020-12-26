
Molpro
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Molpro :cite:`Werner2012` is a comprehensive system of ab initio programs for advanced molecular electronic structure
calculations, designed and maintained by many authors. It comprises efficient and well parallelized
programs for standard computational chemistry applications, such as DFT or many wave function based
methods. Among them, CASSCF method is interfaced with current version of UNI-xMD.

- CASSCF is complete active space self-consistent field method. Molpro supports analytical gradients as
  well as nonadiabatic couplings, thus it can be used for excited state molecular dynamics.

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| CASSCF | o    | o  | o  | o   |
+--------+------+----+----+-----+

CASSCF
"""""""""""""""""""""""""""""""""""""

.. note:: Our interface script is generated with 2015.1 version of Molpro program.
   Here, you should refer to manual of Molpro program if you want to see detailed
   lists for **basis_set** variable.

.. note:: Currently, **guess** variable reads the following two strings.
   One is *'hf'*, which uses HF orbitlas as initial guess of CASSCF method for every time step.
   The other is **read**, which reads wf.wfu file generated from previous step.
   If **guess_file** exists, then wf.wfu file is used as initial guess at t = 0.0 s.

+----------------------+------------------------------------------------+----------------+
| Keywords             | Work                                           | Default        |
+======================+================================================+================+
| **molecule**         | molecular object                               |                |  
| (:class:`Molecule`)  |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **basis_set**        | basis set information                          | *'sto-3g'*     |
| *(string)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **memory**           | allocatable memory in the calculations         | *'500m'*       |
| *(string)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **guess**            | initial guess for MCSCF method                 | *'hf'*         |
| *(string)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **guess_file**       | initial guess file                             | *'./wf.wfu'*   |
| *(string)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **scf_max_iter**     | maximum number of SCF iterations               | *20*           |
| *(integer)*          |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **scf_en_tol**       | energy convergence for SCF iterations          | *1E-8*         |
| *(double)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **scf_rho_tol**      | density convergence for SCF iterations         | *1E-6*         |
| *(double)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **mcscf_max_iter**   | maximum number of MCSCF iterations             | *20*           |
| *(integer)*          |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **mcscf_en_tol**     | energy convergence for MCSCF iterations        | *1E-8*         |
| *(double)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **mcscf_grad_tol**   | gradient convergence for MCSCF iterations      | *1E-6*         |
| *(double)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **mcscf_step_tol**   | step length convergence for MCSCF iterations   | *1E-2*         |
| *(double)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **active_elec**      | number of electrons in active space            | *2*            |
| *(integer)*          |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **active_orb**       | number of orbitals in active space             | *2*            |
| *(integer)*          |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **cpscf_grad_tol**   | gradient tolerance for CP-MCSCF equations      | *1E-7*         |
| *(double)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **qm_path**          | path for QM binary                             | *'./'*         |
| *(string)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **nthreads**         | number of threads in the calculations          | *1*            |
| *(integer)*          |                                                |                |
+----------------------+------------------------------------------------+----------------+
| **version**          | version of Molpro program                      | *2015.1*       |
| *(double)*           |                                                |                |
+----------------------+------------------------------------------------+----------------+

