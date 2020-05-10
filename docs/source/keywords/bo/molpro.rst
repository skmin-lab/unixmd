
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with 2015.1 version of Molpro program.
   Here, you should refer to manual of Molpro program if you want to see detailed
   lists for **basis_set** variable.

.. note:: Currently, **guess** variable reads the following two strings.
   One is **hf**, which uses HF orbitlas as initial guess of CASSCF method for every time step.
   The other is **read**, which reads wf.wfu file generated from previous step.
   If wf.wfu file exists in **guess_path**, then wf.wfu file is used as initial guess at t = 0.0 s.

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| basis_set      | basis set information                          | sto-3g  |
+----------------+------------------------------------------------+---------+
| memory         | allocatable memory in the calculations         | 500m    |
+----------------+------------------------------------------------+---------+
| guess          | initial guess for MCSCF method                 | hf      |
+----------------+------------------------------------------------+---------+
| guess_path     | directory for initial guess file               | ./      |
+----------------+------------------------------------------------+---------+
| scf_max_iter   | maximum number of SCF iterations               | 20      |
+----------------+------------------------------------------------+---------+
| scf_en_tol     | energy convergence for SCF iterations          | 1E-8    |
+----------------+------------------------------------------------+---------+
| scf_rho_tol    | density convergence for SCF iterations         | 1E-6    |
+----------------+------------------------------------------------+---------+
| mcscf_max_iter | maximum number of MCSCF iterations             | 20      |
+----------------+------------------------------------------------+---------+
| mcscf_en_tol   | energy convergence for MCSCF iterations        | 1E-8    |
+----------------+------------------------------------------------+---------+
| mcscf_grad_tol | gradient convergence for MCSCF iterations      | 1E-6    |
+----------------+------------------------------------------------+---------+
| mcscf_step_tol | step length convergence for MCSCF iterations   | 1E-2    |
+----------------+------------------------------------------------+---------+
| active_elec    | number of electrons in active space            | 2       |
+----------------+------------------------------------------------+---------+
| active_orb     | number of orbitals in active space             | 2       |
+----------------+------------------------------------------------+---------+
| cpscf_grad_tol | gradient tolerance for CP-MCSCF equations      | 1E-7    |
+----------------+------------------------------------------------+---------+
| qm_path        | path for QM binary                             | ./      |
+----------------+------------------------------------------------+---------+
| nthreads       | number of threads in the calculations          | 1       |
+----------------+------------------------------------------------+---------+
| version        | version of Molpro program                      | 2015.1  |
+----------------+------------------------------------------------+---------+

