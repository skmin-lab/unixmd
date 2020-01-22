=====================================
Molpro
=====================================

Molpro is a comprehensive system of ab initio programs for advanced molecular electronic structure
calculations, designed and maintained by many authors. It comprises efficient and well parallelized
programs for standard computational chemistry applications, such as DFT or many wave function based
methods. Among them, CASSCF method is interfaced with current version of UNI-XMD.

- CASSCF is complete active space self-consistent field method. It provides analytical gradients as
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

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| basis_set      | basis set information                          | sto-3g  |
+----------------+------------------------------------------------+---------+
| memory         | allocatable memory in the calculations         | 500m    |
+----------------+------------------------------------------------+---------+
| max_iter       | maximum number of MCSCF iterations             | 20      |
+----------------+------------------------------------------------+---------+
| scf_en_tol     | energy convergence for MCSCF iterations        | 1E-8    |
+----------------+------------------------------------------------+---------+
| scf_grad_tol   | gradient convergence for MCSCF iterations      | 1E-6    |
+----------------+------------------------------------------------+---------+
| scf_step_tol   | step length convergence for MCSCF iterations   | 1E-2    |
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

=====================================
Gaussian09
=====================================

version, theory,...

(TD)DFT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
bomd

=====================================
TeraChem
=====================================

TeraChem is general purpose quantum chemistry software designed to run on NVIDIA GPU
architectures under a 64-bit Linux operating system. It includes many functionalities
such as DFT or wave function based methods. Among them, SSR method is interfaced with
current version of UNI-XMD.

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

UNI-XMD automatically determines single-state REKS as BO interfaces for ground state BOMD.
When we include the excited states, SA-REKS or SSR methods can be exploited and these are
determined from the **use_ssr_state** argument.

.. note:: Our interface script is generated with 1.92, 1.93 version of TeraChem program.
   Here, you should refer to manual of TeraChem program if you want to see detailed
   lists for **basis_set**, **functional**, **precision** variables.

+-------------------+------------------------------------------------+---------+
| Keywords          | Work                                           | Default |
+===================+================================================+=========+
| basis_set         | basis set information                          | sto-3g  |
+-------------------+------------------------------------------------+---------+
| functional        | functional in the calculations                 | hf      |
+-------------------+------------------------------------------------+---------+
| precision         | precision in the calculations                  | dynamic |
+-------------------+------------------------------------------------+---------+
| scf_tol           | energy convergence for SCF iterations          | 1E-2    |
+-------------------+------------------------------------------------+---------+
| max_scf_iter      | maximum number of SCF iterations               | 300     |
+-------------------+------------------------------------------------+---------+
| reks22            | use REKS(2,2) calculation?                     | yes     |
+-------------------+------------------------------------------------+---------+
| reks_scf_tol      | energy convergence for REKS SCF iterations     | 1E-6    |
+-------------------+------------------------------------------------+---------+
| reks_max_scf_iter | maximum number of REKS SCF iterations          | 1000    |
+-------------------+------------------------------------------------+---------+
| reks_diis         | DIIS acceleration in REKS SCF iterations       | yes     |
+-------------------+------------------------------------------------+---------+
| shift             | level shifting value in REKS SCF iterations    | 0.3     |
+-------------------+------------------------------------------------+---------+
| use_ssr_state     | calculate SSR state, if not, treat SA-REKS     | 1       |
+-------------------+------------------------------------------------+---------+
| cpreks_max_tol    | gradient tolerance for CP-REKS equations       | 1E-6    |
+-------------------+------------------------------------------------+---------+
| cpreks_max_iter   | maximum number of CP-REKS iterations           | 1000    |
+-------------------+------------------------------------------------+---------+
| qm_path           | path for QM binary                             | ./      |
+-------------------+------------------------------------------------+---------+
| ngpus             | number of GPUs                                 | 1       |
+-------------------+------------------------------------------------+---------+
| gpu_id            | ID of used GPUs                                | 1       |
+-------------------+------------------------------------------------+---------+
| version           | version of TeraChem program                    | 1.92    |
+-------------------+------------------------------------------------+---------+

=====================================
DFTB+
=====================================
version,...

.. note:: naming :(

DFTB-REKS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
bomd

DFTB-SA-REKS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
bomd, sh

DFTB-SI-SA-REKS(2,2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
bomd, sh, eh

=====================================
Model
=====================================
..............?

Shin-Metiu
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
huhuh..

=====================================
Do It Yourself
=====================================

template ~~~ blah blah
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
template..~~? idunno
