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
Terachem
=====================================
availability for each,...

.. note:: naming :(

DFT-REKS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
bomd

DFT-SA-REKS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
bomd, sh

DFT-SI-SA-REKS(2,2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
bomd, sh, eh

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
Turbomole
=====================================

Turbomole (TM) is quantum chemical program package, initially developed in the group of Prof. Dr. Reinhart Ahlrichs at the University of Karlsruhe and at the Forschungszentrum Karlsruhe. 
(TD)DFT method is interfaced with current version of UNI-XMD. 
In UNI-XMD, TM, which is ver. 6.4, is valid for BOMD only

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| TDDFT  | o    | x  | x  | x   |
+--------+------+----+----+-----+

(TD)DFT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with 6.4 version of TM program.
   Here, you should refer to manual of TM program if you want to see detailed
   lists for **basis_set** variable.
+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| functional     | xc functional information                      | b-lyp   |
+----------------+------------------------------------------------+---------+
| basis_set      | basis set information                          | sto-3g  |
+----------------+------------------------------------------------+---------+
| memory         | allocatable memory in the calculations         | 500m    |
+----------------+------------------------------------------------+---------+
| max_iter       | maximum number of SCF iterations               | 20      |
+----------------+------------------------------------------------+---------+
| scf_en_tol     | energy convergence for SCF iterations          | 1E-8    |
+----------------+------------------------------------------------+---------+
| qm_path        | path for QM program                            | ./      |
+----------------+------------------------------------------------+---------+
| qm_bin_path    | path for QM binary                             | ./      |
+----------------+------------------------------------------------+---------+
| qm_scripts_path| path for QM scripts                            | ./      |
+----------------+------------------------------------------------+---------+
| nthreads       | number of threads in the calculations          | 1       |
+----------------+------------------------------------------------+---------+
| version        | version of Turbomole program                   | 6.4     |
+----------------+------------------------------------------------+---------+

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
