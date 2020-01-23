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

Gaussian09 has been a standard program for electronic structure calculations of molecules.
The only BOMD using the DFT option is available with Gaussian09 in the current version of UNI-xMD,
because it doesn't explicitly provide with nonadiabatic coupling vectors. 
Numerical calculation of the coupling elements using the CI overlap is on progress, which allows the other dynamics options.

- (TD)DFT is used to provide with a potential energy and its gradient for a certain adiabatic state.

+---------+------+----+----+-----+
|         | BOMD | SH | Eh | nac |
+=========+======+====+====+=====+
| (TD)DFT | o    | x  | x  | x   |
+---------+------+----+----+-----+


(TD)DFT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with Revision A.02 version of Gaussian09 program.
   Please refer to the manual for the detailed lists for **basis_set** variable.

+----------------+------------------------------------------------+---------------+
| Keywords       | Work                                           | Default       |
+================+================================================+===============+
| basis_set      | basis set information                          | sto-3g        |
+----------------+------------------------------------------------+---------------+
| memory         | allocatable memory in the calculations         | 500m          |
+----------------+------------------------------------------------+---------------+
| functional     | the level of DFT theory                        | BLYP          |
+----------------+------------------------------------------------+---------------+
| G09_root_path  | path for G09 root                              | /opt/gaussian |
+----------------+------------------------------------------------+---------------+
| nthreads       | number of threads in the calculations          | 1             |
+----------------+------------------------------------------------+---------------+
| version        | version of Molpro program                      | Revision A.02 |
+----------------+------------------------------------------------+---------------+

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
