
Columbus:cite:`Columbus` is one of famous free software for high-level ab initio quantum calculation. Similar with
other software, it can do various types of fundamental quantum calculation. However, the major
competitiveness of Columbus compare to other software is, it is mainly designed for calculate
multi-reference calculations on electonic ground and excited states. This feature is indeed well suited
for dynamics in UNI-xMD, it is implemented for various types of dynamics. In the current version of UNI-xMD,
only CASSCF method is available.

- CASSCF is complete active space self-consistent field method. It provides analytical gradients as
  well as nonadiabatic couplings, thus it can be used for excited state molecular dynamics.

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| CASSCF | o    | o  | o  | o   |
+--------+------+----+----+-----+

CASSCF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with 7.0 version of Columbus program.
   Here, you should refer to manual of Columbus program if you want to see detailed
   lists for **basis_set** variable.

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| basis_set      | basis set information                          | 6-31g*  |
+----------------+------------------------------------------------+---------+
| memory         | allocatable memory in the calculations         | 500     |
+----------------+------------------------------------------------+---------+
| active_elec    | number of electrons in active space            | 2       |
+----------------+------------------------------------------------+---------+
| active_orb     | number of orbitals in active space             | 2       |
+----------------+------------------------------------------------+---------+
| qm_path        | path for QM binary                             | ./      |
+----------------+------------------------------------------------+---------+
| nthreads       | number of threads in the calculations          | 1       |
+----------------+------------------------------------------------+---------+
| version        | version of Molpro program                      | 7.0     |
+----------------+------------------------------------------------+---------+

