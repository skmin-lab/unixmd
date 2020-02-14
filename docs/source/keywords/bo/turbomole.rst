
Turbomole (TM) :cite:`TM` is quantum chemical program package, initially developed in the group of Prof. Dr. Reinhart Ahlrichs at the University of Karlsruhe and at the Forschungszentrum Karlsruhe. 
(TD)DFT method is interfaced with current version of UNI-xMD. 

- (TD)DFT provides analytical gradients, thus it can be used born-oppenhiemer molecular dynamics (BOMD).

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| TDDFT  | o    | x  | x  | x   |
+--------+------+----+----+-----+

(TD)DFT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with 6.4 version of TM program.
   Here, you should refer to manual of TM program if you want to see detailed
   lists for **basis_set**, **functional** variable.

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

