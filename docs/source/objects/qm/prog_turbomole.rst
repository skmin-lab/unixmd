
TURBOMOLE
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Turbomole :cite:`Ahlrichs1989` is quantum chemical program package, initially developed
in the group of Prof. Dr. Reinhart Ahlrichs at the University of Karlsruhe and at the Forschungszentrum Karlsruhe.
(TD)DFT method is interfaced with current version of UNI-xMD.

- (TD)DFT provides analytical gradients, thus it can be used born-oppenhiemer molecular dynamics (BOMD).

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| TDDFT  | o    | x  | x  | x   |
+--------+------+----+----+-----+

(TD)DFT
"""""""""""""""""""""""""""""""""""""

.. note:: Our interface script is generated with 6.4 version of TM program.
   Here, you should refer to manual of TM program if you want to see detailed
   lists for **basis_set**, **functional** variable.


+---------------------+-------------------------------------------+----------------+
| Keywords            | Work                                      | Default        |
+=====================+===========================================+================+
| **molecule**        | molecular object                          |                |  
| (:class:`Molecule`) |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **functional**      | xc functional information                 | *'b-lyp'*      |
| *(string)*          |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **basis_set**       | basis set information                     | *'SV(P)'*      |
| *(string)*          |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **memory**          | allocatable memory in the calculations    | *50*           |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **scf_max_iter**    | maximum number of SCF iterations          | *50*           |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **scf_en_tol**      | energy convergence for SCF iterations     | *6*            |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **cis_max_iter**    | maximum number of CIS iterations          | *25*           |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **cis_en_tol**      | energy convergence for CIS iterations     | *6*            |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **qm_path**         | path for QM program                       | *'./'*         |
| *(string)*          |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **nthreads**        | number of threads in the calculations     | *1*            |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **version**         | version of Turbomole program              | *6.4*          |
| *(double)*          |                                           |                |
+---------------------+-------------------------------------------+----------------+

