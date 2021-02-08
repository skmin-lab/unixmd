
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
| **version**         | version of Turbomole program              | *'6.4'*        |
| *(string)*          |                                           |                |
+---------------------+-------------------------------------------+----------------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **functional** *(string)* - Default: *'b-lyp'*

 Functional for DFT calculation. For the detailed list, check TM manual.

\

- **basis_set** *(string)* - Default: *'SV(P)'*

 Basis set for calculation. For the detailed list, check TM manual.

\

- **memory** *(integer)* - Default: *50*

 Total memory used for calculation. unit is 'mb'

\

- **scf_max_iter** *(integer)* - Default: *50*

 Maximum number of SCF iterations.

\

- **scf_en_tol** *(integer)* - Default: *6*

 Energy convergence of SCF iterations. Default value is 1.0E-6.

\

- **cis_max_iter** *(integer)* - Default: *25*

 Maximum number of CIS iterations.

\

- **cis_en_tol** *(integer)* - Default: *6*

 Energy convergence of CIS iterations. Default value is 1.0E-6.

\

- **qm_path** *(string)* - Default: *'./'*

 Path for QM binary. Path must be include binary file itself (ex. /opt/TURBOMOLE/define)

\

- **nthreads** *(integer)* - Default: *1*

 Number of threads for calculation.

\

- **version** *(string)* - Default: *'6.4'*

 Version of Turbomole program. Our interface script is generated with 6.4 version of TM program.
 Here, you should refer to manual of TM program if you want to see detailed lists for **basis_set**, **functional** variable.

\

