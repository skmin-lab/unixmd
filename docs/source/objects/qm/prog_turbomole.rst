
TURBOMOLE
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Turbomole :cite:`Ahlrichs1989` is quantum chemical program package, initially developed
in the group of Prof. Dr. Reinhart Ahlrichs at the University of Karlsruhe and at the Forschungszentrum Karlsruhe.
(TD)DFT method is interfaced with current version of PyUNIxMD.

- (TD)DFT provides analytical gradients, thus it can be used born-oppenhiemer molecular dynamics (BOMD).

+---------+------+--------+----+-----+
|         | BOMD | SH(XF) | Eh | nac |
+=========+======+========+====+=====+
| (TD)DFT | o    | x      | x  | x   |
+---------+------+--------+----+-----+

(TD)DFT
"""""""""""""""""""""""""""""""""""""

+---------------------+-------------------------------------------+----------------+
| Keywords            | Work                                      | Default        |
+=====================+===========================================+================+
| **molecule**        | Molecular object                          |                |
| (:class:`Molecule`) |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **functional**      | XC functional information                 | *'b-lyp'*      |
| *(string)*          |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **basis_set**       | Basis set information                     | *'SV(P)'*      |
| *(string)*          |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **memory**          | Allocatable memory in the calculations    | *50*           |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **scf_max_iter**    | Maximum number of SCF iterations          | *50*           |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **scf_en_tol**      | Energy convergence for SCF iterations     | *6*            |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **cis_max_iter**    | Maximum number of CIS iterations          | *25*           |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **cis_en_tol**      | Energy convergence for CIS iterations     | *6*            |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **qm_path**         | Path for QM binary                        | *'./'*         |
| *(string)*          |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **nthreads**        | Number of threads in the calculations     | *1*            |
| *(integer)*         |                                           |                |
+---------------------+-------------------------------------------+----------------+
| **version**         | Version of Turbomole program              | *'6.4'*        |
| *(string)*          |                                           |                |
+---------------------+-------------------------------------------+----------------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **functional** *(string)* - Default: *'b-lyp'*

  This argument contains functional information about selected QM calculation.
  Not all functionals are supported depending on a QM program, so it is recommended to check a QM program manual for the compatibility with PyUNIxMD.

\

- **basis_set** *(string)* - Default: *'SV(P)'*

  This argument contains basis set information about selected QM calculation.
  Not all basis sets are supported depending on a QM program, so it is recommended to check a QM program manual for the compatibility with PyUNIxMD.

\

- **memory** *(integer)* - Default: *50*

  This argument contains how much memory will be used in a QM calculation. Basically, the unit is MB.

\

- **scf_max_iter** *(integer)* - Default: *50*

  This argument determines maximum number of SCF iterations.

\

- **scf_en_tol** *(integer)* - Default: *6*

  This argument determines energy threshold for SCF iterations. Convergence criteria is :math:`10^{-\textbf{scf_en_tol}}`.

\

- **cis_max_iter** *(integer)* - Default: *25*

  This argument determines maximum number of CIS iterations.

\

- **cis_en_tol** *(integer)* - Default: *6*

  This argument determines energy threshold for CIS iterations. Convergence criteria is :math:`10^{-\textbf{scf_en_tol}}`.

\

- **qm_path** *(string)* - Default: *'./'*

  This argument designates path for QM binary file for the selected QM calculation.
  Path must not include binary file itself. For example, **qm_path** = *'/opt/TURBOMOLE/'*.

\

- **nthreads** *(integer)* - Default: *1*

  This argument contains information of number of threads for QM calculation.

\

- **version** *(string)* - Default: *'6.4'*

  This argument determines version of Turbomole program. Our interface script is generated with 6.4 version of Turbomole program.

\

