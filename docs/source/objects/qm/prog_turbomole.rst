
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

+---------------------+---------------------------------------------+----------------+
| Keywords            | Work                                        | Default        |
+=====================+=============================================+================+
| **molecule**        | Molecule object                             |                |
| (:class:`Molecule`) |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **functional**      | Exchange-correlation functional information | *'b-lyp'*      |
| *(string)*          |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **basis_set**       | Basis set information                       | *'SV(P)'*      |
| *(string)*          |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **memory**          | Allocatable memory in the calculations      | *50*           |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **scf_max_iter**    | Maximum number of SCF iterations            | *50*           |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **scf_en_tol**      | Energy convergence for SCF iterations       | *6*            |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **cis_max_iter**    | Maximum number of CIS iterations            | *25*           |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **cis_en_tol**      | Energy convergence for CIS iterations       | *6*            |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **qm_path**         | Path for QM binary                          | *'./'*         |
| *(string)*          |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **nthreads**        | Number of threads in the calculations       | *1*            |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **version**         | Version of Turbomole program                | *'6.4'*        |
| *(string)*          |                                             |                |
+---------------------+---------------------------------------------+----------------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **functional** *(string)* - Default: *'b-lyp'*

  This argument specifies exchange-correlation functional used in Turbomole calculation.
  These arguments are same as the original arguments of Turbomole.
  It is recommended to check a Turbomole manual for the detailed list of **functional**.

\

- **basis_set** *(string)* - Default: *'SV(P)'*

  This argument specifies basis sets used in Turbomole calculation.
  These arguments are same as the original arguments of Turbomole.
  It is recommended to check a Turbomole manual for the detailed list of **basis_set**.

\

- **memory** *(integer)* - Default: *50*

  This argument determines how much memory will be allocated in a QM calculation. The unit is MB.

\

- **scf_max_iter** *(integer)* - Default: *50*

  This argument determines maximum number of SCF iterations.

\

- **scf_en_tol** *(integer)* - Default: *6*

  This argument determines energy convergence threshold for SCF iterations. Convergence threshold is :math:`10^{-\textbf{scf_en_tol}}`.

\

- **cis_max_iter** *(integer)* - Default: *25*

  This argument determines maximum number of CIS iterations.

\

- **cis_en_tol** *(integer)* - Default: *6*

  This argument determines energy convergence threshold for CIS iterations. Convergence threshold is :math:`10^{-\textbf{scf_en_tol}}`.

\

- **qm_path** *(string)* - Default: *'./'*

  This argument designates path for QM binary file for the Turbomole.
  The `$TURBOMOLE` environment variable determines the directory where Turbomole is installed, not the binary files themselves (For example, `$TURBOMOLE` is '/my_disk/my_name/TURBOMOLE').
  Thus, **qm_path** must be a *'`$TURBOMOLE`'*, not a *'`$TURBOMOLE`/define'*. 

\

- **nthreads** *(integer)* - Default: *1*

  This argument specifies number of threads for QM calculation.

\

- **version** *(string)* - Default: *'6.4'*

  This argument determines version of Turbomole program. PyUNIxMD is currently based on 7.0 version of Turbomole program.

\

