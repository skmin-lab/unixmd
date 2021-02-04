
Gaussian 09
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gaussian 09 :cite:`Frisch2009` has been a standard program for electronic structure calculations of molecules.
The only BOMD using the DFT option is available with Gaussian 09 in the current version of UNI-xMD,
because it doesn't explicitly provide with nonadiabatic coupling vectors.

- (TD)DFT is used to provide with a potential energy and its gradient for a certain adiabatic state.

+---------+------+----+----+-----+
|         | BOMD | SH | Eh | nac |
+=========+======+====+====+=====+
| (TD)DFT | o    | o  | x  | x   |
+---------+------+----+----+-----+

(TD)DFT
"""""""""""""""""""""""""""""""""""""

.. note:: Our interface is tested with Revision A.02 version of Gaussian 09 program.

+-----------------------+----------------------------------------+-------------------+
| Keywords              | Work                                   | Default           |
+=======================+========================================+===================+
| **molecule**          | molecular object                       |                   |  
| (:class:`Molecule`)   |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **functional**        | the level of DFT theory                | *'BLYP'*          |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **basis_set**         | basis set information                  | *'sto-3g'*        |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **memory**            | allocatable memories                   | *'1gb'*           |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **guess**             | initial guess type for SCF iterations  | *'Harris'*        |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **guess_file**        | initial guess file                     | *'./g09.chk'*     |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **G09_root_path**     | path for G09 root                      | *'/opt/gaussian'* |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **nthreads**          | number of threads in the calculations  | *1*               |
| *(integer)*           |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **version**           | version of Gaussian09 program          | *'Revision A.02'* |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+

Detailed description of the arguments
""""""""""""""""""""""""""""""""""""""""""

- **functional** *(string)* - Default: *'BLYP'*

  It specifies the level of DFT theory (exchange-correlation functional). Use the same arguments with the original ones from Gaussian09. For the detailed list, see the manual of the Gaussian 09 program.

\

- **basis_set** *(string)* - Default: *'sto-3g'*

  It specifies the basis set. Use the same arguments with the original ones from Gaussian09. For the detailed list, see the manual of the Gaussian 09 program.

\

- **memory** *(string)* - Default: *'1gb'*

  It specifies allocatable memories in the calculations.

\

- **guess** *(string)* - Default: *'Harris'*

  It determines the initial guess method for SCF iterations.

  + *'Harris'*: Diagonalizes the Harris functional :cite:`Harris1985` for the initial guess. This is the default for all DFT calculations in Gaussian09.
  + *'read'*: Reads the (previous) checkpoint file, of which the path is given by **guess_file**, for the initial guess. If the checkpoint file is not provided at the initial MD step, the *'Harris'* method is used only once instead.

\

- **guess_file** *(string)* - Default: *'./g09.chk'*

  It specifies the path of the checkpoint file for initial guess of SCF iterations.

\

- **G09_root_path** *(string)* - Default: *'/opt/gaussian'*

  It specifies the path for the root directory.

\

- **nthreads** *(string)* - Default: *'1'*

  It specifies the number of threads in the calculations.

\

- **version** *(string)* - Default: *'Revision A.02'*

  It specifies the version of Gaussian 09 program.
