
Gaussian09
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gaussian09 :cite:`Frisch2009` has been a standard program for electronic structure calculations of molecules.
The only BOMD using the DFT option is available with Gaussian09 in the current version of UNI-xMD,
because it doesn't explicitly provide with nonadiabatic coupling vectors.

- (TD)DFT is used to provide with a potential energy and its gradient for a certain adiabatic state.

+---------+------+----+----+-----+
|         | BOMD | SH | Eh | nac |
+=========+======+====+====+=====+
| (TD)DFT | o    | o  | x  | x   |
+---------+------+----+----+-----+

(TD)DFT
"""""""""""""""""""""""""""""""""""""

.. note:: Our interface script is generated with Revision A.02 version of Gaussian09 program.
   Please refer to the manual for the detailed lists for **basis_set** and **functional** variable.

.. note:: Currently, **guess** variable reads the following two strings.
   One is *'Harris'*, which uses the default option (diagonalization of Harris functional) is used 
   for the initial guess of SCF iterations for every time step.
   The other is **read**, which reads the checkpoint file generated from previous step.
   If **guess_file** exists, then the checkpoint file is used as initial guess at t = 0.0 s.

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
| **memory**            | allocatable memory in the calculations | *'1gb'*           |
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

  Sets the level of DFT theory (exchange-correlation functional). These arguments are the same with the original arguments in used in Gaussian09. For the detailed list, see the manual of the Gaussian09 program.

\

- **basis_set** *(string)* - Default: *'sto-3g'*

  Sets the basis set. These arguments are the same with the original arguments in used in Gaussian09. For the detailed list, see the manual of the Gaussian09 program.

\

- **memory** *(string)* - Default: *'1gb'*

  Sets allocatable memory in the calculations.

\

- **guess** *(string)* - Default: *'Harris'*

  Sets the initial guess method for the SCF iteration.

  + 'Harris': Diagonalizes the Harris functional :cite:`Harris1985` for the initial guess. This is the default for all DFT calculations in Gaussian09.

  + 'read': Reads the (previous) checkpoint file, of which the path is given by **guess_file**, for the initial guess. If the checkpoint file is not provided at the initial MD step, the 'Harris' method is used only once instead.

\

- **guess_file** *(string)* - Default: *'./g09.chk'*

  Sets the path of the initial guess file.

\

- **G09_root_path** *(string)* - Default: *'/opt/gaussian'*

  Sets the path for G09 root.

\

- **nthreads** *(string)* - Default: *'1'*

  Sets the number of threads in the calculations.

\

- **version** *(string)* - Default: *'Revision A.02'*

  Sets the version of Gaussian09 program.
