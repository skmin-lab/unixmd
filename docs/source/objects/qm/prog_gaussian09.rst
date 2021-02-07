
Gaussian 09
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gaussian 09 :cite:`Frisch2009` has been a standard program for electronic structure calculations of molecules.

- The TDDFT implementation of Gaussian 09 supports only analytical gradients, not nonadiabatic couplings.
  Instead, nonadiabatic coupling matrix elements (NACME) are calculated by using our wavefunction overlap 
  :cite:`Ryabinkin2015` routines. Thus, it can be used for adiabatic dynamics and surface hopping dynamics.

+---------+------+--------+----+-----+
|         | BOMD | SH(XF) | Eh | nac |
+=========+======+========+====+=====+
| (TD)DFT | o    | o      | x  | x   |
+---------+------+--------+----+-----+

(TD)DFT
"""""""""""""""""""""""""""""""""""""

.. note:: Our interface is tested with Revision A.02 version of Gaussian 09 program.

+-----------------------+----------------------------------------+-------------------+
| Keywords              | Work                                   | Default           |
+=======================+========================================+===================+
| **molecule**          | Molecular object                       |                   |  
| (:class:`Molecule`)   |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **functional**        | The level of DFT theory                | *'BLYP'*          |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **basis_set**         | Basis set information                  | *'sto-3g'*        |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **memory**            | Allocatable memory                     | *'1gb'*           |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **guess**             | Initial guess type for SCF iterations  | *'Harris'*        |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **guess_file**        | Initial guess file                     | *'./g09.chk'*     |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **G09_root_path**     | Path for Gaussian 09 root              | *'/opt/gaussian'* |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **nthreads**          | Number of threads in the calculations  | *1*               |
| *(integer)*           |                                        |                   |
+-----------------------+----------------------------------------+-------------------+
| **version**           | Version of Gaussian 09 program         | *'Revision A.02'* |
| *(string)*            |                                        |                   |
+-----------------------+----------------------------------------+-------------------+

Detailed description of the arguments
""""""""""""""""""""""""""""""""""""""""""

- **functional** *(string)* - Default: *'BLYP'*

  It specifies the level of DFT theory (exchange-correlation functional). Use the same arguments with the original ones from Gaussian 09. For the detailed list, see the manual of the Gaussian 09 program.

\

- **basis_set** *(string)* - Default: *'sto-3g'*

  It specifies the basis set. Use the same arguments with the original ones from Gaussian 09. For the detailed list, see the manual of the Gaussian 09 program.

\

- **memory** *(string)* - Default: *'1gb'*

  It specifies allocatable memory in the calculations.

\

- **guess** *(string)* - Default: *'Harris'*

  It determines the initial guess method for SCF iterations.

  + *'Harris'*: Diagonalizes the Harris functional :cite:`Harris1985` for the initial guess. This is the default for all DFT calculations in Gaussian 09.
  + *'read'*: Reads the (previous) checkpoint file, of which the path is given by **guess_file**, for the initial guess. If the checkpoint file is not provided at the initial MD step, the *'Harris'* method is used only once instead.

\

- **guess_file** *(string)* - Default: *'./g09.chk'*

  It specifies the path of the checkpoint file for initial guess of SCF iterations.

\

- **G09_root_path** *(string)* - Default: *'/opt/gaussian'*

  It specifies the path for the root directory.

\

- **nthreads** *(integer)* - Default: *1*

  It specifies the number of threads in the calculations.

\

- **version** *(string)* - Default: *'Revision A.02'*

  It specifies the version of Gaussian 09 program.
