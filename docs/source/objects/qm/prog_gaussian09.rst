
Gaussian 09
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gaussian 09 :cite:`Frisch2009` has been a standard program for electronic structure calculations of molecules.

- The (TD)DFT implementation of Gaussian 09 supports only analytical gradients, not nonadiabatic couplings.
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

+-----------------------+---------------------------------------------+-------------------+
| Keywords              | Work                                        | Default           |
+=======================+=============================================+===================+
| **molecule**          | Molecule object                             |                   |  
| (:class:`Molecule`)   |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **functional**        | Exchange-correlation functional information | *'BLYP'*          |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **basis_set**         | Basis set information                       | *'sto-3g'*        |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **memory**            | Allocatable memory                          | *'1gb'*           |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **guess**             | Initial guess for SCF iterations            | *'Harris'*        |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **guess_file**        | Initial guess file                          | *'./g09.chk'*     |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **G09_root_path**     | Path for Gaussian 09 root                   | *'./'*            |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **nthreads**          | Number of threads in the calculations       | *1*               |
| *(integer)*           |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **version**           | Version of Gaussian 09 program              | *'Revision A.02'* |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+

Detailed description of the arguments
""""""""""""""""""""""""""""""""""""""""""

- **functional** *(string)* - Default: *'BLYP'*

  This argument specifies exchange-correlation functional used in Gaussian 09 calculation.
  The available options for this argument is same as the original arguments of Gaussian 09.
  It is recommended to check a Gaussian 09 manual for the detailed list of **functional**.

\

- **basis_set** *(string)* - Default: *'sto-3g'*

  This argument specifies a basis set used in Gaussian 09 calculation.
  The available options for this argument is same as the original arguments of Gaussian 09.
  It is recommended to check a Gaussian 09 manual for the detailed list of **basis_set**.
\

- **memory** *(string)* - Default: *'1gb'*

  This argument determines how much memory will be allocated in a QM calculation.

\

- **guess** *(string)* - Default: *'Harris'*

  This argument determines initial guess for (TD)DFT calculations.

  + *'Harris'*: Use the default method of Gaussian 09 (Diagonalizing the Harris functional :cite:`Harris1985`).
  + *'read'*: Use orbitals calculated at the previous time step as the initial guess for the (TD)DFT calculation.

\

- **guess_file** *(string)* - Default: *'./g09.chk'*

  The **guess_file** determines the name of file containing orbitals for the initial guess of orbitals for the (TD)DFT calculation at the first time step.
  This argument is effective only if **guess** = *'read'*.
  If the file does not exist, *'Harris'* option is requested for the initial guess for the (TD)DFT calculation.

\

- **G09_root_path** *(string)* - Default: *'./'*

  This argument designates a path for the Gaussian 09 root directory, that is, the top level directory (for example, '/my_disk/my_name/gaussian09/').

\

- **nthreads** *(integer)* - Default: *1*

  This argument specifies number of threads for QM calculation.

\

- **version** *(string)* - Default: *'Revision A.02'*

  This argument determines version of Gaussian 09 program. PyUNIxMD is currently based on Revision A.02 version of Gaussian 09 program.

\
