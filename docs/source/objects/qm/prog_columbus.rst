
Columbus	
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^	

Columbus :cite:`Lischka2011` is one of famous open-source software for high-level *ab initio*
quantum calculation. Similar with other softwares, it can do various types of fundamental quantum
calculation. However, the major competitiveness of Columbus compared to other softwares is that
it is mainly designed to compute multireference calculations on electonic ground and excited states.
This feature is indeed well suited for dynamics in UNI-xMD, it is implemented for various types of dynamics.
In the current version of UNI-xMD, only CASSCF method is available.

- CASSCF is complete active space self-consistent field method. It provides analytical gradients as
  well as nonadiabatic couplings, thus it can be used for excited state molecular dynamics.

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| CASSCF | o    | o  | o  | o   |
+--------+------+----+----+-----+

CASSCF
"""""""""""""""""""""""""""""""""""""

.. note:: Our interface script is generated with 7.0 version of Columbus program.
   Here, you should refer to manual of Columbus program if you want to see detailed
   lists for **basis_set** variable.

+------------------------+-----------------------------------------------------+----------------+
| Keywords               | Work                                                | Default        |
+========================+=====================================================+================+
| **molecule**           | molecular object                                    |                |
| (:class:`Molecule`)    |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **basis_set**          | basis set information                               | *'6-31g\*'*    |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **memory**             | allocatable memory in the calculations              | *'500'*        |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **guess**              | initial guess for MCSCF method                      | *'hf'*         |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **guess_file**         | initial guess file                                  | *'./mocoef'*   |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **scf_en_tol**         | energy convergence for SCF iterations               | *9*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **scf_max_iter**       | maximum number of SCF iterations                    | *40*           |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **mcscf_en_tol**       | energy convergence for MCSCF iterations             | *8*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **mcscf_max_iter**     | maximum number of MCSCF iterations                  | *100*          |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **cpscf_grad_tol**     | gradient tolerance for CP-MCSCF equations           | *6*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **cpscf_max_iter**     | maximum number of iterations for CP-MCSCF equations | *100*          |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **active_elec**        | number of electrons in active space                 | *2*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **active_orb**         | number of orbitals in active space                  | *2*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **qm_path**            | path for QM binary                                  | *'./'*         |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **nthreads**           | number of threads in the calculations               | *1*            |
| *(integer)*            |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+
| **version**            | version of Molpro program                           | *'7.0'*        |
| *(string)*             |                                                     |                |
+------------------------+-----------------------------------------------------+----------------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **basis_set** *(string)* - Default: *6-31g**

 Basis set for calculation. In PyUNIxMD code, currently 10 basis are supported; {cc-pvdz, cc-pvtz, cc-pvqz, 3-21g*, 3-21+g*, 6-31g, 6-31g*, 6-31+g*, 6-311g*, 6-311+g*}.

\

- **memory** *(string)* - Default: *500*

 Total memory used for calculation. unit is 'm'

\

- **guess** *(string)* - Default: *hf*

 Initial guess method for MCSCF method

 + 'hf': Using HF orbitals as initial guess of CASSCF method for every time step.
 + 'read': Reads mocoef file generated from previous step as initial guess.

\

- **guess_file** *(string)* - Default: *./mocoef*

 Initial molecular orbital file for MCSCF method. It will be used as initial guess for CASSCF calculation in t=0.0s. This can be obtained from other MCSCF calculations.

\


- **scf_en_tol** *(integer)* - Default: *9*

 Energy threshold for SCF iterations. 

\


- **scf_max_iter** *(integer)* - Default: *40*

 Maximum number of SCF iterations.

\


- **mcscf_en_tol** *(integer)* - Default: *8*

 Energy convergence for MCSCF iterations

\

- **mcscf_max_iter** *(integer)* - Default: *100*

 Maximum number of MCSCF iterations

\

- **cpscf_grad_tol** *(integer)* - Default: *6*

 Gradient tolerance for CP-MCSCF equations

\

- **cpscf_max_iter** *(integer)* - Default: *100*

 Maximum number of iterations for CP-MCSCF equations

\

- **active_elec** *(integer)* - Default: *2*

 Number of electrons included in active space. Currently, only closed shell system is supported. 
 Number of electrons included in doubly occupied orbitals are automatically calculated by total number of electrons and active electrons.

\

- **active_orb** *(integer)* - Default: *2*

 Number of orbitals in active space. Currently, orbital matching between different timesteps are not supported, so be careful when deal with degenerated system.

\

- **qm_path** *(string)* - Default: *./*

 Path for QM binary. Path must be include binary file itself (ex. /opt/Columbus7.0/Columbus)

\

- **nthreads** *(integer)* - Default: *1*

 Number of threads for calculation.

\

- **version** *(double)* - Default: *7.0*

 Version of Columbus program. PyUNIxMD Code is currently based on ver 7.0, may not support ver 5.9 or lower.

\

