
DFTB+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DFTB+ :cite:`Hourahine2020` is a fast and efficient versatile quantum mechanical simulation software package.
Using DFTB+ you can carry out quantum mechanical simulations similar to density functional
theory but in an approximate way, typically gaining around two orders of magnitude in
speed. (TD)DFTB and SSR methods are interfaced with current version of UNI-xMD.

- (TD)DFTB is time-dependent density-functional tight-binding method. DFTB+ supports only
  analytical gradients, not nonadiabatic couplings. Instead, nonadiabatic coupling matrix
  elements (NACME) is calculated by using our wavefunction overlap :cite:`Ryabinkin2015` routines. 
  Thus, it can be used for adiabatic dynamics and surface hopping dynamics.
  For Ehrenfest dynamics, it shows no energy conservation since it cannot calculate
  nonadiabatic couplings directly which is used for nuclear propagation.

- In general, spin-restricted ensemble-referenced Kohn-Sham (REKS) method can be classified
  as single-state REKS, state-averaged REKS (SA-REKS) and state-interaction SA-REKS (SSR).
  In single-state REKS, only ground state is calculated and it can treat the multireference
  character. SA-REKS and SSR can calculate the excited state as well as ground state. The
  difference is that the state-interaction term is considered in SSR so that more accurate
  states can be generated. SSR can provide nonadiabatic couplings so it can be used for
  surface hopping or Ehrenfest dynamics.

+-------------------+------+----+----+-----+
|                   | BOMD | SH | Eh | nac |
+===================+======+====+====+=====+
| (TD)DFTB          | o    | o  | x  | x   |
+-------------------+------+----+----+-----+
| single-state REKS | o    | x  | x  | x   |
+-------------------+------+----+----+-----+
| SA-REKS           | o    | x  | x  | x   |
+-------------------+------+----+----+-----+
| SI-SA-REKS (SSR)  | o    | o  | o  | o   |
+-------------------+------+----+----+-----+

.. note:: To run DFTB+ interfacing script, the information about maximum angular momentum is
   needed. In current version of UNI-xMD, the values for maximum angular momentum is included
   in **max_l** dictionary variable of dftbpar.py file. The user can add or modify the key and
   value to **max_l**. The following example shows addition of new information about Si atom.

.. code-block:: python

   from qm.dftbplus.dftbpar import max_l

   max_l["Si"] = "p" # add value of new atom
   max_l["C"] = "s" # modify value of already existing atom

(TD)DFTB
"""""""""""""""""""""""""""""""""""""

.. note:: Our interface script is generated with 19.1/20.1 version of DFTB+ program.
   Here, you should refer to manual of DFTB+ program if you want to see detailed
   lists for **lc_method**, **mixer**, **ex_symmetry** variables.

.. note:: Currently, **guess** variable reads the following two strings.
   One is *'h0'*, which uses zedo charges as initial guess of SCC term for every time step.
   The other is **read**, which reads charges.bin file generated from previous step.
   If **guess_file** exists, then charges.bin file is used as initial guess at t = 0.0 s.

.. note:: For **cell_length** variable, it reads a list variable consisted of 9 float elements,
   which correspond to cell lattice vectors. Similarly **k_point** variable reads a list
   consisted of 3 integer elements, which correspond to number of k points.

+------------------------+------------------------------------------------+--------------------+
| Keywords               | Work                                           | Default            |
+========================+================================================+====================+
| **molecule**           | molecular object                               |                    |  
| (:class:`Molecule`)    |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **scc**                | include self-consistent charge (SCC) scheme    | *True*             |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **scc_tol**            | energy convergence for SCC iterations          | *1E-6*             |
| *(double)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **scc_max_iter**       | maximum number of SCC iterations               | *100*              |
| *(integer)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **ocdftb**             | include onsite correction to SCC term          | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **lcdftb**             | include long-range corrected functional        | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **lc_method**          | algorithms for LC-DFTB                         | *'MatrixBased'*    |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **sdftb**              | include spin-polarisation scheme               | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **unpaired_elec**      | number of unpaired electrons                   | *0.0*              |
| *(double)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **guess**              | initial guess method for SCC scheme            | *'h0'*             |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **guess_file**         | initial guess file for charges                 | *'./charges.bin'*  |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **elec_temp**          | electronic temperature for Fermi-Dirac scheme  | *0.0*              |
| *(double)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **mixer**              | charge mixing method used in DFTB              | *'Broyden'*        |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **ex_symmetry**        | symmetry of excited state in TD-DFTB           | *'singlet'*        |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **k_point**            | number of k-point samplings                    | *3 \* [ 1 ]*       |
| *(integer, list)*      |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **periodic**           | use periodicity in the calculations            | *False*            |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **cell_length**        | the lattice vectors of periodic unit cell      | *9 \* [ 0.0 ]*     |
| *(double, list)*       |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **sk_path**            | path for slater-koster files                   | *'./'*             |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **install_path**       | path for DFTB+ install directory               | *'./'*             |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **mpi**                | use MPI parallelization                        | *False*            |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **mpi_path**           | path for MPI binary                            | *'./'*             |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **nthreads**           | number of threads in the calculations          | *1*                |
| *(integer)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **version**            | version of DFTB+ program                       | *'20.1'*           |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+

SSR
"""""""""""""""""""""""""""""""""""""

PyUNIxMD automatically determines single-state REKS as BO interfaces for ground state BOMD.
When we include the excited states, SA-REKS or SSR methods can be exploited and these are
determined from the **state_interactions** argument.

+------------------------+------------------------------------------------+---------------------+
| Keywords               | Work                                           | Default             |
+========================+================================================+=====================+
| **molecule**           | molecular object                               |                     |
| (:class:`Molecule`)    |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **scc**                | include self-consistent charge (SCC) scheme    | *True*              |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **scc_tol**            | energy convergence for SCC iterations          | *1E-6*              |
| *(double)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **scc_max_iter**       | maximum number of SCC iterations               | *1000*              |
| *(integer)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **ocdftb**             | include onsite correction to SCC term          | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **lcdftb**             | include long-range corrected functional        | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **lc_method**          | algorithms for LC-DFTB                         | *'MatrixBased'*     |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **ssr22**              | use REKS(2,2) calculation?                     | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **ssr44**              | use REKS(4,4) calculation?                     | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **guess**              | initial guess method for SCC scheme            | *'h0'*              |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **guess_file**         | initial guess file for eigenvectors            | *'./eigenvec.bin'*  |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **state_interactions** | include state-interaction terms to SA-REKS     | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **shift**              | level shifting value in SCC iterations         | *0.3*               |
| *(double)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **tuning**             | scaling factor for atomic spin constants       | *None*              |
| *(double, list)*       |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **cpreks_grad_alg**    | algorithms used in CP-REKS equations           | *'pcg'*             |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **cpreks_grad_tol**    | gradient tolerance for CP-REKS equations       | *1E-8*              |
| *(double)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **save_memory**        | save memory in cache used in CP-REKS equations | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **embedding**          | charge-charge embedding options in QM/MM       | *None*              |
| *(string)*             | method                                         |                     |
+------------------------+------------------------------------------------+---------------------+
| **periodic**           | use periodicity in the calculations            | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **cell_length**        | the lattice vectors of periodic unit cell      | *9 \* [ 0.0 ]*      |
| *(double, list)*       |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **sk_path**            | path for slater-koster files                   | *'./'*              |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **install_path**       | path for DFTB+ install directory               | *'./'*              |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **nthreads**           | number of threads in the calculations          | *1*                 |
| *(integer)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **version**            | version of DFTB+ program                       | *'20.1'*            |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+


Detailed description of arguments
''''''''''''''''''''''''''''''''''''
- **scc** *(boolean)*

  Includes self-consistent charge (SCC) scheme. This is a mandatory argument to use SSR calculation.

  + True: Uses a SCC-DFTB method.
  + False: Uses a DFTB method.

\

- **scc_tol** *(double)*

  Energy convergence for SCC iterations.

\

- **scc_max_iter** *(integer)*

  Maximum number of SCC iterations.

\

- **ocdftb** *(boolean)*

  Includes onsite-correction to SCC term in the DFTB method. This is currently experimental feature,
  and not implemented in SSR calculation.

  + True: Uses a OC-DFTB method.
  + False: Uses a SCC-DFTB method.

\

- **lcdftb** *(boolean)*

  Includes long-range corrected functional in the DFTB method. To deal with the excited states properly,
  it is recommended to use LC-DFTB method for SSR calculation.

  + True: Uses a LC-DFTB method.
  + False: Uses a SCC-DFTB method.

\

- **lc_method** *(string)*

  Detailed algorithms used in LC-DFTB. These arguments are same with the original arguments used in DFTB+.

  + 'Thresholded': Screening according to estimated magnitude of terms.
  + 'NeighbourBased': Uses a purely neighbour-list based algorithm.
  + 'MatrixBased': Uses a matrix-matrix multiplication based algorithm.

\

- **ssr22** *(boolean)*

  Uses SSR(2,2) calculation in the context of DFTB method. The detailed type of the REKS calculation is
  determined the number of states and **state_interactions** argument. If the number of states is one,
  the single-state REKS calculation is carried out. When the number of states is larger than one,
  the SA-REKS or SI-SA-REKS calculation is executed according to the **state_interactions** argument.

  + True: Uses a DFTB/SSR(2,2) method.
  + False: Do not use a DFTB/SSR(2,2) method.

\

- **ssr44** *(boolean)*

  Uses SSR(4,4) calculation in the context of DFTB method. The detailed type of the REKS calculation is
  determined the number of states and **state_interactions** argument. If the number of states is one,
  the single-state REKS calculation is carried out. When the number of states is larger than one,
  the SA-REKS or SI-SA-REKS calculation is executed according to the **state_interactions** argument.
  This is currently experimental feature and not implemented.

  + True: Uses a DFTB/SSR(4,4) method.
  + False: Do not use a DFTB/SSR(4,4) method.

\

- **guess** *(string)*

  Initial guess method for the SCC scheme. The 'read' option with DFTB/SSR method is supported in 20.2 version (or newer).

  + 'h0': Initial orbitals are generated from the diagonalization of non-SCC Hamiltonian.
  + 'read': Reads *eigenvec.bin* file generated from previous step. If **guess_file** exists,
then *eigenvec.bin* file is used as initial guess at t = 0.0 s.

\

- **guess_file** *(string)*

  Initial guess file for eigenvectors. Is it vaild when **guess** is 'read' option.

\

- **state_interactions** *(boolean)*

  Includes state-interaction terms to SA-REKS calculation. It is valid when the number of states is larger
  than one. In general, it generates more reliable adiabatic states.

  + True: Uses SI-SA-REKS states.
  + False: Uses SA-REKS states.

\

- **shift** *(double)*

  Level shifting value used in SCC iterations. It can be helpful to increase **Shift** when
  it is hard to converge the SCC iterations.

\

- **tuning** *(double, list)*

  Scaling factor for atomic spin constants. It must be used carefully.
  The list consists of the number of atomic species.

\

- **cpreks_grad_alg** *(string)*

  Algorithms used in CP-REKS equations.

  + 'pcg': Uses a preconditioned conjugate-gradient based algorithm. It is generally faster than other algorithms.
  + 'cg': Uses a conjugate-gradient based algorithm. It is slower than 'pcg', but it can be helpful for systems including a high symmetry.
  + 'direct': Uses a direct matrix-inversion multiplication algorithm.

\

- **cpreks_grad_tol** *(double)*

  Tolerance of the gradient used in CP-REKS equations. This is not used when **cpreks_grad_alg** is 'direct' option.

\

- **save_memory** *(boolean)*

  Saves memory in cache used in CP-REKS equations.

  + True: Some variables which needs large memory allocation are save in the memory. In general, this becomes faster option.
  + False: Do not save in cache. This option is recommended for large systems.

\

- **embedding** *(string)*

  Charge-charge embedding options used in QM/MM method. It is recommended option for the environments showing high polarity.

  + None: Do not use charge-charge embedding in QM/MM method.
  + 'mechanical': Uses a mechanical charge-charge embedding option. The interactions are treated as the energies between MM point charges.
  + 'electrostatic': Uses a electrostatic charge-charge embedding option. Point charges as one-electron terms are included in the Hamiltonian.

\

- **periodic** *(boolean)*

  Uses a periodicity in the calculation. Only :math:`\Gamma`-point calculation is supported with DFTB/SSR method.

  + True: Uses a periodicity in the calculation.
  + False: Consider only cluster in the calculation.

\

- **cell_length** *(double, list)*

  Cell lattice vectors of the periodic unit cell. The list consists of nine elements, which correspond to the :math:`a`, :math:`b`, and :math:`c` vectors, respectively.

\

- **sk_path** *(string)*

  Path for slaker-koster files.

\

- **install_path** *(string)*

  Path for DFTB+ install directory. In general, it becomes '$DFTB/install/', not '$DFTB/install/bin/'.

\

- **nthreads** *(integer)*

  Number of threads in the calculation.

\

- **version** *(string)*

  Version of DFTB+ program. DFTB/SSR method is supported in 20.1 version (or newer).

