
DFTB+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DFTB+ :cite:`Hourahine2020` is a fast and efficient versatile quantum mechanical simulation software package.
Using DFTB+ you can carry out quantum mechanical simulations similar to density functional
theory but in an approximate way, typically gaining around two orders of magnitude in
speed. (TD)DFTB and SSR methods are interfaced with current version of PyUNIxMD.

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
   needed. In current version of PyUNIxMD, the values for maximum angular momentum are included
   in **max_l** dictionary variable in '`$PYUNIXMD`/src/qm/dftb/dftbpar.py' file.
   You can add or modify the **max_l** variable, see the following example.

.. code-block:: python

   from qm.dftbplus.dftbpar import max_l

   max_l["Si"] = "p" # add value of new Si atom
   max_l["C"] = "s" # modify value of already existing C atom

(TD)DFTB
"""""""""""""""""""""""""""""""""""""

+------------------------+------------------------------------------------+--------------------+
| Keywords               | Work                                           | Default            |
+========================+================================================+====================+
| **molecule**           | Molecular object                               |                    |  
| (:class:`Molecule`)    |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **scc**                | Include self-consistent charge (SCC) scheme    | *True*             |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **scc_tol**            | Stopping criteria for the SCC iterations       | *1E-6*             |
| *(double)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **scc_max_iter**       | Maximum number of SCC iterations               | *100*              |
| *(integer)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **ocdftb**             | Include onsite correction to SCC term          | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **lcdftb**             | Include long-range corrected functional        | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **lc_method**          | Algorithms for LC-DFTB                         | *'MatrixBased'*    |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **sdftb**              | Include spin-polarisation scheme               | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **unpaired_elec**      | Number of unpaired electrons                   | *0.0*              |
| *(double)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **guess**              | Initial guess method for SCC scheme            | *'h0'*             |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **guess_file**         | Initial guess file for charges                 | *'./charges.bin'*  |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **elec_temp**          | Electronic temperature for Fermi-Dirac scheme  | *0.0*              |
| *(double)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **mixer**              | Charge mixing method used in DFTB              | *'Broyden'*        |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **ex_symmetry**        | Symmetry of excited state in TD-DFTB           | *'singlet'*        |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **e_window**           | Energy window for TD-DFTB. Increases efficiency| *0.0*              |
| *(double)*             | of NACME calculation                           |                    |
+------------------------+------------------------------------------------+--------------------+
| **k_point**            | Number of k-point samplings                    | *3 \* [ 1 ]*       |
| *(integer, list)*      |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **periodic**           | Use periodicity in the calculations            | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **cell_length**        | The lattice vectors of periodic unit cell      | *9 \* [ 0.0 ]*     |
| *(double, list)*       |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **sk_path**            | Path for slater-koster files                   | *'./'*             |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **install_path**       | Path for DFTB+ install directory               | *'./'*             |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **mpi**                | Use MPI parallelization                        | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **mpi_path**           | Path for MPI binary                            | *'./'*             |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **nthreads**           | Number of threads in the calculations          | *1*                |
| *(integer)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **version**            | Version of DFTB+ program                       | *'20.1'*           |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **scc** *(boolean)* - Default: *True*

  Includes self-consistent charge (SCC) scheme. If this sets to *False*, then the calculation will change to non-SCC DFTB.

\

- **scc_tol** *(double)* - Default: *1E-6*

  Stopping criteria for the SCC iterations.

\

- **scc_max_iter** *(integer)* - Default: *100*

  Maximum number of SCC iterations.

\

- **ocdftb** *(boolean)* - Default: *False*

  Includes onsite-correction (OC) to SCC term in the DFTB method.

\

- **lcdftb** *(boolean)* - Default: *False*

  Includes long-range corrected (LC) functional in the DFTB method.

\

- **lc_method** *(string)* - Default: *'MatrixBased'*

  Detailed algorithms used in LC-DFTB. These arguments are same with the original arguments used in DFTB+.

  + *'Thresholded'*: Screening according to estimated magnitude of terms.
  + *'NeighbourBased'*: Uses a purely neighbour-list based algorithm.
  + *'MatrixBased'*: Uses a matrix-matrix multiplication based algorithm.

\

- **sdftb** *(boolean)* - Default: *False*

  Includes spin-polarisation scheme in the DFTB method. The atomic spin constants is given in '`$PYUNIXMD`/src/qm/dftb/dftbpar.py',
  and the information about hydrogen, carbon, nitrogen, and oxygen atoms is currently included. If you
  want to exploit spin-polarisation scheme with other species, then add the corresponding
  spin constants to '`$PYUNIXMD`/src/qm/dftb/dftbpar.py' file in the source code.

\

- **unpaired_elec** *(double)* - Default: *0.0*

  Number of unpaired electrons. For example, put two into **unpaired_elec** for calculation of triplet ground state.

\

- **guess** *(string)* - Default: *'h0'*

  Initial guess method for the SCC scheme.

  + *'h0'*: Initial charges of SCC term are set to zero for every time step.
  + *'read'*: Reads "charges.bin" file generated from previous step.
    If **guess_file** exists, then "charges.bin" file is used as initial guess at t = 0.0 s.

\

- **guess_file** *(string)* - Default: *'./charges.bin'*

  Initial guess file for charges. It is vaild when **guess** is *'read'* option.

\

- **elec_temp** *(double)* - Default: *0.0*

  Electronic temperature for Fermi-Dirac scheme. The unit is Kelvin.

\

- **mixer** *(string)* - Default: *'Broyden'*

  Mixing method for charges used in DFTB. These arguments are same with the original arguments in used in DFTB+.
  The detailed parameters used in each mixer are set to default values of the DFTB+ program.
  If you want to know the detailed process of each mixer, see the manual of the DFTB+ program.

  + *'Broyden'*: Use Broyden mixer.
  + *'Anderson'*: Use Anderson mixer.
  + *'DIIS'*: Use DIIS mixer.
  + *'Simple'*: Use simple mixer.

\

- **ex_symmetry** *(string)* - Default: *'singlet'*

  Symmetry of excited state used in TD-DFTB. These arguments are same with the original arguments in used in DFTB+.
  Currently, *'triplet'* and *'both'* options are not added in our interfacing script.

  + *'singlet'*: Calculate singlet excited state in Casida formalism.

\

- **e_window** *(double)* - Default: *0.0*

  Energy window for TD-DFTB. It increases the efficiency of NACME evaluation. The value indicates
  the energy range from the highest orbital among the related orbitals for excited states.
  This option must be treated carefully.

\

- **k_point** *(integer, list)* - Default: *3 \* [ 1 ]*

  Number of K-point samplings. The list consists of three elements.
  If the default is used for the periodic cell, the :math:`\Gamma`-point sampling is used.

\

- **periodic** *(boolean)* - Default: *False*

  Uses a periodicity in the calculation.

\

- **cell_length** *(double, list)* - Default: *9 \* [ 0.0 ]*

  Cell lattice vectors of the periodic unit cell. The list consists of nine elements,
  which correspond to the :math:`a`, :math:`b`, and :math:`c` vectors, respectively.

\

- **sk_path** *(string)* - Default: *'./'*

  Path for slaker-koster files.

\

- **install_path** *(string)* - Default: *'./'*

  Path for DFTB+ install directory. The `$DFTB` environment variable determines the directory where DFTB+ program is installed.
  Thus, **install_path** must be '`$DFTB`/install/', not '`$DFTB`/install/bin/'.

\

- **mpi** *(boolean)* - Default: *False*

  Use MPI parallelization for large scale calculation.

\

- **mpi_path** *(string)* - Default: *'./'*

  Path for MPI binary.

\

- **nthreads** *(integer)* - Default: *1*

  Number of threads in the calculation.

\

- **version** *(string)* - Default: *'20.1'*

  Version of DFTB+ program. (TD)DFTB method is supported in 19.1 version (or newer).

SSR
"""""""""""""""""""""""""""""""""""""

PyUNIxMD automatically determines single-state REKS as BO interfaces for ground state BOMD.
When we include the excited states, SA-REKS or SSR methods can be exploited and these are
determined from the **state_interactions** argument.

+------------------------+------------------------------------------------+---------------------+
| Keywords               | Work                                           | Default             |
+========================+================================================+=====================+
| **molecule**           | Molecular object                               |                     |
| (:class:`Molecule`)    |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **scc**                | Include self-consistent charge (SCC) scheme    | *True*              |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **scc_tol**            | Stopping criteria for the SCC iterations       | *1E-6*              |
| *(double)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **scc_max_iter**       | Maximum number of SCC iterations               | *1000*              |
| *(integer)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **ocdftb**             | Include onsite correction to SCC term          | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **lcdftb**             | Include long-range corrected functional        | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **lc_method**          | Algorithms for LC-DFTB                         | *'MatrixBased'*     |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **ssr22**              | Use SSR(2,2) calculation?                      | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **ssr44**              | Use SSR(4,4) calculation?                      | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **guess**              | Initial guess method for SCC scheme            | *'h0'*              |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **guess_file**         | Initial guess file for eigenvectors            | *'./eigenvec.bin'*  |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **state_interactions** | Include state-interaction terms to SA-REKS     | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **shift**              | Level shifting value in SCC iterations         | *0.3*               |
| *(double)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **tuning**             | Scaling factor for atomic spin constants       | *None*              |
| *(double, list)*       |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **cpreks_grad_alg**    | Algorithms used in CP-REKS equations           | *'pcg'*             |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **cpreks_grad_tol**    | Tolerance used in the conjugate-gradient based | *1E-8*              |
| *(double)*             | algorithm                                      |                     |
+------------------------+------------------------------------------------+---------------------+
| **save_memory**        | Save memory in cache used in CP-REKS equations | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **embedding**          | Charge-charge embedding options in QM/MM       | *None*              |
| *(string)*             | method                                         |                     |
+------------------------+------------------------------------------------+---------------------+
| **periodic**           | Use periodicity in the calculations            | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **cell_length**        | The lattice vectors of periodic unit cell      | *9 \* [ 0.0 ]*      |
| *(double, list)*       |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **sk_path**            | Path for slater-koster files                   | *'./'*              |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **install_path**       | Path for DFTB+ install directory               | *'./'*              |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **nthreads**           | Number of threads in the calculations          | *1*                 |
| *(integer)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **version**            | Version of DFTB+ program                       | *'20.1'*            |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **scc** *(boolean)* - Default: *True*

  Includes self-consistent charge (SCC) scheme. This is a mandatory argument to use SSR calculation.
  If this sets to *False*, then the calculation will be died.

\

- **scc_tol** *(double)* - Default: *1E-6*

  Stopping criteria for the SCC iterations.

\

- **scc_max_iter** *(integer)* - Default: *1000*

  Maximum number of SCC iterations.

\

- **ocdftb** *(boolean)* - Default: *False*

  Includes onsite-correction (OC) to SCC term in the DFTB method. This is currently experimental feature,
  and not implemented in SSR calculation.

\

- **lcdftb** *(boolean)* - Default: *False*

  Includes long-range corrected (LC) functional in the DFTB method. To deal with the excited states properly,
  it is recommended to use LC-DFTB method for SSR calculation.

\

- **lc_method** *(string)* - Default: *'MatrixBased'*

  Detailed algorithms used in LC-DFTB. These arguments are same with the original arguments used in DFTB+.

  + *'Thresholded'*: Screening according to estimated magnitude of terms.
  + *'NeighbourBased'*: Uses a purely neighbour-list based algorithm.
  + *'MatrixBased'*: Uses a matrix-matrix multiplication based algorithm.

\

- **ssr22** *(boolean)* - Default: *False*

  Uses SSR(2,2) calculation in the context of DFTB method. When this sets to *True*, detailed type of the REKS calculation is
  automatically determined from the number of states and **state_interactions** argument. If the number of states is one,
  the single-state REKS calculation is carried out. When the number of states is larger than one,
  the SA-REKS or SI-SA-REKS calculation is executed according to the **state_interactions** argument.

\

- **ssr44** *(boolean)* - Default: *False*

  Uses SSR(4,4) calculation in the context of DFTB method. When this sets to *True*, detailed type of the REKS calculation is
  automatically determined from the number of states and **state_interactions** argument. If the number of states is one,
  the single-state REKS calculation is carried out. When the number of states is larger than one,
  the SA-REKS or SI-SA-REKS calculation is executed according to the **state_interactions** argument.
  This is currently experimental feature and not implemented.

\

- **guess** *(string)* - Default: *'h0'*

  Initial guess method for the SCC scheme. The *'read'* option with DFTB/SSR method is supported in 20.2 version (or newer).

  + *'h0'*: Initial orbitals are generated from the diagonalization of non-SCC Hamiltonian.
  + *'read'*: Reads "eigenvec.bin" file generated from previous step.
    If **guess_file** exists, then "eigenvec.bin" file is used as initial guess at t = 0.0 s.

\

- **guess_file** *(string)* - Default: *'./eigenvec.bin'*

  Initial guess file for eigenvectors. It is vaild when **guess** is *'read'* option.

\

- **state_interactions** *(boolean)* - Default: *False*

  Includes state-interaction terms to SA-REKS calculation. If this sets to *True*, the SI-SA-REKS states are calculated.
  Otherwise, the SA-REKS states are obtained. It is valid when the number of states is larger
  than one. In general, it generates more reliable adiabatic states.

\

- **shift** *(double)* - Default: *0.3*

  Level shifting value used in SCC iterations. It can be helpful to increase **shift** when
  it is hard to converge the SCC iterations.

\

- **tuning** *(double, list)* - Default: *None*

  Scaling factor for atomic spin constants. It must be used carefully.
  The list consists of the number of atomic species.
  For example, if you want to calculate an ethylene molecule with scaling of two which includes carbon and hydrogen atom,
  then you can put *[2.0, 2.0]* into **tuning** argument.

\

- **cpreks_grad_alg** *(string)* - Default: *'pcg'*

  Algorithms used in CP-REKS equations.

  + *'pcg'*: Uses a preconditioned conjugate-gradient based algorithm. It is generally faster than other algorithms.
  + *'cg'*: Uses a conjugate-gradient based algorithm. It is slower than *'pcg'*, but it can be helpful for systems including a high symmetry.
  + *'direct'*: Uses a direct matrix-inversion multiplication algorithm.

\

- **cpreks_grad_tol** *(double)* - Default: *1E-8*

  Tolerance used in the conjugate-gradient based algorithm for solving the CP-REKS equations.
  This is not used when **cpreks_grad_alg** is *'direct'* option.

\

- **save_memory** *(boolean)* - Default: *False*

  Saves memory in cache used in CP-REKS equations. If this sets to *True*, some variables
  which needs large memory allocation are save in the memory. In general, this becomes faster option.
  If this sets to *False*, do not save in cache. This option is recommended for large systems.

\

- **embedding** *(string)* - Default: *None*

  Charge-charge embedding options used in QM/MM method. It is recommended option for the environments showing high polarity.
  The **embedding** of the QM object must be same with the **embedding** defined in the MM object.

  + *None*: Do not use charge-charge embedding in QM/MM method.
  + *'mechanical'*: Uses a mechanical charge-charge embedding option.
    The interactions are treated as the energies between MM point charges.
  + *'electrostatic'*: Uses a electrostatic charge-charge embedding option.
    Point charges as one-electron terms are included in the Hamiltonian.

\

- **periodic** *(boolean)* - Default: *False*

  Uses a periodicity in the calculation. Only :math:`\Gamma`-point sampling is supported with DFTB/SSR method when the periodicity is considered.

\

- **cell_length** *(double, list)* - Default: *9 \* [ 0.0 ]*

  Cell lattice vectors of the periodic unit cell. The list consists of nine elements,
  which correspond to the :math:`a`, :math:`b`, and :math:`c` vectors, respectively.

\

- **sk_path** *(string)* - Default: *'./'*

  Path for slaker-koster files.

\

- **install_path** *(string)* - Default: *'./'*

  Path for DFTB+ install directory. The `$DFTB` environment variable determines the directory where DFTB+ program is installed.
  Thus, **install_path** must be '`$DFTB`/install/', not '`$DFTB`/install/bin/'.

\

- **nthreads** *(integer)* - Default: *1*

  Number of threads in the calculation.

\

- **version** *(string)* - Default: *'20.1'*

  Version of DFTB+ program. DFTB/SSR method is supported in 20.1 version (or newer).

