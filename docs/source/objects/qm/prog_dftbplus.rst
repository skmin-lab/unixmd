
DFTB+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DFTB+ :cite:`Hourahine2020` is a fast and efficient versatile quantum mechanical simulation software package.
Using DFTB+ you can carry out quantum mechanical simulations similar to density functional
theory but in an approximate way, typically gaining around two orders of magnitude in
speed. (TD)DFTB and SSR methods are interfaced with the current version of PyUNIxMD.

- (TD)DFTB is time-dependent density-functional tight-binding method. DFTB+ supports only
  analytical gradients, not nonadiabatic couplings. Instead, nonadiabatic coupling matrix
  elements (NACME) is calculated by using our wavefunction overlap :cite:`Ryabinkin2015` routines. 
  Thus, it can be used for adiabatic dynamics and surface hopping dynamics.
  For Ehrenfest dynamics, it shows no energy conservation since it cannot calculate
  nonadiabatic couplings directly which is used for nuclear propagation.

- In general, spin-restricted ensemble-referenced Kohn-Sham (REKS) method can be classified
  as single-state REKS, state-averaged REKS (SA-REKS) and state-interaction SA-REKS (SSR).
  In single-state REKS, only ground state is calculated and it can treat the multireference
  character. SA-REKS and SSR can calculate excited states as well as a ground state. The
  difference is that the state-interaction term is considered in SSR so that more accurate
  states can be generated. SSR can provide nonadiabatic couplings so it can be used for
  surface hopping or Ehrenfest dynamics.

+-------------------+------+--------+----+-----+
|                   | BOMD | SH(XF) | Eh | nac |
+===================+======+========+====+=====+
| (TD)DFTB          | o    | o      | x  | x   |
+-------------------+------+--------+----+-----+
| single-state REKS | o    | x      | x  | x   |
+-------------------+------+--------+----+-----+
| SA-REKS           | o    | x      | x  | x   |
+-------------------+------+--------+----+-----+
| SI-SA-REKS (SSR)  | o    | o      | o  | o   |
+-------------------+------+--------+----+-----+

.. note:: To use the DFTB+ interface, the information about the highest maximum angular momentum for each atom type is
   needed. In the current version of PyUNIxMD, the values for the maximum angular momenta are included
   in **max_l** dictionary variable in '`$PYUNIXMDHOME`/src/qm/dftb/dftbpar.py' file.
   You can add or modify the **max_l** variable, see the following example.

.. code-block:: python

   from qm.dftbplus.dftbpar import max_l

   max_l["Si"] = "p" # add value of new Si atom
   max_l["C"] = "s" # modify value of already existing C atom

(TD)DFTB
"""""""""""""""""""""""""""""""""""""

+------------------------+------------------------------------------------+--------------------+
| Parameters             | Work                                           | Default            |
+========================+================================================+====================+
| **molecule**           | Molecule object                                |                    |  
| (:class:`Molecule`)    |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **l_scc**              | Include self-consistent charge (SCC) scheme    | *True*             |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **scc_tol**            | Stopping criteria for the SCC iterations       | *1E-6*             |
| *(double)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **scc_max_iter**       | Maximum number of SCC iterations               | *100*              |
| *(integer)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **l_onsite**           | Include onsite correction to SCC term          | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **l_range_sep**        | Include long-range corrected functional        | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **lc_method**          | Algorithms for LC-DFTB                         | *'MatrixBased'*    |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **l_spin_pol**         | Include spin-polarisation scheme               | *False*            |
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
| **l_periodic**         | Use periodicity in the calculations            | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **cell_length**        | The lattice vectors of periodic unit cell      | *9 \* [ 0.0 ]*     |
| *(double, list)*       |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **sk_path**            | Path for Slater-Koster files                   | *'./'*             |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **install_path**       | Path for DFTB+ install directory               | *'./'*             |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **l_mpi**              | Use MPI parallelization                        | *False*            |
| *(boolean)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **mpi_path**           | Path for MPI binary                            | *'./'*             |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **nthreads**           | Number of threads in the calculations          | *1*                |
| *(integer)*            |                                                |                    |
+------------------------+------------------------------------------------+--------------------+
| **version**            | Version of DFTB+                               | *'20.1'*           |
| *(string)*             |                                                |                    |
+------------------------+------------------------------------------------+--------------------+

Example input for (TD)DFTB
''''''''''''''''''''''''''''''''''''

.. code-block:: python
   :linenos:

   from molecule import Molecule
   import qm, mqc

   geom = '''
   3
   example
   O  1.14  3.77  0.00  0.00  0.00  0.00
   H  2.11  3.77  0.00  0.00  0.00  0.00
   H  0.81  4.45  0.60  0.00  0.00  0.00
   '''

   mol = Molecule(geometry=geom, ndim=3, nstates=2, unit_pos='angs')

   qm = qm.dftbplus.DFTB(molecule=mol, l_scc=True, unpaired_elec=0, guess='h0', \
       ex_symmetry='singlet', sk_path='./', \
       install_path='/opt/dftbplus-20.1/install-openmp/')

   md = mqc.SHXF(molecule=mol, nsteps=100, nesteps=20, dt=0.5, unit_dt='au', \
       sigma=0.1, istate=1, hop_rescale='energy', hop_reject='keep', elec_object='density')
 
   md.run(qm=qm)

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **l_scc** *(boolean)* - Default: *True*

  When **l_scc** is set to *True*, the self-consistent charge (SCC) scheme is included in DFTB.
  If **l_scc** is *False*, then the calculation will change to the non-SCC DFTB.

\

- **scc_tol** *(double)* - Default: *1E-6*

  The SCC cycles are considered converged when the charge error is less than **scc_tol**.
  It is valid when **l_scc** is *True*.

\

- **scc_max_iter** *(integer)* - Default: *100*

  This parameter determines the maximum number of the SCC iterations.

\

- **l_onsite** *(boolean)* - Default: *False*

  When **l_onsite** is set to *True*, onsite-correction (OC) scheme is added to SCC-DFTB.

\

- **l_range_sep** *(boolean)* - Default: *False*

  When **l_range_sep** is set to *True*, long-range corrected (LC) functional is added to SCC-DFTB.
  In this case, the corresponding Slater-Koster files must be used. Check the **sk_path** carefully.

\

- **lc_method** *(string)* - Default: *'MatrixBased'*

  This parameter specifies the detailed algorithms used in LC-DFTB.
  The available options of the parameter are the same as the original ones of DFTB+.

  + *'Thresholded'*: Screening according to estimated magnitude of terms.
  + *'NeighbourBased'*: Uses a purely neighbour-list based algorithm.
  + *'MatrixBased'*: Uses a matrix-matrix multiplication based algorithm.

\

- **l_spin_pol** *(boolean)* - Default: *False*

  When **l_spin_pol** is set to *True*, the spin-polarisation scheme is added to SCC-DFTB.
  The atomic spin constants are given in '`$PYUNIXMD`/src/qm/dftb/dftbpar.py',
  and the values about hydrogen, carbon, nitrogen, and oxygen atoms are currently included.
  If you want to exploit spin-polarization scheme with other atomic species, then add the
  corresponding spin constants to '`$PYUNIXMD`/src/qm/dftb/dftbpar.py' file in the source code.

\

- **unpaired_elec** *(double)* - Default: *0.0*

  This parameter specifies the number of unpaired electrons. For example,
  put *2.0* into **unpaired_elec** for calculation of triplet ground state.

\

- **guess** *(string)* - Default: *'h0'*

  This parameter determines the initial guess method for the SCC-DFTB calculations.

  + *'h0'*: Initial guess charges for SCC-DFTB calculations are set to zeros.
  + *'read'*: Initial guess charges are read from the 'charges.bin' file which contains the charges calculated at the previous time step.

\

- **guess_file** *(string)* - Default: *'./charges.bin'*

  The **guess_file** determines the name of the file containing orbitals for the initial guess of orbitals for the SCC-DFTB calculation at the first MD step.
  This parameter is effective only if **guess** = *'read'*.
  If the file does not exist, the *'h0'* option is applied for the initial guess for the SCC-DFTB calculation at the first MD step.

\

- **elec_temp** *(double)* - Default: *0.0*

  This parameter determines the electronic temperature in the Fermi-Dirac scheme. The unit is K.

\

- **mixer** *(string)* - Default: *'Broyden'*

  This parameter specifies the mixing method for charges used in SCC-DFTB.
  The available options of the parameter are the same as the original ones of DFTB+.
  The detailed parameters used in each mixer are set to default values of DFTB+.
  If you want to know the detailed process of each mixer, see the manual of DFTB+.
  Following four mixers can be used in the current interface; {*'Broyden'*, *'Anderson'*, *'DIIS'*, *'Simple'*}

\

- **ex_symmetry** *(string)* - Default: *'singlet'*

  This parameter specifies the symmetry of excited states used in TD-DFTB.
  The available options of the parameter are the same as the original ones of DFTB+.
  Currently, *'triplet'* and *'both'* options are not added in our interface.

  + *'singlet'*: Calculate singlet excited states in Casida formalism.

\

- **e_window** *(double)* - Default: *0.0*

  This parameter determines the energy window for TD-DFTB. It increases the efficiency
  of NACME evaluation. **e_window** indicates the energy range above the last transition at the
  highest excitation to be included in the excited state calculation. This option must be treated carefully.

\

- **k_point** *(integer, list)* - Default: *3 \* [ 1 ]*

  This parameter specifies the number of K-point samplings. The list consists of three elements.
  If the default is used for the periodic cell, the :math:`\Gamma`-point sampling is used.

\

- **l_periodic** *(boolean)* - Default: *False*

  When **l_periodic** is set to *True*, periodicity is considered in the calculation.

\

- **cell_length** *(double, list)* - Default: *9 \* [ 0.0 ]*

  This parameter specifies the cell lattice vectors of the periodic cell. The list consists of nine elements,
  which correspond to the :math:`a`, :math:`b`, and :math:`c` vectors, respectively.

\

- **sk_path** *(string)* - Default: *'./'*

  This parameter determines the path for Slaker-Koster files.

\

- **install_path** *(string)* - Default: *'./'*

  This parameter determines the path for DFTB+ install directory. The `$DFTB` environment
  variable determines the directory where DFTB+ is installed
  (For example, `$DFTB` is '/my_disk/my_name/dftbplus-**version**/').
  Thus, **install_path** must be *'`$DFTB`/install/'*, not *'`$DFTB`/install/bin/'*.

\

- **mpi** *(boolean)* - Default: *False*

  When **mpi** is set to *True*, MPI parallelization is used for large scale calculations.
  This option can be used when only ground state is included in the calculations.

\

- **mpi_path** *(string)* - Default: *'./'*

  This parameter determines the path for MPI binaries.

\

- **nthreads** *(integer)* - Default: *1*

  This parameter specifies the number of threads in the calculation.

\

- **version** *(string)* - Default: *'20.1'*

  This parameter determines the version of DFTB+.
  PyUNIxMD is currently based on version 19.1 and 20.1 of DFTB+.

SSR
"""""""""""""""""""""""""""""""""""""

PyUNIxMD automatically determines the single-state REKS as BO interfaces for ground state BOMD.
When we include excited states, the SA-REKS, SSR methods can be exploited and these are
determined from the **l_state_interactions** parameter.

.. note:: In the case of the SSR method, the calculation is possible only when the number
   of states (``molecule.nst``) is smaller than 4 due to the limited active space.
   If you want to treat more excited states, then increase the active space.

+--------------------------+------------------------------------------------+---------------------+
| Parameters               | Work                                           | Default             |
+==========================+================================================+=====================+
| **molecule**             | Molecule object                                |                     |
| (:class:`Molecule`)      |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **l_scc**                | Include self-consistent charge (SCC) scheme    | *True*              |
| *(boolean)*              |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **scc_tol**              | Stopping criteria for the SCC iterations       | *1E-6*              |
| *(double)*               |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **scc_max_iter**         | Maximum number of SCC iterations               | *1000*              |
| *(integer)*              |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **l_onsite**             | Include onsite correction to SCC term          | *False*             |
| *(boolean)*              |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **l_range_sep**          | Include long-range corrected functional        | *False*             |
| *(boolean)*              |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **lc_method**            | Algorithms for LC-DFTB                         | *'MatrixBased'*     |
| *(string)*               |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **active_space**         | Active space for DFTB/SSR calculation          | *2*                 |
| *(integer)*              |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **guess**                | Initial guess method for SCC scheme            | *'h0'*              |
| *(string)*               |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **guess_file**           | Initial guess file for eigenvectors            | *'./eigenvec.bin'*  |
| *(string)*               |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **l_state_interactions** | Include state-interaction terms to SA-REKS     | *False*             |
| *(boolean)*              |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **shift**                | Level shifting value in SCC iterations         | *0.3*               |
| *(double)*               |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **tuning**               | Scaling factor for atomic spin constants       | *None*              |
| *(double, list)*         |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **cpreks_grad_alg**      | Algorithms used in CP-REKS equations           | *'pcg'*             |
| *(string)*               |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **cpreks_grad_tol**      | Tolerance used in the conjugate-gradient based | *1E-8*              |
| *(double)*               | algorithm                                      |                     |
+--------------------------+------------------------------------------------+---------------------+
| **l_save_memory**        | Save memory in cache used in CP-REKS equations | *False*             |
| *(boolean)*              |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **embedding**            | Charge-charge embedding options in QM/MM       | *None*              |
| *(string)*               | method                                         |                     |
+--------------------------+------------------------------------------------+---------------------+
| **l_periodic**           | Use periodicity in the calculations            | *False*             |
| *(boolean)*              |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **cell_length**          | The lattice vectors of periodic unit cell      | *9 \* [ 0.0 ]*      |
| *(double, list)*         |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **sk_path**              | Path for Slater-Koster files                   | *'./'*              |
| *(string)*               |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **install_path**         | Path for DFTB+ install directory               | *'./'*              |
| *(string)*               |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **nthreads**             | Number of threads in the calculations          | *1*                 |
| *(integer)*              |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+
| **version**              | Version of DFTB+                               | *'20.1'*            |
| *(string)*               |                                                |                     |
+--------------------------+------------------------------------------------+---------------------+

Example input for SSR
''''''''''''''''''''''''''''''''''''

.. code-block:: python
   :linenos:

   from molecule import Molecule
   import qm, mqc

   geom = '''
   3
   example
   O  1.14  3.77  0.00  0.00  0.00  0.00
   H  2.11  3.77  0.00  0.00  0.00  0.00
   H  0.81  4.45  0.60  0.00  0.00  0.00
   '''

   mol = Molecule(geometry=geom, ndim=3, nstates=2, unit_pos='angs')

   qm = qm.dftbplus.SSR(molecule=mol, l_scc=True, active_space=2, guess='h0', \
       l_state_interactions=True, shift=0.3, embedding=None, sk_path='./', \
       install_path='/opt/dftbplus-20.1/install-openmp/')

   md = mqc.SHXF(molecule=mol, nsteps=100, nesteps=20, dt=0.5, unit_dt='au', \
       sigma=0.1, istate=1, hop_rescale='energy', elec_object='density')

   md.run(qm=qm)

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **l_scc** *(boolean)* - Default: *True*

  When **l_scc** is set to *True*, the self-consistent charge (SCC) scheme is included in DFTB/SSR.
  If **l_scc** is *False*, then the calculation will be halted since the SCC scheme is a mandatory option.

\

- **scc_tol** *(double)* - Default: *1E-6*

  The SCC cycles are considered converged when the charge error is less than **scc_tol**.
  It is valid when **l_scc** is *True*.

\

- **scc_max_iter** *(integer)* - Default: *1000*

  This parameter determines the maximum number of the SCC iterations.

\

- **l_onsite** *(boolean)* - Default: *False*

  When **l_onsite** is set to *True*, onsite-correction (OC) scheme is added to DFTB/SSR.
  It is currently experimental feature, and not implemented in the SSR calculation.

\

- **l_range_sep** *(boolean)* - Default: *False*

  When **l_range_sep** is set to *True*, long-range corrected (LC) functional is added to DFTB/SSR.
  To deal with the excited states properly, it is recommended to use LC funtionals for the DFTB/SSR calculations.
  In this case, the corresponding Slater-Koster files must be used. Check the **sk_path** carefully.

\

- **lc_method** *(string)* - Default: *'MatrixBased'*

  This parameter specifies the detailed algorithms used in LC-DFTB.
  The available options of the parameter are the same as the original ones of DFTB+.

  + *'Thresholded'*: Screening according to estimated magnitude of terms.
  + *'NeighbourBased'*: Uses a purely neighbour-list based algorithm.
  + *'MatrixBased'*: Uses a matrix-matrix multiplication based algorithm.

\

- **active_space** *(integer)* - Default: *2*

  This parameter specifies the active space for DFTB/SSR calculation. Detailed types of the REKS calculation are
  automatically determined by ``molecule.nst`` and **l_state_interactions** parameters. If ``molecule.nst`` is *1*,
  the single-state REKS calculation is carried out. When ``molecule.nst`` is larger than *1*,
  the SA-REKS or the SI-SA-REKS calculation is executed according to the **l_state_interactions** parameter.
  Currently, only (2,2) space is available for DFTB/SSR calculation.

  + *2*: The numbers of electrons and orbitals are 2 and 2, respectively.

\

- **guess** *(string)* - Default: *'h0'*

  This parameter determines the initial guess method for the DFTB/SSR method.
  The *'read'* option with the DFTB/SSR method is supported in version 20.2 (or newer).

  + *'h0'*: Initial guess orbitals for the DFTB/SSR method are generated from the diagonalization of the non-SCC Hamiltonian.
  + *'read'*: Initial guess orbitals are read from the 'eigenvec.bin' file which contains the orbitals calculated at the previous time step.

\

- **guess_file** *(string)* - Default: *'./eigenvec.bin'*

  The **guess_file** determines the name of the file containing orbitals for the initial guess of orbitals for the DFTB/SSR calculation at the first MD step.
  This parameter is effective only if **guess** = *'read'*.
  If the file does not exist, *'h0'* option is applied for the initial guess for the DFTB/SSR calculation at the first MD step.

\

- **l_state_interactions** *(boolean)* - Default: *False*

  When **l_state_interactions** is set to *True*, state-interaction terms are included so that the SI-SA-REKS states are generated.
  Otherwise, the SA-REKS states are obtained. It is valid when ``molecule.nst`` is larger
  than *1*. In general, it generates more reliable adiabatic states.

\

- **shift** *(double)* - Default: *0.3*

  This parameter specifies the level shifting value used in the SCC iterations. It can be helpful to increase **shift** when
  it is hard to converge the SCC iterations.

\

- **tuning** *(double, list)* - Default: *None*

  This parameter specifies the scaling factor for atomic spin constants. It must be used carefully.
  The list consists of the number of atomic species.
  For example, if you want to calculate an ethylene molecule with scaling factor of two which includes carbon and hydrogen atom,
  then you can put *[2.0, 2.0]* into **tuning** parameter.

\

- **cpreks_grad_alg** *(string)* - Default: *'pcg'*

  This parameter specifies the detailed algorithms used to solve the CP-REKS equations.

  + *'pcg'*: Uses a preconditioned conjugate-gradient based algorithm. It is generally faster than other algorithms.
  + *'cg'*: Uses a conjugate-gradient based algorithm. It is slower than *'pcg'*, but it can be helpful for systems including a high symmetry.
  + *'direct'*: Uses a direct matrix-inversion multiplication algorithm. It requires large memory allocation.

\

- **cpreks_grad_tol** *(double)* - Default: *1E-8*

  This parameter determines the tolerance used in the conjugate-gradient based algorithm for solving the CP-REKS equations.
  This is not used when **cpreks_grad_alg** is *'direct'*.

\

- **l_save_memory** *(boolean)* - Default: *False*

  This parameter controls whether to save memory used in the CP-REKS equations in cache or not.
  If **l_save_memory** sets to *True*, some variables which needs large memory allocation are saved in the memory.
  In general, this becomes a faster option. If **l_save_memory** sets to *False*, not saved in the cache.
  This option is recommended for large systems.

\

- **embedding** *(string)* - Default: *None*

  This parameter specifies the charge-charge embedding option used in the QM/MM method.
  It is recommended option for the environments showing high polarity.
  The **embedding** of the QM object must be same with the **embedding** defined in the MM object.
  If this parameter is *None*, the charge-charge embedding is not included in the QM/MM calculation.

  + *'mechanical'*: Uses a mechanical charge-charge embedding option.
    The interactions are treated as the energies between MM point charges.
  + *'electrostatic'*: Uses a electrostatic charge-charge embedding option.
    Point charges as one-electron terms are included in the Hamiltonian.

\

- **l_periodic** *(boolean)* - Default: *False*

  When **l_periodic** is set to *True*, periodicity is considered in the calculation.
  Only :math:`\Gamma`-point sampling is supported with the DFTB/SSR method when the periodicity is considered.

\

- **cell_length** *(double, list)* - Default: *9 \* [ 0.0 ]*

  This parameter specifies the cell lattice vectors of the periodic cell. The list consists of nine elements,
  which correspond to the :math:`a`, :math:`b`, and :math:`c` vectors, respectively.

\

- **sk_path** *(string)* - Default: *'./'*

  This parameter determines the path for Slaker-Koster files.

\

- **install_path** *(string)* - Default: *'./'*

  This parameter determines the path for DFTB+ install directory. The `$DFTB` environment
  variable determines the directory where DFTB+ is installed
  (For example, `$DFTB` is '/my_disk/my_name/dftbplus-**version**/').
  Thus, **install_path** must be *'`$DFTB`/install/'*, not *'`$DFTB`/install/bin/'*.

\

- **nthreads** *(integer)* - Default: *1*

  This parameter specifies the number of threads in the calculation.

\

- **version** *(string)* - Default: *'20.1'*

  This parameter determines the version of DFTB+.
  PyUNIxMD is currently based on version 20.1 of DFTB+.

