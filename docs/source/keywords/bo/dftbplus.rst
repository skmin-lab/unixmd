
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

(TD)DFTB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with 19.1/20.1 version of DFTB+ program.
   Here, you should refer to manual of DFTB+ program if you want to see detailed
   lists for ``lc_method``, ``mixer``, ``ex_symmetry`` variables.

.. note:: Currently, ``guess`` variable reads the following two strings.
   One is ``h0``, which uses zedo charges as initial guess of SCC term for every time step.
   The other is ``read``, which reads charges.bin file generated from previous step.
   If ``guess_file`` exists, then charges.bin file is used as initial guess at t = 0.0 s.

.. note:: For ``cell_length`` variable, it reads a list variable consisted of 9 elements,
   which correspond to cell lattice vectors.

+-------------------+------------------------------------------------+---------------------+
| Keywords          | Work                                           | Default             |
+===================+================================================+=====================+
| ``scc``           | include self-consistent charge (SCC) scheme    | ``True``            |
+-------------------+------------------------------------------------+---------------------+
| ``scc_tol``       | energy convergence for SCC iterations          | ``1E-6``            |
+-------------------+------------------------------------------------+---------------------+
| ``scc_max_iter``  | maximum number of SCC iterations               | ``100``             |
+-------------------+------------------------------------------------+---------------------+
| ``ocdftb``        | include onsite correction to SCC term          | ``False``           |
+-------------------+------------------------------------------------+---------------------+
| ``lcdftb``        | include long-range corrected functional        | ``False``           |
+-------------------+------------------------------------------------+---------------------+
| ``lc_method``     | algorithms for LC-DFTB                         | ``MatrixBased``     |
+-------------------+------------------------------------------------+---------------------+
| ``sdftb``         | include spin-polarisation scheme               | ``False``           |
+-------------------+------------------------------------------------+---------------------+
| ``unpaired_elec`` | number of unpaired electrons                   | ``0.0``             |
+-------------------+------------------------------------------------+---------------------+
| ``guess``         | initial guess method for SCC scheme            | ``h0``              |
+-------------------+------------------------------------------------+---------------------+
| ``guess_file``    | initial guess file for charges                 | ``./charges.bin``   |
+-------------------+------------------------------------------------+---------------------+
| ``elec_temp``     | electronic temperature for Fermi-Dirac scheme  | ``0.0``             |
+-------------------+------------------------------------------------+---------------------+
| ``mixer``         | charge mixing method used in DFTB              | ``Broyden``         |
+-------------------+------------------------------------------------+---------------------+
| ``ex_symmetry``   | symmetry of excited state in TD-DFTB           | ``singlet``         |
+-------------------+------------------------------------------------+---------------------+
| ``periodic``      | use periodicity in the calculations            | ``False``           |
+-------------------+------------------------------------------------+---------------------+
| ``cell_length``   | the lattice vectors of periodic unit cell      | ``[ 9 * 0.0 ]``     |
+-------------------+------------------------------------------------+---------------------+
| ``sk_path``       | path for slater-koster files                   | ``./``              |
+-------------------+------------------------------------------------+---------------------+
| ``install_path``  | path for DFTB+ install directory               | ``./``              |
+-------------------+------------------------------------------------+---------------------+
| ``mpi``           | use MPI parallelization                        | ``False``           |
+-------------------+------------------------------------------------+---------------------+
| ``mpi_path``      | path for MPI binary                            | ``./``              |
+-------------------+------------------------------------------------+---------------------+
| ``nthreads``      | number of threads in the calculations          | ``1``               |
+-------------------+------------------------------------------------+---------------------+
| ``version``       | version of DFTB+ program                       | ``20.1``            |
+-------------------+------------------------------------------------+---------------------+

SSR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

UNI-xMD automatically determines single-state REKS as BO interfaces for ground state BOMD.
When we include the excited states, SA-REKS or SSR methods can be exploited and these are
determined from the ``state_interactions`` argument.

.. note:: Our interface script is generated with 20.1 version of DFTB+ program.
   Here, you should refer to manual of DFTB+ program if you want to see detailed
   lists for ``lc_method`` variable.

.. note:: Currently, ``ocdftb`` is not implemented in current version of DFTB+.

.. note:: Currently, ``guess`` variable reads the following two strings.
   One is ``h0``, which uses zedo charges as initial guess of SCC term for every time step.
   The other is ``read``, which reads charges.bin file generated from previous step.
   If ``guess_file`` exists, then charges.bin file is used as initial guess at t = 0.0 s.

.. note:: For ``cell_length`` variable, it reads a list variable consisted of 9 elements,
   which correspond to cell lattice vectors. Similarly ``tuning`` variable reads a list
   with as many as the number of atomic species.

+------------------------+------------------------------------------------+---------------------+
| Keywords               | Work                                           | Default             |
+========================+================================================+=====================+
| ``scc``                | include self-consistent charge (SCC) scheme    | ``True``            |
+------------------------+------------------------------------------------+---------------------+
| ``scc_tol``            | energy convergence for SCC iterations          | ``1E-6``            |
+------------------------+------------------------------------------------+---------------------+
| ``scc_max_iter``       | maximum number of SCC iterations               | ``1000``            |
+------------------------+------------------------------------------------+---------------------+
| ``ocdftb``             | include onsite correction to SCC term          | ``False``           |
+------------------------+------------------------------------------------+---------------------+
| ``lcdftb``             | include long-range corrected functional        | ``False``           |
+------------------------+------------------------------------------------+---------------------+
| ``lc_method``          | algorithms for LC-DFTB                         | ``MatrixBased``     |
+------------------------+------------------------------------------------+---------------------+
| ``ssr22``              | use REKS(2,2) calculation?                     | ``False``           |
+------------------------+------------------------------------------------+---------------------+
| ``ssr44``              | use REKS(4,4) calculation?                     | ``False``           |
+------------------------+------------------------------------------------+---------------------+
| ``guess``              | initial guess method for SCC scheme            | ``h0``              |
+------------------------+------------------------------------------------+---------------------+
| ``guess_file``         | initial guess file for eigenvectors            | ``./eigenvec.bin``  |
+------------------------+------------------------------------------------+---------------------+
| ``state_interactions`` | include state-interaction terms to SA-REKS     | ``False``           |
+------------------------+------------------------------------------------+---------------------+
| ``shift``              | level shifting value in SCC iterations         | ``0.3``             |
+------------------------+------------------------------------------------+---------------------+
| ``tuning``             | scaling factor for atomic spin constants       | ``None``            |
+------------------------+------------------------------------------------+---------------------+
| ``cpreks_grad_alg``    | algorithms used in CP-REKS equations           | ``1``               |
+------------------------+------------------------------------------------+---------------------+
| ``cpreks_grad_tol``    | gradient tolerance for CP-REKS equations       | ``1E-8``            |
+------------------------+------------------------------------------------+---------------------+
| ``save_memory``        | save memory in cache used in CP-REKS equations | ``False``           |
+------------------------+------------------------------------------------+---------------------+
| ``periodic``           | use periodicity in the calculations            | ``False``           |
+------------------------+------------------------------------------------+---------------------+
| ``cell_length``        | the lattice vectors of periodic unit cell      | ``[ 9 * 0.0 ]``     |
+------------------------+------------------------------------------------+---------------------+
| ``sk_path``            | path for slater-koster files                   | ``./``              |
+------------------------+------------------------------------------------+---------------------+
| ``install_path``       | path for DFTB+ install directory               | ``./``              |
+------------------------+------------------------------------------------+---------------------+
| ``nthreads``           | number of threads in the calculations          | ``1``               |
+------------------------+------------------------------------------------+---------------------+
| ``version``            | version of DFTB+ program                       | ``20.1``            |
+------------------------+------------------------------------------------+---------------------+

