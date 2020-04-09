
DFTB+ :cite:`Aradi2007,Hourahine2020` is a fast and efficient versatile quantum mechanical simulation software package.
Using DFTB+ you can carry out quantum mechanical simulations similar to density functional
theory but in an approximate way, typically gaining around two orders of magnitude in
speed. (TD)DFTB and SSR methods are interfaced with current version of UNI-xMD.

- (TD)DFTB is time-dependent density-functional tight-binding method. DFTB+ supports only
  analytical gradients, not nonadiabatic couplings. Thus, it can be used for only adiabatic dynamics.

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
| (TD)DFTB          | o    | x  | x  | x   |
+-------------------+------+----+----+-----+
| single-state REKS | o    | x  | x  | x   |
+-------------------+------+----+----+-----+
| SA-REKS           | o    | x  | x  | x   |
+-------------------+------+----+----+-----+
| SI-SA-REKS (SSR)  | o    | o  | o  | o   |
+-------------------+------+----+----+-----+

(TD)DFTB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with 19.1 version of DFTB+ program.
   Here, you should refer to manual of DFTB+ program if you want to see detailed
   lists for **e_temp**, **mixer**, **ex_symmetry** variables.

+-------------------+------------------------------------------------+---------+
| Keywords          | Work                                           | Default |
+===================+================================================+=========+
| scc               | include SCC scheme                             | True    |
+-------------------+------------------------------------------------+---------+
| scc_tol           | energy convergence for SCC iterations          | 1E-6    |
+-------------------+------------------------------------------------+---------+
| max_scc_iter      | maximum number of SCC iterations               | 100     |
+-------------------+------------------------------------------------+---------+
| sdftb             | include spin-polarisation scheme               | False   |
+-------------------+------------------------------------------------+---------+
| unpaired_e        | number of unpaired electrons                   | 0.0     |
+-------------------+------------------------------------------------+---------+
| e_temp            | electronic temperature for Fermi-Dirac scheme  | 0.0     |
+-------------------+------------------------------------------------+---------+
| mixer             | charge mixing method used in SCC-DFTB          | Broyden |
+-------------------+------------------------------------------------+---------+
| ex_symmetry       | symmetry (singlet or triplet) in TD-DFTB       | S       |
+-------------------+------------------------------------------------+---------+
| sk_path           | path for slater-koster files                   | ./      |
+-------------------+------------------------------------------------+---------+
| periodic          | use periodicity in the calculations            | False   |
+-------------------+------------------------------------------------+---------+
| a(b, c)_axis      | the length of cell lattice                     | 0.0     |
+-------------------+------------------------------------------------+---------+
| qm_path           | path for QM binary                             | ./      |
+-------------------+------------------------------------------------+---------+
| nthreads          | number of threads in the calculations          | 1       |
+-------------------+------------------------------------------------+---------+
| mpi               | use MPI parallelization                        | False   |
+-------------------+------------------------------------------------+---------+
| mpi_path          | path for MPI binary                            | ./      |
+-------------------+------------------------------------------------+---------+
| version           | version of DFTB+ program                       | 19.1    |
+-------------------+------------------------------------------------+---------+

SSR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

UNI-xMD automatically determines single-state REKS as BO interfaces for ground state BOMD.
When we include the excited states, SA-REKS or SSR methods can be exploited and these are
determined from the **use_ssr_state** argument.

.. note:: Our interface script is generated with 19.1 version of DFTB+ program.
   Here, you should refer to manual of DFTB+ program if you want to see detailed
   lists for **lc_method**, **state_l**, **guess**, **grad_level**, **mem_level** variables.

+-------------------+------------------------------------------------+---------+
| Keywords          | Work                                           | Default |
+===================+================================================+=========+
| scc               | include SCC scheme                             | True    |
+-------------------+------------------------------------------------+---------+
| scc_tol           | energy convergence for REKS SCC iterations     | 1E-6    |
+-------------------+------------------------------------------------+---------+
| max_scc_iter      | maximum number of REKS SCC iterations          | 1000    |
+-------------------+------------------------------------------------+---------+
| sdftb             | include spin-polarisation parameters           | True    |
+-------------------+------------------------------------------------+---------+
| lcdftb            | include long-range corrected functional        | True    |
+-------------------+------------------------------------------------+---------+
| lc_method         | algorithms for LC-DFTB                         | NB      |
+-------------------+------------------------------------------------+---------+
| ocdftb            | include onsite correction (test option)        | False   |
+-------------------+------------------------------------------------+---------+
| ssr22             | use REKS(2,2) calculation?                     | True    |
+-------------------+------------------------------------------------+---------+
| use_ssr_state     | calculate SSR state, if not, treat SA-REKS     | 1       |
+-------------------+------------------------------------------------+---------+
| state_l           | set L-th microstate as taget state             | 0       |
+-------------------+------------------------------------------------+---------+
| guess             | initial guess setting for eigenvectors         | 1       |
+-------------------+------------------------------------------------+---------+
| shift             | level shifting value in REKS SCF iterations    | 0.3     |
+-------------------+------------------------------------------------+---------+
| tuning            | scaling factor for atomic spin constants       | 1.0     |
+-------------------+------------------------------------------------+---------+
| grad_level        | algorithms to calculate gradients              | 1       |
+-------------------+------------------------------------------------+---------+
| grad_tol          | gradient tolerance for CP-REKS equations       | 1E-8    |
+-------------------+------------------------------------------------+---------+
| mem_level         | memory allocation setting, 2 is recommended    | 2       |
+-------------------+------------------------------------------------+---------+
| sk_path           | path for slater-koster files                   | ./      |
+-------------------+------------------------------------------------+---------+
| periodic          | use periodicity in the calculations            | False   |
+-------------------+------------------------------------------------+---------+
| a(b, c)_axis      | the length of cell lattice                     | 0.0     |
+-------------------+------------------------------------------------+---------+
| qm_path           | path for QM binary                             | ./      |
+-------------------+------------------------------------------------+---------+
| script_path       | path for DFTB+ python script (dptools)         | ./      |
+-------------------+------------------------------------------------+---------+
| nthreads          | number of threads in the calculations          | 1       |
+-------------------+------------------------------------------------+---------+
| version           | version of DFTB+ program                       | 19.1    |
+-------------------+------------------------------------------------+---------+

