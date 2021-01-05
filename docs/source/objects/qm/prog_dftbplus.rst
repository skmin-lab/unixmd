
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

UNI-xMD automatically determines single-state REKS as BO interfaces for ground state BOMD.
When we include the excited states, SA-REKS or SSR methods can be exploited and these are
determined from the **state_interactions** argument.

.. note:: Our interface script is generated with 20.1 version of DFTB+ program.
   Here, you should refer to manual of DFTB+ program if you want to see detailed
   lists for **lc_method** variable.

.. note:: Currently, **ocdftb** is not implemented in current version of DFTB+.

.. note:: Currently, **guess** variable reads the following two strings.
   One is *'h0'*, which uses zedo charges as initial guess of SCC term for every time step.
   The other is **read**, which reads charges.bin file generated from previous step.
   If **guess_file** exists, then charges.bin file is used as initial guess at t = 0.0 s.

.. note:: For **cell_length** variable, it reads a list variable consisted of 9 elements,
   which correspond to cell lattice vectors. Similarly **tuning** variable reads a list
   with as many as the number of atomic species.

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
| *(double)*             |                                                |                     |
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
| **embedding**          | charge embedding options; electrostatic,       | *None*              |
| *(string)*             | mechanical                                     |                     |
+------------------------+------------------------------------------------+---------------------+
| **periodic**           | use periodicity in the calculations            | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **cell_length**        | the lattice vectors of periodic unit cell      | *9 \* [ 0.0 ]*      |
| *(double, list)*       |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **sk_path**            | path for slater-koster files                   | *'./'*              |
| *(string)*             | mechanical                                     |                     |
+------------------------+------------------------------------------------+---------------------+
| **install_path**       | path for DFTB+ install directory               | *'./'*              |
| *(string)*             | mechanical                                     |                     |
+------------------------+------------------------------------------------+---------------------+
| **nthreads**           | number of threads in the calculations          | *1*                 |
| *(integer)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **version**            | version of DFTB+ program                       | *'20.1'*            |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+

