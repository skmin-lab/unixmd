=====================================
Columbus
=====================================

Columbus[ref] is one of famous free software for high-level ab initio quantum calculation. Similar with
other software, it can do various types of fundamental quantum calculation. However, the major
competitiveness of Columbus compare to other software is, it is mainly designed for calculate
multi-reference calculations on electonic ground and excited states. This feature is indeed well suited
for dynamics in UNI-xMD, it is implemented for various types of dynamics. In the current version of UNI-xMD,
only CASSCF method is available.

- CASSCF is complete active space self-consistent field method. It provides analytical gradients as
  well as nonadiabatic couplings, thus it can be used for excited state molecular dynamics.

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| CASSCF | o    | o  | o  | o   |
+--------+------+----+----+-----+

CASSCF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with 7.0 version of Columbus program.
   Here, you should refer to manual of Columbus program if you want to see detailed
   lists for **basis_set** variable.

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| basis_set      | basis set information                          | 6-31g*  |
+----------------+------------------------------------------------+---------+
| memory         | allocatable memory in the calculations         | 500     |
+----------------+------------------------------------------------+---------+
| active_elec    | number of electrons in active space            | 2       |
+----------------+------------------------------------------------+---------+
| active_orb     | number of orbitals in active space             | 2       |
+----------------+------------------------------------------------+---------+
| qm_path        | path for QM binary                             | ./      |
+----------------+------------------------------------------------+---------+
| nthreads       | number of threads in the calculations          | 1       |
+----------------+------------------------------------------------+---------+
| version        | version of Molpro program                      | 7.0     |
+----------------+------------------------------------------------+---------+


=====================================
DFTBplus
=====================================

DFTB+ is a fast and efficient versatile quantum mechanical simulation software package.
Using DFTB+ you can carry out quantum mechanical simulations similar to density functional
theory but in an approximate way, typically gaining around two orders of magnitude in
speed. (TD)DFTB and SSR methods are interfaced with current version of UNI-XMD.

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

UNI-XMD automatically determines single-state REKS as BO interfaces for ground state BOMD.
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
| nthreads          | number of threads in the calculations          | 1       |
+-------------------+------------------------------------------------+---------+
| version           | version of DFTB+ program                       | 19.1    |
+-------------------+------------------------------------------------+---------+

=====================================
Molpro
=====================================

Molpro is a comprehensive system of ab initio programs for advanced molecular electronic structure
calculations, designed and maintained by many authors. It comprises efficient and well parallelized
programs for standard computational chemistry applications, such as DFT or many wave function based
methods. Among them, CASSCF method is interfaced with current version of UNI-XMD.

- CASSCF is complete active space self-consistent field method. Molpro supports analytical gradients as
  well as nonadiabatic couplings, thus it can be used for excited state molecular dynamics.

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| CASSCF | o    | o  | o  | o   |
+--------+------+----+----+-----+

CASSCF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with 2015.1 version of Molpro program.
   Here, you should refer to manual of Molpro program if you want to see detailed
   lists for **basis_set** variable.

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| basis_set      | basis set information                          | sto-3g  |
+----------------+------------------------------------------------+---------+
| memory         | allocatable memory in the calculations         | 500m    |
+----------------+------------------------------------------------+---------+
| max_iter       | maximum number of MCSCF iterations             | 20      |
+----------------+------------------------------------------------+---------+
| scf_en_tol     | energy convergence for MCSCF iterations        | 1E-8    |
+----------------+------------------------------------------------+---------+
| scf_grad_tol   | gradient convergence for MCSCF iterations      | 1E-6    |
+----------------+------------------------------------------------+---------+
| scf_step_tol   | step length convergence for MCSCF iterations   | 1E-2    |
+----------------+------------------------------------------------+---------+
| active_elec    | number of electrons in active space            | 2       |
+----------------+------------------------------------------------+---------+
| active_orb     | number of orbitals in active space             | 2       |
+----------------+------------------------------------------------+---------+
| cpscf_grad_tol | gradient tolerance for CP-MCSCF equations      | 1E-7    |
+----------------+------------------------------------------------+---------+
| qm_path        | path for QM binary                             | ./      |
+----------------+------------------------------------------------+---------+
| nthreads       | number of threads in the calculations          | 1       |
+----------------+------------------------------------------------+---------+
| version        | version of Molpro program                      | 2015.1  |
+----------------+------------------------------------------------+---------+




=====================================
Gaussian09
=====================================

Gaussian09 has been a standard program for electronic structure calculations of molecules.
The only BOMD using the DFT option is available with Gaussian09 in the current version of UNI-xMD,
because it doesn't explicitly provide with nonadiabatic coupling vectors. 
Numerical calculation of the coupling elements using the CI overlap is on progress, which allows the other dynamics options.

- (TD)DFT is used to provide with a potential energy and its gradient for a certain adiabatic state.

+---------+------+----+----+-----+
|         | BOMD | SH | Eh | nac |
+=========+======+====+====+=====+
| (TD)DFT | o    | x  | x  | x   |
+---------+------+----+----+-----+


(TD-)DFT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with Revision A.02 version of Gaussian09 program.
   Please refer to the manual for the detailed lists for **basis_set** and **functional** variable.

+----------------+------------------------------------------------+---------------+
| Keywords       | Work                                           | Default       |
+================+================================================+===============+
| basis_set      | basis set information                          | sto-3g        |
+----------------+------------------------------------------------+---------------+
| memory         | allocatable memory in the calculations         | 500m          |
+----------------+------------------------------------------------+---------------+
| functional     | the level of DFT theory                        | BLYP          |
+----------------+------------------------------------------------+---------------+
| G09_root_path  | path for G09 root                              | /opt/gaussian |
+----------------+------------------------------------------------+---------------+
| nthreads       | number of threads in the calculations          | 1             |
+----------------+------------------------------------------------+---------------+
| version        | version of Gaussian09 program                  | Revision A.02 |
+----------------+------------------------------------------------+---------------+

=====================================
TeraChem
=====================================

TeraChem is general purpose quantum chemistry software designed to run on NVIDIA GPU
architectures under a 64-bit Linux operating system. It includes many functionalities
such as DFT or wave function based methods. Among them, SSR method is interfaced with
current version of UNI-XMD.

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
| single-state REKS | o    | x  | x  | x   |
+-------------------+------+----+----+-----+
| SA-REKS           | o    | x  | x  | x   |
+-------------------+------+----+----+-----+
| SI-SA-REKS (SSR)  | o    | o  | o  | o   |
+-------------------+------+----+----+-----+

SSR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

UNI-XMD automatically determines single-state REKS as BO interfaces for ground state BOMD.
When we include the excited states, SA-REKS or SSR methods can be exploited and these are
determined from the **use_ssr_state** argument.

.. note:: Our interface script is generated with 1.92, 1.93 version of TeraChem program.
   Here, you should refer to manual of TeraChem program if you want to see detailed
   lists for **basis_set**, **functional**, **precision** variables.

+-------------------+------------------------------------------------+---------+
| Keywords          | Work                                           | Default |
+===================+================================================+=========+
| basis_set         | basis set information                          | sto-3g  |
+-------------------+------------------------------------------------+---------+
| functional        | functional in the calculations                 | hf      |
+-------------------+------------------------------------------------+---------+
| precision         | precision in the calculations                  | dynamic |
+-------------------+------------------------------------------------+---------+
| scf_tol           | energy convergence for SCF iterations          | 1E-2    |
+-------------------+------------------------------------------------+---------+
| max_scf_iter      | maximum number of SCF iterations               | 300     |
+-------------------+------------------------------------------------+---------+
| reks22            | use REKS(2,2) calculation?                     | yes     |
+-------------------+------------------------------------------------+---------+
| reks_scf_tol      | energy convergence for REKS SCF iterations     | 1E-6    |
+-------------------+------------------------------------------------+---------+
| reks_max_scf_iter | maximum number of REKS SCF iterations          | 1000    |
+-------------------+------------------------------------------------+---------+
| reks_diis         | DIIS acceleration in REKS SCF iterations       | yes     |
+-------------------+------------------------------------------------+---------+
| shift             | level shifting value in REKS SCF iterations    | 0.3     |
+-------------------+------------------------------------------------+---------+
| use_ssr_state     | calculate SSR state, if not, treat SA-REKS     | 1       |
+-------------------+------------------------------------------------+---------+
| cpreks_max_tol    | gradient tolerance for CP-REKS equations       | 1E-6    |
+-------------------+------------------------------------------------+---------+
| cpreks_max_iter   | maximum number of CP-REKS iterations           | 1000    |
+-------------------+------------------------------------------------+---------+
| qm_path           | path for QM binary                             | ./      |
+-------------------+------------------------------------------------+---------+
| ngpus             | number of GPUs                                 | 1       |
+-------------------+------------------------------------------------+---------+
| gpu_id            | ID of used GPUs                                | 1       |
+-------------------+------------------------------------------------+---------+
| version           | version of TeraChem program                    | 1.92    |
+-------------------+------------------------------------------------+---------+

=====================================
Turbomole
=====================================

Turbomole (TM) is quantum chemical program package, initially developed in the group of Prof. Dr. Reinhart Ahlrichs at the University of Karlsruhe and at the Forschungszentrum Karlsruhe. 
(TD)DFT method is interfaced with current version of UNI-XMD. 

- (TD)DFT provides analytical gradients, thus it can be used born-oppenhiemer molecular dynamics (BOMD).

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| TDDFT  | o    | x  | x  | x   |
+--------+------+----+----+-----+

(TD)DFT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with 6.4 version of TM program.
   Here, you should refer to manual of TM program if you want to see detailed
   lists for **basis_set**, **functional** variable.

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| functional     | xc functional information                      | b-lyp   |
+----------------+------------------------------------------------+---------+
| basis_set      | basis set information                          | sto-3g  |
+----------------+------------------------------------------------+---------+
| memory         | allocatable memory in the calculations         | 500m    |
+----------------+------------------------------------------------+---------+
| max_iter       | maximum number of SCF iterations               | 20      |
+----------------+------------------------------------------------+---------+
| scf_en_tol     | energy convergence for SCF iterations          | 1E-8    |
+----------------+------------------------------------------------+---------+
| qm_path        | path for QM program                            | ./      |
+----------------+------------------------------------------------+---------+
| qm_bin_path    | path for QM binary                             | ./      |
+----------------+------------------------------------------------+---------+
| qm_scripts_path| path for QM scripts                            | ./      |
+----------------+------------------------------------------------+---------+
| nthreads       | number of threads in the calculations          | 1       |
+----------------+------------------------------------------------+---------+
| version        | version of Turbomole program                   | 6.4     |
+----------------+------------------------------------------------+---------+

=====================================
Model
=====================================
BO interface for a few model sytems are provided in UNI-xMD package.

Shin-Metiu
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1D charge transfer model system proposed by Shin and Metiu (J. Chem. Phys. 1995, 102, 9285) is implemented.

.. math::

   \hat{H}_{BO}(r;R) = -\frac{1}{2}\frac{\partial^2}{\partial r^2} 
   +\frac{1}{|\frac{L}{2}-R|}&+\frac{1}{|\frac{L}{2}+R|}\nonumber\\
   -\frac{\text{erf}\left(|R-r|/R_c\right)}{|R-r|}
   -\frac{\text{erf}\left(|r-\frac{L}{2}|/R_r\right)}{|r-\frac{L}{2}|}
   &-\frac{\text{erf}\left(|r+\frac{L}{2}|/R_r\right)}{|r+\frac{L}{2}|}

Parameters :math:`L`, :math:`R_l`, :math:`R_r`, and :math:`R_c` are set to 19.0, 3.1, 4.0, 
and 5.0 in atomic units, respectively. Only lowest two BO states are calculated at a given nuclear configuration. 

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| qm_path        | path for QM binary                             | ./      |
+----------------+------------------------------------------------+---------+

=====================================
Do It Yourself
=====================================

template ~~~ blah blah
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
template..~~? idunno
