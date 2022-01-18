
Molpro
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Molpro :cite:`Werner2012` is a comprehensive system of ab initio programs for advanced molecular electronic structure
calculations, designed and maintained by many authors. It comprises efficient and well-parallelized
programs for standard computational chemistry applications, such as DFT or many wave function based
methods. Among them, the state-averaged complete active space self consistent field (SA-CASSCF) method is interfaced to the current version of PyUNIxMD.

- (SA-)CASSCF is the most famous multi-configurational SCF (MCSCF) method.
  Since Molpro supports calculations of analytical gradients of (SA-)CASSCF states as well as nonadiabatic coupling vectors (NACVs) among them with coupled-perturbed CASSCF (CP-CASSCF) equations, excited state molecular dynamics simulations are possible.

+-------------+------+--------+----+-----+
|             | BOMD | SH(XF) | Eh | nac |
+=============+======+========+====+=====+
| (SA-)CASSCF | o    | o      | o  | o   |
+-------------+------+--------+----+-----+

(SA-)CASSCF
"""""""""""""""""""""""""""""""""""""

+----------------------+----------------------------------------------------------------+----------------+
| Parameters           | Work                                                           | Default        |
+======================+================================================================+================+
| **molecule**         | Molecule object                                                |                |  
| (:class:`Molecule`)  |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **basis_set**        | Basis set information                                          | *'sto-3g'*     |
| *(string)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **memory**           | Allocatable memory in the calculations                         | *'500m'*       |
| *(string)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **guess**            | Initial guess for (SA-)CASSCF calculations                     | *'hf'*         |
| *(string)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **guess_file**       | File containing initial guesses for (SA-)CASSCF calculations   | *'./wf.wfu'*   |
| *(string)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **scf_max_iter**     | Maximum number of HF iterations                                | *20*           |
| *(integer)*          |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **scf_en_tol**       | Energy convergence threshold for HF iterations                 | *1E-8*         |
| *(double)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **scf_rho_tol**      | Density convergence threshold for HF iterations                | *1E-6*         |
| *(double)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **mcscf_max_iter**   | Maximum number of (SA-)CASSCF iterations                       | *20*           |
| *(integer)*          |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **mcscf_en_tol**     | Energy convergence threshold for (SA-)CASSCF iterations        | *1E-8*         |
| *(double)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **mcscf_grad_tol**   | Gradient convergence threshold for (SA-)CASSCF iterations      | *1E-6*         |
| *(double)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **mcscf_step_tol**   | Step length convergence threshold for (SA-)CASSCF iterations   | *1E-2*         |
| *(double)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **active_elec**      | Number of electrons in active space                            | *2*            |
| *(integer)*          |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **active_orb**       | Number of orbitals in active space                             | *2*            |
| *(integer)*          |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **cpscf_grad_tol**   | Gradient convergence threshold for CP-CASSCF equations         | *1E-7*         |
| *(double)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **qm_path**          | Path for QM binary                                             | *'./'*         |
| *(string)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **nthreads**         | Number of threads in the calculations                          | *1*            |
| *(integer)*          |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+
| **version**          | Version of Molpro                                              | *'2015.1'*     |
| *(string)*           |                                                                |                |
+----------------------+----------------------------------------------------------------+----------------+

Example input for CASSCF
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

   qm = qm.molpro.CASSCF(molecule=mol, basis_set='6-31g*', guess='hf', \
       active_elec=2, active_orb=2, qm_path='/opt/molpro2015.1/bin/')

   md = mqc.SHXF(molecule=mol, nsteps=100, nesteps=20, dt=0.5, unit_dt='au', \
       sigma=0.1, istate=1, hop_rescale='energy', elec_object='density')

   md.run(qm=qm)

Detailed description of the parameters
''''''''''''''''''''''''''''''''''''''''''

.. note:: Please refer Molpro (version 2015.1) manual for available options and more detailed descripttion about the input of Molpro.

- **basis_set** *(string)* - Default: *'sto-3g'*

  The **basis_set** determines the basis set for atomic orbitals.

\

- **memory** *(string)* - Default: *'500m'*

  This parameter determines the memory to be allocated for the Molpro calculation.

\

- **guess** *(string)* - Default: *'hf'*

  This parameter determines the initial guess method for the (SA-)CASSCF calculations. 

  + *'hf'*: Initial guess orbitals for the (SA-)CASSCF calculations are generated from the HF calculations.
  + *'read'*: Initial guesses of orbitals and CI coefficients are read from the 'wf.wfu' file which contains the orbitals and CI coefficients calculated at the previous time step.

\

- **guess_file** *(string)* - Default: *'./wf.wfu'*
   
  The **guess_file** determines the name of the file containing orbitals for the initial guess of orbitals for the (SA-)CASSCF calculation at the first MD step.
  This parameter is effective only if **guess** = *'read'*.
  If the file does not exist, *'hf'* option is applied for the initial guess for the (SA-)CASSCF calculation at the first MD step.

\

- **scf_max_iter** *(integer)* - Default: *20*

  This parameter determines the maximum number of the HF iterations.
  
\

- **scf_en_tol** *(double)* - Default: *1E-8*

  This parameter determines the convergence threshold for the HF energy.
  
\

- **scf_rho_tol** *(double)* - Default: *1E-6*

  This parameter determines the convergence threshold for the HF density matrix.
  
\

- **mcscf_max_iter** *(integer)* - Default: *20*

  This parameter determines the maximum number of the (SA-)CASSCF interations.
  
\

- **mcscf_en_tol** *(integer)* - Default: *1E-8*

  This parameter determines the convergence threshold for the (SA-)CASSCF energy.
  
\

- **mcscf_grad_tol** *(integer)* - Default: *1E-6*

  This parameter determines the convergence threshold for the (SA-)CASSCF gradient.
  
\

- **mcscf_step_tol** *(integer)* - Default: *1E-2*

  This parameter determines the convergence threshold for (SA-)CASSCF step length.
  
\

- **active_elec** *(integer)* - Default: *2*

  This parameter determines the number of electrons to be included in the active space of the (SA-)CASSCF calculations.

\

- **active_orb** *(integer)* - Default: *2*
  
  This parameter determines the number of orbitals to be included in the active space of the (SA-)CASSCF calculations.

\

- **cpscf_grad_tol** *(double)*  - Default: *1E-7*

  This parameter determines the convergence threshold for the accuracy of the CP-MCSCF equations for the analytical gradients and NACVs of the (SA-)CASSCF states.

\

- **qm_path** *(string)* - Default: *'./'*
  
  This parameter determines the path to be searched by PyUNIxMD for the Molpro executable file, 'molpro' for the (SA-)CASSCF calculations.

\

- **nthreads** *(integer)* - Default: *1*
  
  This parameter determines the number of thread for parallel execution of Molpro.

\

- **version** *(string)* - Default: *'2015.1'*
  
  This parameter indicates the version of Molpro to be executed.
  Currently, only version 2015.1 is interfaced.

