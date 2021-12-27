
CTMQC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Coupled-trajectory mixed quantum-classical (CTMQC) method is included in PyUNIxMD package.
Electronic equation of motion in CTMQC contains "decoherence term" which is derived from exact factorization,
in addition to Eherenfest term. The propagation of nuclei and electrons is done as follows.

.. math::

   \mathbf{F}_{\nu}^{(I)}=-\sum_{i} \nabla_{\nu}\mathbf{E}_{i}^{(I)}(t) + \sum_{i\neq j} C_{i}^{(I)\ast}(t)C_{j}^{(I)}(t)(E_{i}^{(I)}(t)-E_{j}^{(I)}(t))\mathbf{d}_{ij\nu}^{(I)}(t)
   - \sum_{i}|C_{i}^{(I)}(t)|^2\left(\sum^{N_n}_{\nu'=1}\frac{2}{\hbar M_{\nu'}}\mathbf{Q}^{(I)}_{\nu'}(t)\cdot\mathbf{f}^{(I)}_{i,\nu'}(t)\right)
   \left[\sum_{j}|C_{j}^{(I)}(t)|^2\mathbf{f}_{j,\nu}^{(I)}(t)-\mathbf{f}_{i,\nu}^{(I)}(t)\right],

.. math::

    \dot C^{(I)}_k(t) = -\frac{i}{\hbar}E^{(I)}_k(t)C^{(I)}_k(t)
    - \sum_j\sum_{\nu}{\bf d}^{(I)}_{kj\nu}(t)\cdot\dot{\bf R}^{(I)}_\nu(t)C^{(I)}_j(t)
    - \sum_{\nu}^{N_n}\frac{\mathbf{Q}^{(I)}_{\nu}(t)}{\hbar M_{\nu}}\cdot\left[\sum_{j}|C^{(I)}_{j}(t)|^2\mathbf{f}^{(I)}_{j,\nu}(t)-\mathbf{f}^{(I)}_{k,\nu}(t)\right]C^{(I)}_{k}.

Detailed description of CTMQC method is in the :cite:`Agostini2016`.

.. note:: We strongly recommend that users see detailed description of parameter.

+--------------------------------+------------------------------------------------+-----------------+
| Parameters                     | Work                                           | Default         |
+================================+================================================+=================+
| **molecules**                  | Molecule object                                |                 |
| *(*:class:`Molecule`, *list)*  |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **thermostat**                 | Thermostat object                              | *None*          |
| (:class:`Thermostat`)          |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **istates**                    | Initial states                                 | *None*          |
| *(integer, list)*              |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **dt**                         | Time interval                                  | *0.5*           |
| *(double)*                     |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **nsteps**                     | Total step of nuclear propagation              | *1000*          |
| *(integer)*                    |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **nesteps**                    | Total step of electronic propagation           | *20*            |
| *(integer)*                    |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **elec_object**                | Electronic equation of motions                 | *'coefficient'* |
| *(string)*                     |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **propagator**                 | Electronic propagator                          | *'rk4'*         |
| *(string)*                     |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **l_print_dm**                 | Logical to print BO population and coherence   | *True*          |
| *(boolean)*                    |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **l_adj_nac**                  | Logical to adjust nonadiabatic coupling        | *True*          |
| *(boolean)*                    |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **rho_threshold**              | Electronic density threshold for decoherence   | *0.01*          |
| *(double)*                     | term calculation                               |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **init_coefs**                 | Initial BO coefficients                        | *None*          |
| *(double/complex, 2D list)*    |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **dist_parameter**             | Distance parameter to determine quantum        | *10.0*          |
| *(double)*                     | momentum center                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **min_sigma**                  | Minimum sigma value                            | *0.3*           |
| *(double)*                     |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **const_dist_cutoff**          | Constant distance cutoff to determine          | *None*          |
| *(double)*                     | variance                                       |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **const_center_cutoff**        | Constant distance cutoff to determine          | *None*          |
| *(double)*                     | center of quantum momentum                     |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **unit_dt**                    | Unit of time step                              | *'fs'*          |
| *(string)*                     |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **out_freq**                   | Frequency of printing output                   | *1*             |
| *(integer)*                    |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **verbosity**                  | Verbosity of output                            | *0*             | 
| *(integer)*                    |                                                |                 |
+--------------------------------+------------------------------------------------+-----------------+

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

.. note:: CTMQC requires multiple trajectories, different to other MQC methods. Therefore, The data type of **molecules**, and **istates** must be list. 
   Also, the date type of elements in **init_coefs** must be list.

- **molecules** 
  
  This parameter defines molecular information for coupled trajectories.
  The data type of **molecules** must be list of which the elements are objects for :ref:`Molecule <Objects Molecule>`.
  For example, if the number of coupled trajectories is 3, then **molecules** can be given *[mol1, mol2, mol3]* 
  where *mol1*, *mol2*, and *mol3* are objects for :ref:`Molecule <Objects Molecule>`.

\

- **istates** *(integer, list)* - Default: *None*

  The BO coefficients and BO density matrices for coupled trajectories are initialized according to this parameter. 
  The data type of **istates** must be list of which the elements are integer.
  Hence, the number of elements in **istates** must be same to the number of coupled trajectories.
  The possible range for element in **istates** is from *0* to ``molecule.nst - 1``.
  For example, if the number of coupled trajectories is 3, then **istates** can be *[0, 0, 0]* 
  which means BO coefficients of ground state for all trajectories are initially set to 1.0.

\

- **dt** *(double)* - Default: *0.5*

  This parameter determines the time interval of the nuclear time steps.
  You can select the unit of time for the dynamics with the **unit_dt** parameter.

\

- **nsteps** *(integer)* - Default: *1000*

  This parameter determines the total number of the nuclear time steps.

\

- **nesteps** *(integer)* - Default: *20*

  This parameter determines the number of electronic time steps between one nuclear time step for the integration of the electronic equation of motion.
  The electronic equation of motion is more sensitive to the time interval than the nuclear equation of motion since the electrons are much lighter than the nuclei.
  Therefore, the nuclear time step is further divided and electronic equation of motion is integrated with smaller time step.

\

- **elec_object** *(string)* - Default: *'coefficient'*

  The **elec_object** parameter determines the representation for the electronic state.
  Now, CTMQC is only vaild for *'coefficient'*.

  + *'coefficient'*: Propagates the coefficients, i.e., :math:`\{C_{i}^{(I)}(t)\}`

\

- **propagator** *(string)* - Default: *'rk4'*

  This parameter determines the numerical integration method for the electronic equation of motion.
  Currently, only the RK4 algorithm (*'rk4'*) is available.

\

- **l_print_dm** *(boolean)* - Default: *True*

  This parameter determines whether to write output files for the density matrix elements ('BOPOP', 'BOCOH') or not.
  If this option is set to *True*, then the 'BOPOP' and 'BOCOH' files are written during the dynamics.
  This option is effective only if the parameter **obj** is set to *'coefficient'* or ignored otherwise.

\

- **l_adj_nac** *(boolean)* - Default: *True*

  If this parameter is set to *True*, the signs of the NACVs are adjusted to match the phases to the previous time step during the dynamics.

\

- **rho_threshold** *(double)* - Default: *0.01*

  This parameter defines the numerical density threshold for the coherence. 
  Specifically, if the populations of two or more states are larger than this value, the electronic state is 'coherent' and the decoherence term is calculated.

\

- **init_coefs** *(double/complex, 2D list)* - Default: *None*

  This parameter defines the initial BO coefficients. 
  The data type of element in this parameter must be list of which the elements are either real or complex values.
  The length of list, which is element of **init_coefs**, should be same to ``molecule.nst``.
  For example, if ``molecule.nst`` = *2* and the number of coupled trajectories is 3, **init_coefs** can be given *[[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]]*.
  If the parameter is not given, the BO coefficients and the density matrix are initialized according to **istates**.

\

- **dist_parameter** *(double)* - Default: *10.0*

  This parameter defines distance parameter to determine variance and position of quantum momentum center.
  The distance cutoff is determined by **dist_parameter** :math:`\times` :math:`\sigma(t)`, where :math:`\sigma(t)` is time-dependent variance.
  The :math:`\sigma(t)` is constructed by trajectories in the distance cutoff.
  If a center of quantum momentum is located out of the distance cutoff, the quanum momentum is set to zero.

\

- **min_sigma** *(double)* - Default: *0.3*

  This parameter defines minimum value for time-dependent variance. 

\

- **const_dist_cutoff** *(double)* - Default: *None*

  This parameter defines constant distance cutoff to determine variance.
  The :math:`\sigma(t)` is constructed by trajectories in the distance cutoff.

\

- **const_center_cutoff** *(double)* - Default: *None*

  This parameter defines constant distance cutoff to determine center of quantum momentum.
  If a center of quantum momentum is located out of the distance cutoff, the quanum momentum is set to zero.

\

- **unit_dt** *(string)* - Default: *'fs'*

  This parameter determines the unit of time for the simulation.

  + *'fs'*: femtosecond
  + *'au'*: atomic unit

\

- **out_freq** *(integer)* - Default: *1*

  PyUNIxMD prints and writes the dynamics information at every **out_freq** time steps.

\

- **verbosity** *(integer)* - Default: *0*

  This parameter determines the verbosity of the output files and stream.  

  + **verbosity** :math:`\geq` *1*: Prints potential energy of all BO states.
  + **verbosity** :math:`\geq` *2*: Writes the NACVs ('NACV\_\ :math:`i`\_\ :math:`j`').
