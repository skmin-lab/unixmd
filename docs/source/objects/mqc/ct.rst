
CTMQC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Coupled-trajectory mixed quantum-classical (CTMQC) method is included in PyUNIxMD package.
Electronic equation of motion in CTMQC contains "decoherence term" which is derived from exact factorization,
in addition to Eherenfest term, i.e.

Here, add equation of motions for nuclei.

< Nuclei >

Here, add equation of motions for electrons.

< Electrons >

Detailed description of CTMQC method is in the ref.?.

.. note:: We recommend that users see detailed description for usage of paramters in CTMQC.

+--------------------------------+------------------------------------------------+-----------------+
| Parameters                     | Work                                           | Default         |
+================================+================================================+=================+
| **molecules**                  | Molecule object                                |                 |
| *(:class:`Molecule`, list)*    |                                                |                 |
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
| **elec_object**                | Representation for electronic state            | *'coefficient'* |
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
| **dist_cutoff**                | Distance cutoff for quantum momentum           | *0.5*           |
| *(double)*                     | calculation                                    |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **sigma_threshold**            | Sigma thershold for quantum momentum           | *0.25*          |
| *(double)*                     | calculation                                    |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **dist_parameter**             | Distance parameter to determine quantum        | *10.0*          |
| *(double)*                     | momentum center                                |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **sigma**                      | Sigma to determine qunatum momentum            | *0.3*           |
| *(double)*                     | center                                         |                 |
+--------------------------------+------------------------------------------------+-----------------+
| **init_coefs**                 | Initial BO coefficients                        | *None*          |
| *(double/complex, 2D list)*    |                                                |                 |
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
  The data type of **molecules** must be list of which the elements are instances for :ref:`Molecule <Objects Molecule>`.
  For example, if the number of coupled trajectories is 3, then **molecules** can be given *[mol1, mol2, mol3]*.

\

- **istates** *(integer, list)* - Default: *None*

  The BO coefficients and BO density matrices for coupled trajectories are initialized according to **istate**, implying that the BO coefficient of **istate** is initially set to 1.0. 
  The data type of this parameter must be list of which the elements are integer.
  Hence, the number of elements in **istates** must be the number of coupled trajectories.
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

- **elec_obj** *(string)* - Default: *'coefficient'*

  The **elec_object** parameter determines the representation for the electronic state.

  + *'coefficient'*: Propagates the coefficients, i.e., :math:`\{C_{i}^{(I)}(t)\}`

  Now, CTMQC is only vaild for *'coefficient'*.

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

- **dist_cutoff** *(double)* - Default: *0.5*

  This parameter defines the distance cutoff to construct Gaussian wavepacket from coupled trajectories.
  The wavepacket for an atom :math:`\nu` in a given trajectory is constructed by using trajectories in which an atom :math:`\nu'` is in **dist_cutoff** 
  from the atom :math:`\nu`.

\

- **sigma_threshold** *(double)* - Default: *0.25*

  This parameter defines the sigma threshold for quantum momentum calculation.

\

- **dist_parameter** *(double)* - Default: *10.0*

  This parameter defines distance parameter to determine position of quantum momentum center.
  if a position difference between an atom :math:`\nu` in quantum momentum center and a given trajectory is larger than **dist_parameter** :math:`\times` **sigma**, quantum momentum is set to *0.0*

\

- **sigma** *(double)* - Default: *0.3*

  This parameter defines sigma to determine position of quantum momentum center. 
  if a position difference between an atom :math:`\nu` in quantum momentum center and a given trajectory is larger than **dist_parameter** :math:`\times` **sigma**, quantum momentum is set to *0.0*

\

- **init_coefs** *(double/complex, list, list)* - Default: *None*

  This parameter defines the initial BO coefficients.
  The data type of element in this parameter must be list of which the elements are either real or complex values which means the initial coefficient for each trajecory.
  The length of list, which is element of **init_coefs**, should be same to ``molecule.nst``.
  For example, if ``molecule.nst`` = *2* and the number of coupled trajectories is 3, **init_coefs** can be given *[[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]]*.
  If the parameter is not given, the BO coefficients and the density matrix are initialized according to **istates**.

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
