
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

+----------------------------+------------------------------------------------+-----------------+
| Parameters                 | Work                                           | Default         |
+============================+================================================+=================+
| **molecule**               | Molecule object                                |                 |
| (:class:`Molecule`)        |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **thermostat**             | Thermostat object                              | *None*          |
| (:class:`Thermostat`)      |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **istates**                | Initial states                                 | *None*          |
| *(integer, list)*          |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **dt**                     | Time interval                                  | *0.5*           |
| *(double)*                 |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **nsteps**                 | Total step of nuclear propagation              | *1000*          |
| *(integer)*                |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **nesteps**                | Total step of electronic propagation           | *20*            |
| *(integer)*                |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **rho_threshold**          | Electronic density threshold for decoherence   | *0.01*          |
| *(double)*                 | term calculation                               |                 |
+----------------------------+------------------------------------------------+-----------------+
| **obj**                    | Representation for electronic state            | *'coefficient'* |
| *(string)*                 |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **propagator**             | Electronic propagator                          | *'rk4'*         |
| *(string)*                 |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **dist_threshold**         | Electronic density threshold for decoherence   | *0.25*          |
| *(double)*                 | term calculation                               |                 |
+----------------------------+------------------------------------------------+-----------------+
| **dist_cutoff**            | Electronic density threshold for decoherence   | *0.5*           |
| *(double)*                 | term calculation                               |                 |
+----------------------------+------------------------------------------------+-----------------+
| **dist_parameter**         | Electronic density threshold for decoherence   | *10.0*          |
| *(double)*                 | term calculation                               |                 |
+----------------------------+------------------------------------------------+-----------------+
| **sigma**                  | Electronic density threshold for decoherence   | *0.3*           |
| *(double)*                 | term calculation                               |                 |
+----------------------------+------------------------------------------------+-----------------+
| **l_print_dm**             | Logical to print BO population and coherence   | *True*          |
| *(boolean)*                |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **l_adj_nac**              | Logical to adjust nonadiabatic coupling        | *True*          |
| *(boolean)*                |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **init_coefs**             | Initial BO coefficients                        | *None*          |
| *(double/complex, list)*   |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **unit_dt**                | Unit of time step                              | *'fs'*          |
| *(string)*                 |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **out_freq**               | Frequency of printing output                   | *1*             |
| *(integer)*                |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+
| **verbosity**              | Verbosity of output                            | *2*             | 
| *(integer)*                |                                                |                 |
+----------------------------+------------------------------------------------+-----------------+

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **istates** *(integer)* - Default: *0* (Ground state)

  This parameter specifies the initial running state. The possible range is from *0* to ``molecule.nst - 1``.

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

- **rho_threshold** *(double)* - Default: *0.01*

  This parameter defines the numerical threshold for the coherence. 
  Specifically, if the populations of two or more states are larger than this value, the electronic state is 'coherent' and the decoherence term is calculated.

\

- **obj** *(string)* - Default: *'density'*

  The **obj** parameter determines the representation for the electronic state.

  + *'density'*: Propagates the density matrix elements, i.e., :math:`\{\rho_{ij}^{(I)}(t)\}`
  + *'coefficient'*: Propagates the coefficients, i.e., :math:`\{C_{i}^{(I)}(t)\}`

\

- **propagator** *(string)* - Default: *'rk4'*

  This parameter determines the numerical integration method for the electronic equation of motion.
  Currently, only the RK4 algorithm (*'rk4'*) is available.

\

- **dist_threshold** *(double)* - Default: *0.25*

  This parameter defines the numerical threshold for the coherence. 
  Specifically, if the populations of two or more states are larger than this value, the electronic state is 'coherent' and the decoherence term is calculated.

\

- **dist_cutoff** *(double)* - Default: *0.5*

  This parameter defines the numerical threshold for the coherence. 
  Specifically, if the populations of two or more states are larger than this value, the electronic state is 'coherent' and the decoherence term is calculated.

\

- **dist_parameter** *(double)* - Default: *10.0*

  This parameter defines the numerical threshold for the coherence. 
  Specifically, if the populations of two or more states are larger than this value, the electronic state is 'coherent' and the decoherence term is calculated.

\

- **sigma** *(double)* - Default: *0.3*

  This parameter defines the numerical threshold for the coherence. 
  Specifically, if the populations of two or more states are larger than this value, the electronic state is 'coherent' and the decoherence term is calculated.

\

- **l_print_dm** *(boolean)* - Default: *True*

  This parameter determines whether to write output files for the density matrix elements ('BOPOP', 'BOCOH') or not.
  If this option is set to *True*, then the 'BOPOP' and 'BOCOH' files are written during the dynamics.
  This option is effective only if the parameter **obj** is set to *'coefficient'* or ignored otherwise.

\

- **l_adj_nac** *(boolean)* - Default: *True*

  If this parameter is set to *True*, the signs of the NACVs are adjusted to match the phases to the previous time step during the dynamics.

\

- **init_coefs** *(double/complex, list)* - Default: *None*

  This parameter defines the initial BO coefficients.
  The elements can be either real or complex values.
  If the parameter is not given, the BO coefficients and the density matrix are initialized according to **istate**.

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
