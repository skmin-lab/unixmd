
DISH-XF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Decoherence induced surface hopping based on exact factorization (DISH-XF) :cite:`Ha2018` method is included in UNI-xMD package.
Electronic equation of motion in DISH-XF contains "decoherence term" which is derived from exact factorization,
in addition to Eherenfest term, i.e.

.. math::

    \dot C^{(I)}_K(t) =& -\frac{i}{\hbar}E^{(I)}_K(t)C^{(I)}_K(t)
    - \sum_J\sum_\nu{\bf d}^{(I)}_{KJ\nu}(t)\cdot\dot{\bf R}^{(I)}_\nu(t)C^{(I)}_J(t) \nonumber\\
    &+\sum_J\sum_\nu\frac{1}{M_\nu}\frac{\nabla_\nu|\chi|}{|\chi|}\Bigg|_{\underline{\underline{\bf R}}^{(I)}(t)}
    \cdot\left\{{\bf f}^{(I)}_{J\nu}(t)-{\bf f}^{(I)}_{K\nu}(t)\right\}|C^{(I)}_J(t)|^2 C^{(I)}_K(t)

Detailed description of DISH-XF method is in :cite:`Ha2018`

+----------------------------+------------------------------------------------------+--------------+
| Keywords (type)            | Work                                                 | Default      |
+============================+======================================================+==============+
| **molecule**               | Molecule object                                      |              |
| (:class:`Molecule`)        |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **thermostat**             | Thermostat object                                    | *None*       |
| (:class:`Thermostat`)      |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **istate**                 | Initial state                                        | *0*          |
| *(integer)*                |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **dt**                     | Time interval                                        | *0.5*        |
| *(double)*                 |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **nsteps**                 | Total step of nuclear propagation                    | *1000*       |
| *(integer)*                |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **nesteps**                | Total step of electronic propagation                 | *20*         |
| *(integer)*                |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **propagation**            | Propagation scheme                                   | *'density'*  |
| *(string)*                 |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **solver**                 | Propagation solver                                   | *'rk4'*      |
| *(string)*                 |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **l_pop_print**            | Logical to print BO population and coherence         | *False*      |
| *(boolean)*                |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **l_adjnac**               | Adjust nonadiabatic coupling to align the phases     | *True*       |
| *(boolean)*                |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **vel_rescale**            | Velocity rescaling method after successful hop       | *'momentum'* |
| *(string)*                 |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **vel_reject**             | Velocity rescaling method after frustrated hop       | *'reverse'*  |
| *(string)*                 |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **threshold**              | Electronic density threshold for decoherence term    | *0.01*       |
| *(double)*                 |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **wsigma**                 | Width of nuclear wave packet of auxiliary trajectory | *None*       |
| *(double/(double, list))*  | for auxiliary trajectories                           |              |
+----------------------------+------------------------------------------------------+--------------+
| **coefficient**            | Initial BO coefficient                               | *None*       |
| *(double/complex, list)*   |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **l_state_wise**           | Logical to use state-wise total energies             | *False*      |
| *(boolean)*                | for auxiliary trajectories                           |              |
+----------------------------+------------------------------------------------------+--------------+
| **unit_dt**                | Unit of time interval                                | *'fs'*       |
| *(string)*                 |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **out_freq**               | Frequency of printing output                         | *1*          |
| *(integer)*                |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **verbosity**              | Verbosity of output                                  | *0*          | 
| *(integer)*                |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+


Detailed description of the arguments
""""""""""""""""""""""""""""""""""""""""""

- **istate** *(integer)* - Default: *0* (Ground state)
  
  This argument specifies the initial running state. The possible range of the argument is from *0* to ``molecule.nst-1``.
   
\

- **dt** *(double)* - Default: *0.5*
  
  This argument determines the time interval of the nuclear time steps.
  You can select the unit of time for the dynamics with the argument **unit_dt**.

\

- **nsteps** *(integer)* - Default: *1000*

  This argument determines the total number of the nuclear time steps.

\

- **nesteps** *(integer)* - Default: *20*
  
  This argument determines the number of electronic time steps between one nuclear time step for the integration of the electronic equation of motion.
  The electronic equation of motion is more sensitive to the time interval than the nuclear equation of motion since the electrons are much lighter than the nuclei.
  Therefore, the nuclear time step is further divided and electronic equation of motion is integrated with smaller time step.

\

- **propagation** *(string)*- Default: *'density'*
  
  The **propagation** argument determines the representation for the electronic state.
   
  + *'density'*: Propagates the density matrix elements, i.e., :math:`\{\rho_{ij}\}`
  + *'coefficient'*: Propagates the coefficients, i.e., :math:`\{C_{i}\}`

\

- **solver** *(string)* - Default: *'rk4'*

  This argument determines the numerical integration method for the electronic equation of motion.
  Currently, only the RK4 algorithm (*'rk4'*) is available.

\

- **l_pop_print** *(boolean)* - Default: *False*
  
  This argument determines whether to write output files for density matrix elements ('BOPOP', 'BOCOH') or not.
  If this option is set to *True*, then the 'BOPOP' and 'BOCOH' files are written during the dynamics.
  This option is effective only if the argument **propagation** is set to *'coefficient'* or ignored otherwise.

\

- **l_adjnac** *(boolean)* - Default: *True* 

  If this argument is set to *True*, the signs of the NACVs are adjusted to match the phases to the previous time step during the dynamics.

\

- **vel_rescale** *(string)* - Default: *'momentum'*

  This argument determines the direction of the momentum to be adjusted after a hop to conserve the total energy.
  If there is not enough kinetic energy in this direction, the hop is rejected and the running state is switched back to the original state.
  
  + *'energy'*: Simply rescale the nuclear velocities.
  + *'momentum'*: Adjust the momentum in the direction of the NACV.
  + *'augment'*: First, the hop is evaluated as the *'momentum'*. 
    If the kinetic energy is not enough, then the hop is evaluated again as the *'energy'*. 

\
   
- **vel_reject** *(string)* - Default: *'reverse'*
  
  This argument determines the momentum rescaling method when a hop is rejected.
  
  + *'keep'*: Do nothing, keeps the nuclear velocities.
  + *'reverse'*: Reverse the momentum along the NACV.

\

- **threshold** *(double)* - Default: *0.01*

  This argument defines the numerical threshold for the coherence. 
  Specifically, if the populations of two or more states are larger than this value, the electronic state is 'coherent' and the decoherence term is calculated.

\

- **wsigma** *(double/(double, list))* - Default: *None*

  This argument defines the width of the frozen gaussian wave packet on the auxiliary trajectories.
  If a scalar value is given, all nuclei share the same width.
  Or, if a list with the length of the number of the atoms is given, atom-wise width is used.
  In this case, the order of the atoms is same as the order of the xyz format string when the molecule object is created (``molecule.symbols``).

\

- **coefficient** *(double/complex, list)* - Default: *None*

  This argument defines the initial BO coefficients.
  The elements can be either real or complex values.
  If the argument is not given, the BO coefficients and density matrix are initialized according to the initial running state.

\

- **l_state_wise** *(boolean)* - Default: *False*

  This argument determines whether the total energies of the auxiliary trajectories are different or identical.
  If this is set to *True*, auxiliary trajectories have differnt total energy, or they all have same total energy.

\

- **unit_dt** *(string)* - Default: *'fs'*

  This argument determines the unit of time for the simulation.
  
  + *'fs'*: Femtosecond
  + *'au'*: Atomic unit

\

- **out_freq** *(integer)* - Default: *1*
  
  PyUNIxMD prints and writes the dynamics information at every **out_freq** time step.

\

- **verbosity** *(integer)* - Default: *0*

  This argument determines the verbosity of the output files and stream.

  + **verbosity** :math:`\geq` *1*: Prints potential energy of all BO states.
  + **verbosity** :math:`\geq` *2*: Prints accumulated hopping probabilities and writes the NACVs ('NACV\_\ :math:`i`\_\ :math:`j`'), qauntum momentum (QMOM), 
    phase terms ('AUX_PHASE\_\ :math:`i`'), and atomic postions and velocities of the auxiliary trajectories ('AUX_MOVIE\_\ :math:`i`.xyz') where :math:`i` and :math:`j` represent BO states.
