
Surface Hopping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Surface hopping dynamics, often called as Tully's fewest switches surface hopping dynamics (FSSH) is basic method
for propagate of artificial wavepackets through time. It was introduced by Tully, J. C. in 1990 :cite:`Tully1990`, and many other
augmented has been introduced up to now. The basic algorithm of FSSH has been implemented in the UNI-xMD with
following equations.

.. math::

   M_{v}\ddot{R}^{(I)}_{v}(t) = -\nabla_{v}E^{(I)}_{L}

Nuclear equation of motion is expressed by Newtonian equation, F = ma. It is expressed fully with classical.
However, the electronic degrees of freedom are represented just as follows.

.. math::

   \dot{C}^{(I)}_{K}(t) = -{{i}\over{\hbar}}E^{(I)}_K(t)C^{(I)}_{K}(t)-\sum_{J}\sum_{v}d^{(I)}_{KJv}\cdot\dot{R}^{(I)}
   _v(t)C^{(I)}_J(t)

BO coefficient for electronic propagator is derived from force acting on each surface and nonadiabatic coupling
vector d. Using this coefficient we can structure hopping probability express as follows.

.. math::

   P^{(I)}_{L{\rightarrow}K}[t,t+{\Delta}t] = {{2\Re[\rho^{(I)}_{LK}(t)\sum_vd^{(I)}_{LKv}\cdot\dot{R}^{(I)}_v(t)]
   {\Delta}t}\over{\rho^{(I)}_{LL}(t)}}, \rho^{(I)}_{LK}=C^{(I)}_L{\cdot}C^{(I)}_K

:math:`{\rho}` represents electronic density matrix. In this algorithm, hopping probability
to running state to all other states are considered (including running state) and roll a random dice to select next
running state. If coupling is strong enough to transit to other state, the probability will be increase, and the overall
trajectories will be transit to that state in stochastical behavior.

+----------------------------+------------------------------------------------+----------------+
| Keywords                   | Work                                           | Default        |
+============================+================================================+================+
| **molecule**               | molecular object                               |                |
| (:class:`Molecule`)        |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **thermostat**             | thermostat type                                | *None*         |
| (:class:`Thermostat`)      |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **istate**                 | initial state                                  | *0*            |
| *(integer)*                |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **dt**                     | time interval (fs)                             | *0.5*          |
| *(double)*                 |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **nsteps**                 | Total step of nuclear propagation              | *1000*         |
| *(integer)*                |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **nesteps**                | Total step of electronic propagation           | *10000*        |
| *(integer)*                |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **propagation**            | propagation scheme                             | *'density'*    |
| *(string)*                 |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **solver**                 | propagation solver                             | *'rk4'*        |
| *(string)*                 |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **l_pop_print**            | logical to print BO population and coherence   | *False*        |
| *(boolean)*                |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **l_adjnac**               | adjust nonadiabatic coupling                   | *True*         |
| *(boolean)*                |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **vel_rescale**            | velocity rescaling method after successful hop | *'momentum'*   |
| *(string)*                 |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **vel_reject**             | velocity rescaling method after frustrated hop | *'reverse'*    |
| *(string)*                 |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **coefficient**            | initial BO coefficient                         | *None*         |
| *(double/complex, list)*   |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **deco_correction**        | simple decoherence correction schemes          | *None*         |
| *(string)*                 |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **edc_parameter**          | energy constant for rescaling coefficients     | *0.1*          |
| *(double)*                 | in edc                                         |                |
+----------------------------+------------------------------------------------+----------------+
| **unit_dt**                | unit of time step (fs = femtosecond,           | *'fs'*         |
| *(double)*                 | au = atomic unit)                              |                |
+----------------------------+------------------------------------------------+----------------+
| **out_freq**               | frequency of printing output                   | *1*            |
| *(integer)*                |                                                |                |
+----------------------------+------------------------------------------------+----------------+
| **verbosity**              | verbosity of output                            | *0*            | 
| *(integer)*                |                                                |                |
+----------------------------+------------------------------------------------+----------------+


Detailed description of the arguments
""""""""""""""""""""""""""""""""""""""""""

- **istate** *(integer)* - Default: *0* (Ground state)
  
  Initial running state. The possible range of the argument is from *0* to ``molecule.nstate-1``.
   
\

- **dt** *(double)* - Default: *0.5*
  
  This argument determines the time interval of the nuclear time steps.
  You can select the unit of time for the dynamics with the argument **unit_dt**.

\

- **nsteps** *(integer)* - Default: *1000*

  This argument determines the total number of the nuclear time steps.

\

- **nesteps** *(integer)* - Default: *20*
  
  Number of electronic time steps between one nuclear time step for the integration of the electronic equation of motion.
  The electronic equation of motion is more sensitive to the time interval than the nuclear equation of motion since the electrons are much lighter than the nuclei.
  Therefore, the nuclear time step is further divided and electronic equation of motion is integrated with smaller time step.

\

- **propagation** *(string)*- Default: *'density'*
  
  The **propagation** argument determines the representation for the electronic state.
   
  + *'density'*: Propagates the density matrix elements, i.e., :math:`\{\rho_{ij}\}`
  + *'coefficient'*: Propagates the coefficients, i.e., :math:`\{C_{i}\}`

\

- **solver** *(string)* - Default: *'rk4'*

  Numerical integration method for the electronic equation of motion.
  Currently, only the RK4 algorithm (*'rk4'*) is available.

\

- **l_pop_print** *(boolean)* - Default: *'False'*
  
  Determine whether write output files for density matrix elements (BOPOP, BOCOH) or not.
  If this option is set to *True*, then the BOPOP and BOCOH files are written during the dynamics.
  This option is effective only if the argument **propagation** is set to *'coefficient'* or ignored otherwise.

\

- **l_adjnac** *(boolean)* - Default: *True* 

  If this argument is set to *True*, the signs of the NACVs are adjusted to match the phases to the previous time step during the dynamics.

\

- **vel_rescale** *(string)* - Default: *'momentum'*

  Determines the direction of the momentum to be adjusted after a hop to conserve the total energy.
  If there is not enough kinetic energy in this direction, the hop is rejected and the running state is switched back to the original state.
  
  + *'energy'*: Simply rescale the nuclear velocities.
  + *'momentum'*: Adjust the momentum in the direction of the NACV.
  + *'augment'*: First, the hop is evaluated as the  *'momentum'*. 
    If the kinetic energy is not enough, then the hop is evaluated again as the *'energy'*. 

\
   
- **vel_reject** *(string)* - Default: *'reverse'*
  
  Determines the momentum rescaling method when a hop is rejected.
  
  + *'keep'*: Do nothing, keeps the nuclear velocities.
  + *'reverse'*: Reverse the momentum along the NACV.

\

- **coefficient** *(double/complex, list)* - Default: *None*

  Defines the initial density matrix.
  The elements can be either real or complex values.
  If the argument is not given, the density matrix is initialized according to the initial running state.

\

- **deco_correction** *(string)* - Default: *None*

  Determines the decoherence correction method.

  + *'edc'*: Energy based decoherence correction (EDC) scheme of Granucci et al :cite:`Granucci2010`. 
  + *'idc'*: Instantaneous decoherence correction scheme

\

- **edc_parameter** *(double)* - Default: *0.1*

  Energy parameter in the EDC equation.

\

- **unit_dt** *(string)* - Default: *'fs'*

  This argument determines the unit of time for the simulation.
  
  + *'fs'*: Femtosecond
  + *'au'*: Atomic unit

\

- **out_freq** *(integer)* - Default: *1*
  
  PyUNIxMD prints and writes the dynamics information at every **out_freq** time steps.

\

- **verbosity** *(integer)* - Default: *0*

  Determines the verbosity of the output files and stream.

  + **verbosity** :math:`\geq` 1: Prints potential energy of all BO states.
  + **verbosity** :math:`\geq` 2: Prints accumulated hopping probabilities and writes the NACVs (NACV\_\ *ist*\_\ *jst*).
