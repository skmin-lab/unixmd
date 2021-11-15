
DISH-XF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Decoherence induced surface hopping based on exact factorization (DISH-XF) :cite:`Ha2018` method is included in PyUNIxMD package.
Electronic equation of motion in DISH-XF contains "decoherence term" which is derived from exact factorization,
in addition to Eherenfest term, i.e.

.. math::

    \dot C^{(I)}_k(t) =& -\frac{i}{\hbar}E^{(I)}_k(t)C^{(I)}_k(t)
    - \sum_j\sum_\nu{\bf d}^{(I)}_{kj\nu}(t)\cdot\dot{\bf R}^{(I)}_\nu(t)C^{(I)}_j(t) \nonumber\\
    &+\sum_j\sum_\nu\frac{1}{M_\nu}\frac{\nabla_\nu|\chi|}{|\chi|}\Bigg|_{\underline{\underline{\bf R}}^{(I)}(t)}
    \cdot\left\{{\bf f}^{(I)}_{j\nu}(t)-{\bf f}^{(I)}_{k\nu}(t)\right\}|C^{(I)}_j(t)|^2 C^{(I)}_k(t)

Detailed description of DISH-XF method is in :cite:`Ha2018`

+----------------------------+------------------------------------------------------+--------------+
| Parameters                 | Work                                                 | Default      |
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
| **elec_object**            | Electronic equation of motions                       | *'density'*  |
| *(string)*                 |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **propagator**             | Electronic propagator                                | *'rk4'*      |
| *(string)*                 |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **l_print_dm**             | Logical to print BO population and coherence         | *True*       |
| *(boolean)*                |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **l_adj_nac**              | Adjust nonadiabatic coupling to align the phases     | *True*       |
| *(boolean)*                |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **hop_rescale**            | Velocity rescaling method after successful hop       | *'augment'*  |
| *(string)*                 |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **hop_reject**             | Velocity rescaling method after frustrated hop       | *'reverse'*  |
| *(string)*                 |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **rho_threshold**          | Electronic density threshold for decoherence term    | *0.01*       |
| *(double)*                 | calculation                                          |              |
+----------------------------+------------------------------------------------------+--------------+
| **sigma**                  | Width of nuclear wave packet of auxiliary trajectory | *None*       |
| *(double/(double, list))*  |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **init_coef**              | Initial BO coefficient                               | *None*       |
| *(double/complex, list)*   |                                                      |              |
+----------------------------+------------------------------------------------------+--------------+
| **l_econs_state**          | Logical to use identical total energies              | *True*       |
| *(boolean)*                | for all auxiliary trajectories                       |              |
+----------------------------+------------------------------------------------------+--------------+
| **aux_econs_viol**         | How to treat trajectories violating the total energy | *'fix'*      |
| *(string)*                 | conservation                                         |              |
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

Detailed description of the parameters
""""""""""""""""""""""""""""""""""""""""""

- **istate** *(integer)* - Default: *0* (Ground state)

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

- **elec_object** *(string)*- Default: *'density'*

  The **elec_object** parameter determines the representation for the electronic state.

  + *'density'*: Propagates the density matrix elements, i.e., :math:`\{\rho_{ij}^{(I)}(t)\}`
  + *'coefficient'*: Propagates the coefficients, i.e., :math:`\{C_{i}^{(I)}(t)\}`

\

- **propagator** *(string)* - Default: *'rk4'*

  This parameter determines the numerical integration method for the electronic equation of motion.
  Currently, only the RK4 algorithm (*'rk4'*) is available.

\

- **l_print_dm** *(boolean)* - Default: *True*

  This parameter determines whether to write output files for density matrix elements ('BOPOP', 'BOCOH') or not.
  If this option is set to *True*, then the 'BOPOP' and 'BOCOH' files are written during the dynamics.
  This option is effective only if the **elec_object** parameter is set to *'coefficient'* or ignored otherwise.

\

- **l_adj_nac** *(boolean)* - Default: *True* 

  If this parameter is set to *True*, the signs of the NACVs are adjusted to match the phases to the previous time step during the dynamics.

\

- **hop_rescale** *(string)* - Default: *'augment'*

  This parameter determines the direction of the momentum to be adjusted after a hop to conserve the total energy.
  If there is not enough kinetic energy in this direction, the hop is rejected and the running state is switched back to the original state.

  + *'energy'*: Simply rescale the nuclear velocities.
  + *'momentum'*: Adjust the momentum in the direction of the NACV.
  + *'augment'*: First, the hop is evaluated as the *'momentum'*. 
    If the kinetic energy is not enough, then the hop is evaluated again as the *'energy'*. 

\
   
- **hop_reject** *(string)* - Default: *'reverse'*

  This parameter determines the momentum rescaling method when a hop is rejected.

  + *'keep'*: Do nothing, keeps the nuclear velocities.
  + *'reverse'*: Reverse the momentum along the NACV.

\

- **rho_threshold** *(double)* - Default: *0.01*

  This parameter defines the numerical threshold for the coherence. 
  Specifically, if the populations of two or more states are larger than this value, the electronic state is 'coherent' and the decoherence term is calculated.

\

- **sigma** *(double/(double, list))* - Default: *None*

  This parameter defines the width (:math:`\sigma_\nu`) of the frozen Gaussian nuclear densities (:math:`|\chi_k|^2`) 
  on the auxiliary trajectories (:math:`\underline{\underline{\textbf{R}}}_{k}`) where 
  the total nuclear denisity (:math:`|\chi|^2`) is a linear combination of the densities on the auxiliary trajectories as follows,

  .. math::
     |\chi|^2 = \sum_{k}|\chi_{k}|^2 = \sum_{k}N_{k}\prod^{N_{atom}}_\nu 
              \exp\left(-\dfrac{|\textbf{R}^{(I)}_\nu-\textbf{R}_{k,\nu}|^2}{2\sigma^2_{\nu}}\right).

  If a scalar value is given, all nuclei share the same width.
  Or, if a list of values with the length of the number of the atoms is given, an atom-wise width is used.
  In this case, the order of the atoms is the same as the order of the XYZ format string when the molecule object is created (``molecule.symbols``).

\

- **init_coef** *(double/complex, list)* - Default: *None*

  This parameter defines the initial BO coefficients.
  The elements can be either real or complex values.
  The length of this paramter should be same to ``molecule.nst``.
  If the argument is not given, the BO coefficients and density matrix are initialized according to the initial running state.

\

- **l_econs_state** *(boolean)* - Default: *True*

  This parameter determines whether the total energies of all auxiliary trajectories are identical or not.
  If this is set to *True*, auxiliary trajectories have same total energy, or they all have different total energy.
  In various system, *True* is recommended for **l_econs_state**.

\

- **aux_econs_viol** *(string)* - Default: *'fix'*

  This parameter determines how to deal with auxiliary trajectories violating the total energy conservation law.
  The velocity of an auxiliary trajectory is given as the velocity of the true nuclear trajectory multiplied by a factor determined from the total energy conservation condition, i.e.
  
  .. math::
     \underline{\underline{\dot{\textbf{R}}}}_{k} = \underline{\underline{\dot{\textbf{R}}}}^{(I)}
       \sqrt{\dfrac{E^k_{tot}-E_k^{(I)}}{\sum_{\nu}\frac{1}{2}M_{\nu}|\dot{\textbf{R}}^{(I)}_{\nu}|^2}}

  When :math:`E_{tot}^k-E^{(I)}_k < 0`, the auxiliary trajectory is either fixed or destroyed, depending on the given value of this parameter.

  + *'fix'*: Fix the auxiliary trajectory until decoherence.
  + *'collapse'*: Destroy the auxiliary trajectory, collapse the corresponding coefficient/density to zero, and renormalize. 
  
\

- **unit_dt** *(string)* - Default: *'fs'*

  This parameter determines the unit of time for the simulation.

  + *'fs'*: Femtosecond
  + *'au'*: Atomic unit

\

- **out_freq** *(integer)* - Default: *1*

  PyUNIxMD prints and writes the dynamics information at every **out_freq** time step.

\

- **verbosity** *(integer)* - Default: *0*

  This parameter determines the verbosity of the output files and stream.

  + **verbosity** :math:`\geq` *1*: Prints accumulated hopping probabilities and random numbers,
    and writes decoherence terms in time-derivative of BO populations to DOTPOPDEC.
  + **verbosity** :math:`\geq` *2*: Writes the NACVs ('NACV\_\ :math:`i`\_\ :math:`j`'), qauntum momentum (QMOM), 
    phase terms ('AUX_PHASE\_\ :math:`i`'), and atomic postions and velocities of the auxiliary trajectories
    ('AUX_MOVIE\_\ :math:`i`.xyz') where :math:`i` and :math:`j` represent BO states.

