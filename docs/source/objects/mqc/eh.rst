
Ehrenfest
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In Ehrenfest dynamics :cite:`Prezhdo1999`, which is mean-field dynamics, evolves nuclei on an averaged potential energy surface,

.. math::

   E(\underline{\underline{\bf R}}^{(I)}(t))=\sum_{i}\vert C_{i}^{(I)}(t) \vert^2E_{i}^{(I)}(t),

where :math:`C_{i}^{(I)}(t)` and :math:`E_{i}^{(I)}(t)` is the :math:`i`-th BO coefficient and the adiabatic energy respectively from the :math:`I`-th trajectory. The equations of the electronic coefficients are given by introducing the BO basis expansion, :math:`\Phi^{(I)}(t)=\sum_{i}C_{i}^{(I)}(t)\Phi_{i}^{(I)}(t)`, to the electronic Schroedinger equation, :math:`i \hbar \partial_{t} \Phi^{(I)}(t)=\hat{H}_{\mathrm{BO}}\Phi^{(I)}(t)`:

.. math::

    \dot C^{(I)}_k(t) = -\frac{i}{\hbar}E^{(I)}_k(t)C^{(I)}_k(t)
    - \sum_j\sum_{\nu}{\bf d}^{(I)}_{kj\nu}(t)\cdot\dot{\bf R}^{(I)}_\nu(t)C^{(I)}_j(t).

Thus, the driving force is given by

.. math::

   \mathbf{F}_{\nu}^{(I)}=-\sum_{i} \nabla_{\nu}\mathbf{E}_{i}^{(I)}(t) + \sum_{i\neq j} C_{i}^{(I)\ast}(t)C_{j}^{(I)}(t)(E_{i}^{(I)}(t)-E_{j}^{(I)}(t))\mathbf{d}_{ij\nu}^{(I)}(t),

where :math:`\mathbf{d}_{ij\nu}^{(I)}(t) = \int d \underline{\underline{\mathbf{r}}}\Phi_{i}(\underline{\underline{\mathbf{r}}};\underline{\underline{\mathbf{R}}}^{(I)}(t))\nabla_{\nu}\Phi_{j}(\underline{\underline{\mathbf{r}}};\underline{\underline{\mathbf{R}}}^{(I)}(t))` is the nonadiabatic couping vector between the :math:`i`-th and the :math:`j`-th adiabatic state.

+----------------------------+------------------------------------------------+-------------+
| Parameters                 | Work                                           | Default     |
+============================+================================================+=============+
| **molecule**               | Molecule object                                |             |
| (:class:`Molecule`)        |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **thermostat**             | Thermostat object                              | *None*      |
| (:class:`Thermostat`)      |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **istate**                 | Initial state                                  | *0*         |
| *(integer)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **dt**                     | Time interval                                  | *0.5*       |
| *(double)*                 |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **nsteps**                 | Total step of nuclear propagation              | *1000*      |
| *(integer)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **nesteps**                | Total step of electronic propagation           | *20*        |
| *(integer)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **elec_object**            | Electronic equation of motions                 | *'density'* |
| *(string)*                 |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **propagator**             | Electronic propagator                          | *'rk4'*     |
| *(string)*                 |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **l_print_dm**             | Logical to print BO population and coherence   | *True*      |
| *(boolean)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **l_adj_nac**              | Logical to adjust nonadiabatic coupling        | *True*      |
| *(boolean)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **init_coef**              | Initial BO coefficient                         | *None*      |
| *(double/complex, list)*   |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **unit_dt**                | Unit of time step                              | *'fs'*      |
| *(string)*                 |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **out_freq**               | Frequency of printing output                   | *1*         |
| *(integer)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **verbosity**              | Verbosity of output                            | *0*         | 
| *(integer)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **istate** *(integer)* - Default: *0* (Ground state)

  The BO coefficient and BO density matrix are initialized according to **istate**, implying that the BO coefficient of **istate** is initially set to 1.0. 
  The possible range is from *0* to ``molecule.nst - 1``.

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

- **elec_object** *(string)* - Default: *'density'*

  The **elec_object** parameter determines the representation for the electronic state.

  + *'density'*: Propagates the density matrix elements, i.e., :math:`\{\rho_{ij}^{(I)}(t)\}`
  + *'coefficient'*: Propagates the coefficients, i.e., :math:`\{C_{i}^{(I)}(t)\}`

\

- **propagator** *(string)* - Default: *'rk4'*

  This parameter determines the numerical integration method for the electronic equation of motion.
  Currently, only the RK4 algorithm (*'rk4'*) is available.

\

- **l_print_dm** *(boolean)* - Default: *True*

  This parameter determines whether to write output files for the density matrix elements ('BOPOP', 'BOCOH') or not.
  If this option is set to *True*, then the 'BOPOP' and 'BOCOH' files are written during the dynamics.
  This option is effective only if the parameter **elec_object** is set to *'coefficient'* or ignored otherwise.

\

- **l_adj_nac** *(boolean)* - Default: *True*

  If this parameter is set to *True*, the signs of the NACVs are adjusted to match the phases to the previous time step during the dynamics.

\

- **init_coef** *(double/complex, list)* - Default: *None*

  This parameter defines the initial BO coefficients.
  The elements can be either real or complex values.
  The length of this paramter should be same to ``molecule.nst``.
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

  + **verbosity** :math:`\geq` *2*: Writes the NACVs ('NACV\_\ :math:`i`\_\ :math:`j`').

