
Ehrenfest
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ehrenfest dynamics :cite:`Prezhdo1999`, which is mean-field dynamics, evolves nuclei on averaging potential energy surfaces,

.. math::

   E(\underline{\underline{\bf R}}(t))=\sum_{i}\vert c_i \vert^2E_i,

where :math:`E_i` is :math:`i`-th adiabatic energy and
the driving force is given by:

.. math::

   \vec{F}=\sum_{i} \vec{F}_i + \sum_{i\neq j} c_ic_j(E_i-E_j)d_{ij},

where :math:`d_{ij}` is nonadiabatic couping between :math:`i`-th and :math:`j`-th adiabatic state.

+----------------------------+------------------------------------------------+-------------+
| Keywords                   | Work                                           | Default     |
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
| **propagation**            | Propagation scheme                             | *'density'* |
| *(string)*                 |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **solver**                 | Propagation solver                             | *'rk4'*     |
| *(string)*                 |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **l_pop_print**            | Logical to print BO population and coherence   | *False*     |
| *(boolean)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **l_adjnac**               | Logical to adjust nonadiabatic coupling        | *True*      |
| *(boolean)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **coefficient**            | Initial BO coefficient                         | *None*      |
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

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

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

- **propagation** *(string)* - Default: *'density'*

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

- **coefficient** *(double/complex, list)* - Default: *None*

  This argument defines the initial BO coefficients.
  The elements can be either real or complex values.
  If the argument is not given, the BO coefficients and density matrix are initialized according to **istate**.

\

- **unit_dt** *(string)* - Default: *'fs'*

  This argument determines the unit of time for the simulation.

  + *'fs'*: femtosecond
  + *'au'*: atomic unit

\

- **out_freq** *(integer)* - Default: *1*

  PyUNIxMD prints and writes the dynamics information at every **out_freq** time steps.

\

- **verbosity** *(integer)* - Default: *0*

  This argument determines the verbosity of the output files and stream.  

  + **verbosity** :math:`\geq` *1*: Prints potential energy of all BO states.
  + **verbosity** :math:`\geq` *2*: Writes the NACVs ('NACV\_\ :math:`i`\_\ :math:`j`').
