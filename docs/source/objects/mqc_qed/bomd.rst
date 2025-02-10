
BOMD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Born-Oppenheimer molecular dynamics (BOMD) is the most simplified form of the mixed
quantum-classical dynamics, considering only a *single* potential energy surface.

.. math::

   M_{\nu}\ddot{\mathbf{R}}^{(I)}_{\nu}(t) = -\nabla_{\nu}E^{(I)}_{r}(t),

In :class:`BOMD` class, the hopping concept is used to check the trivial crossing,
which can occur due to the strong light-matter interaction in polariton dynamics.
Thus, the nonadiabatic coupling vectors or transition dipole moments would be needed
according to the level of force evaluation; see :ref:`MQC_QED <Module MQC_QED>`.

+----------------------------+--------------------------------------------------+------------------+
| Parameters                 | Work                                             | Default          |
+============================+==================================================+==================+
| **polariton**              | Polariton object                                 |                  |
| (:class:`Polariton`)       |                                                  |                  |
+----------------------------+--------------------------------------------------+------------------+
| **thermostat**             | Thermostat object                                | *None*           |
| (:class:`Thermostat`)      |                                                  |                  |
+----------------------------+--------------------------------------------------+------------------+
| **istate**                 | Initial state                                    | *0*              |
| *(integer)*                |                                                  |                  |
+----------------------------+--------------------------------------------------+------------------+
| **dt**                     | Time interval                                    | *0.5*            |
| *(double)*                 |                                                  |                  |
+----------------------------+--------------------------------------------------+------------------+
| **nsteps**                 | Total step of nuclear propagation                | *1000*           |
| *(integer)*                |                                                  |                  |
+----------------------------+--------------------------------------------------+------------------+
| **l_adj_nac**              | Adjust nonadiabatic coupling to align the phases | *True*           |
| *(boolean)*                |                                                  |                  |
+----------------------------+--------------------------------------------------+------------------+
| **l_adj_tdp**              | Adjust transition dipole moments to align        | *True*           |
| *(boolean)*                | the phases                                       |                  |
+----------------------------+--------------------------------------------------+------------------+
| **unit_dt**                | Unit of time interval                            | *'fs'*           |
| *(string)*                 |                                                  |                  |
+----------------------------+--------------------------------------------------+------------------+
| **out_freq**               | Frequency of printing output                     | *1*              |
| *(integer)*                |                                                  |                  |
+----------------------------+--------------------------------------------------+------------------+
| **verbosity**              | Verbosity of output                              | *0*              | 
| *(integer)*                |                                                  |                  |
+----------------------------+--------------------------------------------------+------------------+


Detailed description of the parameters
""""""""""""""""""""""""""""""""""""""""""

- **istate** *(integer)* - Default: *0* (Ground state)
  
  This parameter specifies the initial running state. The possible range is from *0* to ``polariton.pst - 1``.
   
\

- **dt** *(double)* - Default: *0.5*
  
  This parameter determines the time interval of the nuclear time steps.
  You can select the unit of time for the dynamics with the **unit_dt** parameter.

\

- **nsteps** *(integer)* - Default: *1000*

  This parameter determines the total number of the nuclear time steps.

\

- **l_adj_nac** *(boolean)* - Default: *True* 

  If this parameter is set to *True*, the signs of the NACVs are adjusted to match the phases to the previous time step during the dynamics.

\

- **l_adj_tdp** *(boolean)* - Default: *True* 

  If this parameter is set to *True*, the signs of the TDPs are adjusted to match the phases to the previous time step during the dynamics.

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
  BOMD does not print or write any additional data.

