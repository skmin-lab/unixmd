
Thermostat
-------------------------------------------

A thermostat object should be provided to set temperature conditions of the dynamics (even in the
case that no thermostating process is needed). The following are the options.

**Ex.** Making thermostat object

.. code-block:: python

   from thermostat import *

   bathT = rescale1(temperature=300.0, nrescale=20) # velocity rescaling thermostat

Rescale1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`Rescale1` thermostat rescales the velocities periodically during the dynamics.
The target temperature and the number of MD steps between rescalings can be specified.

+---------------------+----------------------------------------------------+-----------+
| Keywords            | Work                                               | Default   |
+=====================+====================================================+===========+
| **temperature**     | Target temperature (K) of the thermostat           | *300.0*   |
| *(double)*          |                                                    |           |
+---------------------+----------------------------------------------------+-----------+
| **nrescale**        | The number of MD steps between rescalings          | *20*      |
| *(integer)*         |                                                    |           |
+---------------------+----------------------------------------------------+-----------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **temperature** *(double)* - Default: *300.0*

  The target temperature (K) of the thermostat to be used.

\

- **nrescale** *(integer)* - Default: *20*

  The velocities are rescaled at each **nrescale** step.

Rescale2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`Rescale2` thermostat rescales the velocities when the difference between the current temperature
at a specific time step and the target temperature is beyond a specified threshold.
The target temperature and the temperature difference threshold can be specified.

+------------------+----------------------------------------------------+-----------+
| Keywords         | Work                                               | Default   |
+==================+====================================================+===========+
| **temperature**  | Target temperature (K) of the thermostat           | *300.0*   |
| *(double)*       |                                                    |           |
+------------------+----------------------------------------------------+-----------+
| **dtemperature** | Threshold temperature difference (K)               | *100.0*   |
| *(double)*       |                                                    |           |
+------------------+----------------------------------------------------+-----------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **temperature** *(double)* - Default: *300.0*

  The target temperature (K) of the thermostat to be used.

\

- **dtemperature** *(double)* - Default: *100.0*

  The velocities are rescaled when temperature difference exceeds **dtemperature**.

Berendsen
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`Berendsen` thermostat :cite:`Berendsen`  rescales the velocities by mimicking weak coupling with first-order kinetics
to an external heat bath with given temperature :math:`T_0`. The scaling factor, :math:`\alpha`, is given by:

.. math::

   \alpha^2 = 1 + \frac{dt}{\tau} (T_0 - T(t)),

where :math:`dt` and :math:`\tau` is the MD step and coupling parameter, respectively. 
Coupling parameter determines how strongly the system and the bath have coupling.

+------------------------+----------------------------------------------------+-----------+
| Keywords               | Work                                               | Default   |
+========================+====================================================+===========+
| **temperature**        | Target temperature (K) of the thermostat           | *300.0*   |
| *(double)*             |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **coupling_parameter** | The characteristic time (fs) to damp               | *None*    |
| *(double)*             | temperature toward target temperature              |           |
+------------------------+----------------------------------------------------+-----------+
| **coupling_strength**  | Dimensionless coupling strength for the thermostat | *None*    |
| *(double)*             |                                                    |           |
+------------------------+----------------------------------------------------+-----------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **temperature** *(double)* - Default: *300.0*

  The target temperature (K) of the thermostat is used in dynamics.

\

- **coupling_parameter** *(double)* - Default: *None*

  The **coupling_parameter**, :math:`\tau`, is characteristic time to damp temperature toward target temperature.
  It can be set directly as the characteristic length of time to damp temperature. The unit is femtosecond.
  Either **coupling_parameter** or **coupling_strength** should be set and only **coupling_parameter** or **coupling_strength** can be set.

\

- **coupling_strength** *(double)* - Default: *None*

  Dimensionless coupling strength for the thermostat is given by :math:`\frac{dt}{\tau}`, where :math:`dt` is the MD step :math:`\tau` is **coupling_parameter**.
  Either **coupling_parameter** or **coupling_strength** should be set and only **coupling_parameter** or **coupling_strength** can be set.

NHC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`NHC` thermostat, which is Nosé-Hoover chain thermostat :cite:`NHC`, rescales the velocities by using friction factor which comes from imaginary particles. 

+------------------------+----------------------------------------------------+-----------+
| Keywords               | Work                                               | Default   |
+========================+====================================================+===========+
| **temperature**        | Target temperature (K) of the thermostat           | *300.0*   |
| *(double)*             |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **coupling_strength**  | The coupling strength (cm\ :sup:`-1`\) for the     | *None*    |
| *(double)*             | thermostat                                         |           |
+------------------------+----------------------------------------------------+-----------+
| **time_scale**         | The coupling time scale (fs)                       | *None*    |
| *(double)*             |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **chain_length**       | The number of imaginary particles in the thermostat| *3*       |
| *(integer)*            | chain                                              |           |
+------------------------+----------------------------------------------------+-----------+
| **order**              | The order of the evolution operator                | *3*       |
| *(integer)*            |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **nsteps**             | NHC propagation step                               | *1*       |
| *(integer)*            |                                                    |           |
+------------------------+----------------------------------------------------+-----------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **temperature** *(double)* - Default: *300.0*

  The target temperature (K) of the thermostat to be used.

\

- **coupling_strength** *(double)* - Default: *None*

  The coupling strength is used in thermostat.
  This indicates frequency of oscillation of the thermostating particles.
  This is typically related to the highest vibrational mode frequency of given system.
  **coupling_strength** or **time_scale** should be set and only **coupling_strength** or **time_scale** can be set. 

\

- **time_scale** *(double)* - Default: *None*

  The coupling time scale is used in thermostat. The unit is femtosecond.
  When **time_scale** is given as :math:`t`, **coupling_strength** set to :math:`1/t`.
  **coupling_strength** or **time_scale** should be set and only **coupling_strength** or **time_scale** can be set. 

\

- **chain_length** *(integer)* - Default: *3*

  The number of imaginary particles in the thermostat chain is used in dynamics.

\

- **order** *(integer)* - Default: *3*

  The order of the evolution operator is used in dynamics. **order** should be *3* or *5*.

\

- **nsteps** *(integer)* - Default: *3*

  The propagation step in NHC thermostat. 
