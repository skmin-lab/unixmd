.. _Objects Thermostat:

Thermostat
-------------------------------------------

A thermostat object should be provided to set temperature conditions of the dynamics (even in the
case that no thermostating process is needed). The following are the options.

**Ex.** Making a thermostat object

.. code-block:: python

   import thermostat

   bathT = thermostat.Rescale1(temperature=300.0, nrescale=20) # velocity rescaling thermostat

Rescale1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`Rescale1` thermostat rescales the velocities periodically during the dynamics.
The target temperature and the number of MD steps between rescalings can be specified.

+---------------------+----------------------------------------------------+-----------+
| Parameters          | Work                                               | Default   |
+=====================+====================================================+===========+
| **temperature**     | Target temperature (K) of the thermostat           | *300.0*   |
| *(double)*          |                                                    |           |
+---------------------+----------------------------------------------------+-----------+
| **nrescale**        | The number of MD steps between rescalings          | *20*      |
| *(integer)*         |                                                    |           |
+---------------------+----------------------------------------------------+-----------+

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **temperature** *(double)* - Default: *300.0*

  This parameter specifies the target temperature. The unit is K.
  **temperature** is stored as ``thermostat.temp`` in the PyUNIxMD, thus ``thermostat.temp`` will be used throughout the code and manual.

\

- **nrescale** *(integer)* - Default: *20*

  This parameter specifies the period to rescale the velocities, that is,
  the velocities are rescaled at each **nrescale** step.

Rescale2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`Rescale2` thermostat rescales the velocities when the difference between the current temperature and the target temperature is 
beyond a specified threshold. The temperature difference threshold can be specified.

+------------------+----------------------------------------------------+-----------+
| Parameters       | Work                                               | Default   |
+==================+====================================================+===========+
| **temperature**  | Target temperature (K) of the thermostat           | *300.0*   |
| *(double)*       |                                                    |           |
+------------------+----------------------------------------------------+-----------+
| **dtemperature** | Threshold temperature difference (K)               | *100.0*   |
| *(double)*       |                                                    |           |
+------------------+----------------------------------------------------+-----------+

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **temperature** *(double)* - Default: *300.0*

  This parameter specifies the target temperature in unit of K.
  **temperature** is stored as ``thermostat.temp`` in the PyUNIxMD, thus ``thermostat.temp`` will be used throughout the code and manual.

\

- **dtemperature** *(double)* - Default: *100.0*

  This parameter specifies the deviation to rescale the velocities, that is,
  the velocities are rescaled when the temperature difference exceeds **dtemperature**. The unit is K.
  **dtemperature** is stored as ``thermostat.dtemp`` in the PyUNIxMD, thus ``thermostat.dtemp`` will be used throughout the code and manual.

Berendsen
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`Berendsen` thermostat :cite:`Berendsen`  rescales the velocities by mimicking weak coupling with first-order kinetics
to an external heat bath with given temperature :math:`T_0`. The scaling factor, :math:`\alpha`, is given by:

.. math::

   \alpha^2 = 1 + \frac{dt}{\tau} (T_0 - T(t)),

where :math:`dt` and :math:`\tau` is the MD step and the coupling parameter, respectively. 
The coupling parameter determines how strongly the system and the bath have coupling.

+------------------------+----------------------------------------------------+-----------+
| Parameters             | Work                                               | Default   |
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

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **temperature** *(double)* - Default: *300.0*

  This parameter specifies the target temperature in the unit of K.
  **temperature** is stored as ``thermostat.temp`` in the PyUNIxMD, thus ``thermostat.temp`` will be used throughout the code and manual.

\

- **coupling_parameter** *(double)* - Default: *None*

  The **coupling_parameter**, :math:`\tau`, is characteristic time to damp the temperature toward the target temperature.
  It can be set directly as the characteristic length of time to damp the temperature. The unit is fs.
  Either **coupling_parameter** or **coupling_strength** should be set and only **coupling_parameter** or **coupling_strength** can be set.
  **coupling_parameter** is stored as ``thermostat.coup_prm`` in the PyUNIxMD, thus ``thermostat.coup_prm`` will be used throughout the code and manual.

\

- **coupling_strength** *(double)* - Default: *None*

  This parameter is the dimensionless coupling strength, **coupling_strength**, for the thermostat given by :math:`\frac{dt}{\tau}`, 
  where :math:`dt` is the time interval of the dynamics and :math:`\tau` is **coupling_parameter**.
  Either **coupling_parameter** or **coupling_strength** should be set and only **coupling_parameter** or **coupling_strength** can be set.
  **coupling_strength** is stored as ``thermostat.coup_str`` in the PyUNIxMD, thus ``thermostat.coup_str`` will be used throughout the code and manual.

NHC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`NHC` thermostat, which is Nosé-Hoover chain thermostat :cite:`NHC`, rescales the velocities by using a friction factor which comes from imaginary particles. 

+------------------------+----------------------------------------------------+-----------+
| Parameters             | Work                                               | Default   |
+========================+====================================================+===========+
| **temperature**        | Target temperature (K) of the thermostat           | *300.0*   |
| *(double)*             |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **coupling_strength**  | Coupling strength (cm\ :sup:`-1`\) for the         | *None*    |
| *(double)*             | thermostat                                         |           |
+------------------------+----------------------------------------------------+-----------+
| **time_scale**         | Coupling time scale (fs)                           | *None*    |
| *(double)*             |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **chain_length**       | The number of particles in the NHC                 | *3*       |
| *(integer)*            |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **order**              | The order of the evolution operator                | *3*       |
| *(integer)*            |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **nsteps**             | NHC propagation step                               | *1*       |
| *(integer)*            |                                                    |           |
+------------------------+----------------------------------------------------+-----------+

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **temperature** *(double)* - Default: *300.0*

  This parameter specifies the target temperature in the unit of K.
  **temperature** is stored as ``thermostat.temp`` in the PyUNIxMD, thus ``thermostat.temp`` will be used throughout the code and manual.

\

- **coupling_strength** *(double)* - Default: *None*

  This parameter specifies the coupling strength which indicates the oscillation frequency of the bath particles.
  The coupling strength is typically related to the highest vibrational mode frequency of a given system. The unit is cm :sup:`-1`.
  **coupling_strength** or **time_scale** should be set and only **coupling_strength** or **time_scale** can be set.
  **coupling_strength** is stored as ``thermostat.coup_str`` in the PyUNIxMD, thus ``thermostat.coup_str`` will be used throughout the code and manual.

\

- **time_scale** *(double)* - Default: *None*

  This parameter specifies the coupling time scale in the unit of fs.
  When **time_scale** is given as :math:`t`, **coupling_strength** set to :math:`1/t`.
  **coupling_strength** or **time_scale** should be set and only **coupling_strength** or **time_scale** can be set.

\

- **chain_length** *(integer)* - Default: *3*

  This parameter specifies the number of imaginary particles in the thermostat chain is used in the dynamics.

\

- **order** *(integer)* - Default: *3*

  This parameter specifies the order of the evolution operator. 
  **order** should be *3* or *5*.

\

- **nsteps** *(integer)* - Default: *1*

  This parameter specifies the propagation step in the NHC thermostat.
