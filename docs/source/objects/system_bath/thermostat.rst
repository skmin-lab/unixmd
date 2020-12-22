
Thermostat
-------------------------------------------

A thermostat object should be provided to set temperature conditions of the dynamics (even in the
case that no thermostating process is needed). The following are the options.

Rescale1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`Rescale1` thermostat rescales the velocities periodically during the dynamics.
The target temperature and the number of MD steps between rescalings can be specified.

+-----------------+----------------------------------------------------+-----------+
| Keywords        | Work                                               | Default   |
+=================+====================================================+===========+
| ``temperature`` | target temperature (K) of the thermostat           | ``300.0`` |
+-----------------+----------------------------------------------------+-----------+
| ``nrescale``    | the number of MD steps between rescalings          | ``20``    |
+-----------------+----------------------------------------------------+-----------+

Rescale2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`Rescale2` thermostat rescales the velocities when the difference between the current temperature
at a specific time step and the target temperature is beyond a specified thereshold.
The target temperature and the temperature difference threshold can be specified.

+------------------+----------------------------------------------------+-----------+
| Keywords         | Work                                               | Default   |
+==================+====================================================+===========+
| ``temperature``  | target temperature (K) of the thermostat           | ``300.0`` |
+------------------+----------------------------------------------------+-----------+
| ``dtemperature`` | threshold temperature difference (K)               | ``100.0`` |
+------------------+----------------------------------------------------+-----------+

Berendsen
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`Berendsen` thermostat :cite:`Berendsen`  rescales the velocities by mimicking weak coupling with first-order kinetics
to an external heat bath with given temperature.

.. note:: Either ``coupling_parameter`` or ``coupling_strength`` should be set and only ``coupling_parameter`` or ``coupling_strength`` can be set.

+------------------------+----------------------------------------------------+-----------+
| Keywords               | Work                                               | Default   |
+========================+====================================================+===========+
| ``temperature``        | target temperature (K) of the thermostat           | ``300.0`` |
+------------------------+----------------------------------------------------+-----------+
| ``coupling_parameter`` | the characteristic time (fs) to damp               | ``None``  |
|                        | temperature toward target temperature              |           |
+------------------------+----------------------------------------------------+-----------+
| ``coupling_strength``  | dimensionless coupling strength for the thermostat | ``None``  |
+------------------------+----------------------------------------------------+-----------+

NHC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:class:`NHC` thermostat, which is Nosé-Hoover chain thermostat :cite:`NHC` , rescales the velocities by using frition factor, which comes from imaginary particles. 

.. note:: Either ``coupling_strength`` or ``time_scale`` should be set and only ``coupling_strength`` or ``time_scale`` can be set. 
   ``order`` should be ``3`` or ``5``.
   
+------------------------+----------------------------------------------------+-----------+
| Keywords               | Work                                               | Default   |
+========================+====================================================+===========+
| ``temperature``        | target temperature (K) of the thermostat           | ``300.0`` |
+------------------------+----------------------------------------------------+-----------+
| ``coupling_strength``  | the coupling strength (cm\ :sup:`-1`\) for the     | ``None``  |
|                        | thermostat                                         |           |
+------------------------+----------------------------------------------------+-----------+
| ``time_scale``         | the coupling time scale (fs)                       | ``None``  |
+------------------------+----------------------------------------------------+-----------+
| ``chain_length``       | the number of imaginary particles in the thermostat| ``3``     |
|                        | chain                                              |           |
+------------------------+----------------------------------------------------+-----------+
| ``order``              | the order of the evolution operator                | ``3``     |
+------------------------+----------------------------------------------------+-----------+
| ``nsteps``             | NHC propagation step                               | ``1``     |
+------------------------+----------------------------------------------------+-----------+

**Ex.** Making thermostat objects

.. code-block:: python

   from thermostat import *

   bathT = rescale1(temperature=300.0, nrescale=20) # velocity rescaling thermostat

