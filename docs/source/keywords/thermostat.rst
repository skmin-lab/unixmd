
A thermostat object should be provided to set temperature conditions of the dynamics (even in the
case that no thermostating process is needed). The following are the options.

rescale1
-------------------------------------
:class:`rescale1` thermostat rescales the velocities periodically during the dynamics.
The target temperature and the number of MD steps between rescalings can be specified.

+-----------------+----------------------------------------------------+-----------+
| Keywords        | Work                                               | Default   |
+=================+====================================================+===========+
| ``temperature`` | target temperature (K) of the thermostat           | ``300.0`` |
+-----------------+----------------------------------------------------+-----------+
| ``nrescale``    | the number of MD steps between rescalings          | ``20``    |
+-----------------+----------------------------------------------------+-----------+

rescale2
-------------------------------------
:class:`rescale2` thermostat rescales the velocities when the difference between the current temperature
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
-------------------------------------

NHC
-------------------------------------

**Ex.** Making thermostat objects

.. code-block:: python

   from thermostat import *

   bathT = rescale1(temperature=300.0, nrescale=20) # velocity rescaling thermostat

