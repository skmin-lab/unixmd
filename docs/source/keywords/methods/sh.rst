
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

   P^{(I)}_{L{\rightarrow}K}[t,t+{\Delta}t] = {{2H[\rho^{(I)}_{LK}(t)\sum_vd^{(I)}_{LKv}\cdot\dot{R}^{(I)}_v(t)]
   {\Delta}t}\over{\rho^{(I)}_{LL}(t)}}, \rho^{(I)}_{LK}=C^{(I)}_L{\cdot}C^{(I)}_K

:math:`{H}` is Heaviside function and :math:`{\rho}` represents electronic density matrix. In this algorithm, hopping probability
to running state to all other states are considered (including running state) and roll a random dice to select next
running state. If coupling is strong enough to transit to other state, the probability will be increase, and the overall
trajectories will be transit to that state in stochastical behavior.

+--------------------+------------------------------------------------+--------------+
| Keywords           | Work                                           | Default      |
+====================+================================================+==============+
| ``istate``         | initial state                                  | ``0(GS)``    |
+--------------------+------------------------------------------------+--------------+
| ``dt``             | time interval (fs)                             | ``0.5``      |
+--------------------+------------------------------------------------+--------------+
| ``nsteps``         | Total step of nuclear propagation              | ``1000``     |
+--------------------+------------------------------------------------+--------------+
| ``nesteps``        | Total step of electronic propagation           | ``10000``    |
+--------------------+------------------------------------------------+--------------+
| ``propagation``    | propagation scheme                             | ``density``  |
+--------------------+------------------------------------------------+--------------+
| ``solver``         | propagation solver                             | ``rk4``      |
+--------------------+------------------------------------------------+--------------+
| ``l_pop_print``    | logical to print BO population and coherence   | ``False``    |
+--------------------+------------------------------------------------+--------------+
| ``l_adjnac``       | adjust nonadiabatic coupling                   | ``True``     |
+--------------------+------------------------------------------------+--------------+
| ``vel_rescale``    | velocity rescaling method after successful hop | ``momentum`` |
+--------------------+------------------------------------------------+--------------+
| ``vel_reject``     | velocity rescaling method after frustrated hop | ``reverse``  |
+--------------------+------------------------------------------------+--------------+
| ``vel_augment``    | force hop occurs with simple rescaling         | ``False``    |
+--------------------+------------------------------------------------+--------------+
| ``coefficient``    | initial BO coefficient                         | ``None``     |
+--------------------+------------------------------------------------+--------------+
| ``unit_dt``        | unit of time step (fs = femtosecond,           | ``fs``       |
|                    | au = atomic unit)                              |              |
+--------------------+------------------------------------------------+--------------+

