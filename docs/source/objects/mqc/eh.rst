
Ehrenfest
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ehrenfest dynamics :cite:`Prezhdo1999`, which is mean-field dynamics, evolves nuclei on averaging potential energy surfaces,

.. math::

   E(\underline{\underline{\bf R}}(t))=\sum_{i}\vert c_i \vert^2E_i,

where :math:`E_i` is i-th adiabatic energy and
the driving force is given by:

.. math::

   \vec{F}=\sum_{i} \vec{F}_i + \sum_{i\neq j} c_ic_j(E_i-E_j)d_{ij},

where :math:`d_{ij}` is nonadiabatic couping between i-th and j-th adiabatic state.

+--------------------+------------------------------------------------+-------------+
| Keywords           | Work                                           | Default     |
+====================+================================================+=============+
| ``molecule``       | molecular object                               |             |
+--------------------+------------------------------------------------+-------------+
| ``thermostat``     | thermostat type                                | ``None``    |
+--------------------+------------------------------------------------+-------------+
| ``istate``         | initial state                                  | ``0(GS)``   |
+--------------------+------------------------------------------------+-------------+
| ``dt``             | time interval (fs)                             | ``0.5``     |
+--------------------+------------------------------------------------+-------------+
| ``nsteps``         | Total step of nuclear propagation              | ``1000``    |
+--------------------+------------------------------------------------+-------------+
| ``nesteps``        | Total step of electronic propagation           | ``10000``   |
+--------------------+------------------------------------------------+-------------+
| ``propagation``    | propagation scheme                             | ``density`` |
+--------------------+------------------------------------------------+-------------+
| ``solver``         | propagation solver                             | ``rk4``     |
+--------------------+------------------------------------------------+-------------+
| ``l_pop_print``    | logical to print BO population and coherence   | ``False``   |
+--------------------+------------------------------------------------+-------------+
| ``l_adjnac``       | adjust nonadiabatic coupling                   | ``True``    |
+--------------------+------------------------------------------------+-------------+
| ``coefficient``    | initial BO coefficient                         | ``None``    |
+--------------------+------------------------------------------------+-------------+
| ``unit_dt``        | unit of time step (fs = femtosecond,           | ``fs``      |
|                    | au = atomic unit)                              |             |
+--------------------+------------------------------------------------+-------------+
| ``out_freq``       | frequency of printing output                   | ``1``       |
+--------------------+------------------------------------------------+-------------+
| ``verbosity``      | verbosity of output                            | ``0``       | 
+--------------------+------------------------------------------------+-------------+
