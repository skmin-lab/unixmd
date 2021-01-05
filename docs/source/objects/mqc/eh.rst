
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

+----------------------------+------------------------------------------------+-------------+
| Keywords                   | Work                                           | Default     |
+============================+================================================+=============+
| **molecule**               | molecular object                               |             |
| (:class:`Molecule`)        |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **thermostat**             | thermostat type                                | *None**     |
| (:class:`Thermostat`)      |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **istate**                 | initial state                                  | *0*         |
| *(integer)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **dt**                     | time interval (fs)                             | *0.5*       |
| *(double)*                 |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **nsteps**                 | Total step of nuclear propagation              | *1000*      |
| *(integer)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **nesteps**                | Total step of electronic propagation           | *10000*     |
| *(integer)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **propagation**            | propagation scheme                             | *'density'* |
| *(string)*                 |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **solver**                 | propagation solver                             | *'rk4'*     |
| *(string)*                 |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **l_pop_print**            | logical to print BO population and coherence   | *False*     |
| *(boolean)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **l_adjnac**               | adjust nonadiabatic coupling                   | *True*      |
| *(boolean)*                | au = atomic unit)                              |             |
+----------------------------+------------------------------------------------+-------------+
| **coefficient**            | initial BO coefficient                         | *None*      |
| *(double/complex, list)*   |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **unit_dt**                | unit of time step (fs = femtosecond,           | *'fs'*      |
| *(string)*                 | au = atomic unit)                              |             |
+----------------------------+------------------------------------------------+-------------+
| **out_freq**               | frequency of printing output                   | *1*         |
| *(integer)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
| **verbosity**              | verbosity of output                            | *0*         | 
| *(integer)*                |                                                |             |
+----------------------------+------------------------------------------------+-------------+
