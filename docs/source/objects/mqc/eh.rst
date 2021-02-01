
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
| **l_adjnac**               | logical to adjust nonadiabatic coupling        | *True*      |
| *(boolean)*                |                                                |             |
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

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **istate** *(integer)* - Default: *0*

The initial adiabatic state is used in dynamics. 
*0* means ground state.

\

- **dt** *(double)* - Default: *0.5*

The time interval is used in dynamics.

\

- **nsteps** *(integer)* - Default: *1000*

The number of nuclei propagation steps is used in dynamics.

\

- **nesteps** *(integer)* - Default: *20*

The number of electronic propagation steps is used in dynamics.

\

- **propagation** *(string)* - Default: *'density'*

The electronic propagation scheme is used in dynamics.

+ 'coefficient': Solve time-dependent schrodinger equation in terms of adiabatic coefficient.
+ 'density': Solve time-dependent schrodinger equation in terms of adaiabatic population.

\

- **solver** *(string)* - Default: *'rk4'*

The algorithms to solve time-dependent schrodinger equation is used in dynamics.

+ 'rk4': Use Runge-Kutta 4th order.

\

- **l_pop_print** *(boolean)* - Default: *False*

Whether print adiabatic population or not when **propagation** is set to 'coefficient'.

+ True: Print adiabatic population
+ False: Do not print adiabatic population

\

- **l_adjnac** *(boolean)* - Default: *True*

Whether adjust non-adiabatic coupling (NAC) or not.

+ True: Adjust NAC
+ False: Do not adjust NAC

\

- **coefficient** *(double/complex, list)* - Default: *None*

The initial adiabatic coefficient is used in dynamics.
This should be given by list which has length same to total adiabatic state in dynamics.
The components in this can be both double and complex.

\

- **unit_dt** *(string)* - Default: *'fs'*

the time unit is used in dynamics

+ 'fs': femto second
+ 'au': atomic unit

\

- **out_freq** *(integer)* - Default: *1*

the output is printed at each **out_freq** steps.

\

- **verbosity** *(integer)* - Default: *0*

verbosity

\
