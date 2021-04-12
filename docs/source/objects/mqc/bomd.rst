
BOMD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Born-Oppenheimer molecular dynamics (BOMD) is the most simplified form of the mixed
quantum-classical dynamics, considering only a *single* potential energy surface.
Briefly, the molecular wave function is approximated to a product between an eigenfunction of the
Born-Oppenheimer Hamiltonian and its corresponding nuclear wave function as follows.

.. math::

   \Psi(\underline{\underline{\mathbf{r}}},\underline{\underline{\mathbf{R}}},t) \approx
   \chi_{i}(\underline{\underline{\mathbf{R}}},t) \Phi_{i}(\underline{\underline{\mathbf{r}}};
   \underline{\underline{\mathbf{R}}}(t)),

where :math:`\underline{\underline{\mathbf{r}}}` is the collective index for electrons,
:math:`\underline{\underline{\mathbf{R}}}` is for nuclei, and the subscript :math:`i`
stands for an index of the eigenfunctions. The electronic wave function is parametrically
dependent on :math:`\underline{\underline{\mathbf{R}}}` given at each time t. Plugging
the above product into the full molecular time-dependent Schroedinger equation and
approximating the nuclear equation further within the quantum hydrodynamic formulation
:cite:`Madelung1927` :cite:`Sakurai1994`, the coupled equations are deduced. The BO quantities are given only at positions belonging to a classical trajectory :math:`\left\{\underline{\underline{\mathbf{R}}}^{(I)}(t) \right\}` from the dynamics.

.. math::

   M_{\nu} \ddot{\mathbf{R}}_{\nu}^{(I)}(t) = - \nabla_{\nu}E_{i}^{(I)}(t),

.. math::

   \hat{H}_{\mathrm{BO}}\Phi_{i}^{(I)}(t)
    = E_{i}^{(I)}(t) \Phi_{i}^{(I)}(t),

where :math:`E_{i}^{(I)}(t)` is the :math:`i`-th BO energy eigenvalue, :math:`\nu` is the index for each nucleus, and the superscript :math:`(I)` stands for the given trajectory :math:`\left\{\underline{\underline{\mathbf{R}}}^{(I)}(t) \right\}`.

In this implementation, the eigenvalue gradient (the negative BO force) of a given
target state is provided by external electronic structure softwares
or customized Hamiltonians, and the nuclear propagation is done by the Velocity-Verlet algorithm.

+------------------------+------------------------------------------------+------------+
| Parameters             | Work                                           | Default    |
+========================+================================================+============+
| **molecule**           | Molecule object                                |            |
| (:class:`Molecule`)    |                                                |            |
+------------------------+------------------------------------------------+------------+
| **thermostat**         | Thermostat object                              | *None*     |
| (:class:`Thermostat`)  |                                                |            |
+------------------------+------------------------------------------------+------------+
| **istate**             | Electronic state                               | *0*        |
| *(integer)*            |                                                |            |
+------------------------+------------------------------------------------+------------+
| **dt**                 | Time interval                                  | *0.5*      |
| *(double)*             |                                                |            |
+------------------------+------------------------------------------------+------------+
| **nsteps**             | Total step of nuclear propagation              | *1000*     |
| *(integer)*            |                                                |            |
+------------------------+------------------------------------------------+------------+
| **unit_dt**            | Unit of time interval                          | *'fs'*     |
| *(string)*             |                                                |            |
+------------------------+------------------------------------------------+------------+
| **out_freq**           | Frequency of printing output                   | *1*        |
| *(integer)*            |                                                |            |
+------------------------+------------------------------------------------+------------+
| **verbosity**          | Verbosity of output                            | *0*        | 
| *(integer)*            |                                                |            |
+------------------------+------------------------------------------------+------------+


Detailed description of the parameters
""""""""""""""""""""""""""""""""""""""""""

- **istate** *(integer)* - Default: *0* (Ground state)
  
  This parameter specifies the adiabatic state which provides the potential energy surface nuclei follow. The possible range is from *0* to ``molecule.nst - 1``.
   
\

- **dt** *(double)* - Default: *0.5*

  This parameter determines the time interval of the nuclear time steps.
  You can select the unit of time for the dynamics with the **unit_dt** parameter.

\

- **nsteps** *(integer)* - Default: *1000*

  This parameter determines the total number of the nuclear time steps.

\

- **unit_dt** *(string)* - Default: *'fs'*

  This parameter determines the unit of time for the simulation.
  
  + *'fs'*: Femtosecond
  + *'au'*: Atomic unit

\

- **out_freq** *(integer)* - Default: *1*

  PyUNIxMD prints and writes the dynamics information at every **out_freq** time steps.

\

- **verbosity** *(integer)* - Default: *0*

  This parameter determines the verbosity of the output files and stream.
  BOMD does not print or write any additional data.

