
The Born-Oppenheimer molecular dynamics (BOMD) is the most simplified form of mixed
quantum-classical dynamics, considering only a *single* potential energy surface.
Briefly, the molecular wave function is approximated to a product between an eigenfunction of the
Born-Oppenheimer Hamiltonian and its corresponding nuclear wave function as follows.

.. math::

   \Psi(\underline{\underline{\mathbf{r}}},\underline{\underline{\mathbf{R}}},t) \approx
   \chi_{i}(\underline{\underline{\mathbf{R}}},t) \Phi_{i}(\underline{\underline{\mathbf{r}}};
   \underline{\underline{\mathbf{R}}}),

where :math:`\underline{\underline{\mathbf{r}}}` is the collective index for electrons,
:math:`\underline{\underline{\mathbf{R}}}` is for nuclei, and the subscript :math:`i`
stands for an index of the eigenfunctions. The electronic wave function is parametrically
dependent on :math:`\underline{\underline{\mathbf{R}}}` given at each time t. Plugging
the above product into the full molecular time-dependent Schrodinger equation and
approximating the nuclear equation further within the quantum hydrodynamics formulation
:cite:`Madelung1927` :cite:`Sakurai1994`, the coupled equations are deduced. The BO quantities are given only at positions belonging a classical trajectory :math:`\left\{\underline{\underline{\mathbf{R}}}^{(I)}(t) \right\}` from the dynamics.

.. math::

   M_{\nu} \ddot{\mathbf{R}}_{\nu}^{(I)} = - \nabla_{\nu}E_{i}^{(I)},

.. math::

   \hat{H}_{\mathrm{BO}}\Phi_{i}^{(I)}
    = E_{i}^{(I)} \Phi_{i}^{(I)},

where :math:`E_{i}^{(I)}` is the :math:`i`-th BO energy eigenvalue, :math:`\nu` is index for each nucleus, and the superscript :math:`(I)` shows dependency on the given trajectory :math:`\left\{\underline{\underline{\mathbf{R}}}^{(I)}(t) \right\}`.

In this implementation, the eigenvalue gradient (the negative BO force) of a given
target state is provided with an external electronic structure
package or a customized Hamiltonian, and nuclear propagation is done by the Velocity-Verlet algorithm.

+----------------+------------------------------------------------+------------+
| Keywords       | Work                                           | Default    |
+================+================================================+============+
| ``istate``     | initial state                                  | ``0`` (GS) |
+----------------+------------------------------------------------+------------+
| ``dt``         | time interval (fs)                             | ``0.5``    |
+----------------+------------------------------------------------+------------+
| ``nsteps``     | Total step of nuclear propagation              | ``1000``   |
+----------------+------------------------------------------------+------------+
| ``unit_dt``    | unit of time step (fs = femtosecond,           | ``fs``     |
|                | au = atomic unit)                              |            |
+----------------+------------------------------------------------+------------+

