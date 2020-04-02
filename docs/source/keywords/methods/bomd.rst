
The Born-Oppenheimer molecular dynamics (BOMD) is the most simplified form of mixed
quantum-classical dynamics, considering only a *single* potential energy surface.
Briefly, the molecular wave function is approximated to a product between an eigenfunction of the
Born-Oppenheimer Hamiltonian and its corresponding nuclear wave function as follows.

.. math::

   \Psi(\underline{\underline{\mathbf{r}}},\underline{\underline{\mathbf{R}}},t) \approx
   \chi_{I}(\underline{\underline{\mathbf{R}}},t) \Phi_{I}(\underline{\underline{\mathbf{r}}},
   \underline{\underline{\mathbf{R}}}),

where :math:`\underline{\underline{\mathbf{r}}}` is the collective index for electrons,
:math:`\underline{\underline{\mathbf{R}}}` is for nuclei, and the subscript :math:`I`
stands for an index of the eigenfunctions. The electronic wave function is paramagnetically
dependent on :math:`\underline{\underline{\mathbf{R}}}` given at each time t. Plugging
the above product into the full molecular time-dependent Schrodinger equation and
approximating the nuclear equation further within the quantum hydrodynamics formulation
:cite:`Madelung1927` :cite:`Sakurai1994`, the coupled equation is derived as follows.

.. math::

   M_{\nu} \ddot{\mathbf{R}}_{\nu} = - \nabla_{\nu}E_{I}(\underline{\underline{\mathbf{R}}}),

.. math::

   H_{BO}(\underline{\underline{\mathbf{r}}},\underline{\underline{\mathbf{R}}})\Phi_{I}
   (\underline{\underline{\mathbf{r}}},\underline{\underline{\mathbf{R}}}) = E_{I}
   (\underline{\underline{\mathbf{R}}})\Phi_{I}(\underline{\underline{\mathbf{r}}},\underline{\underline{\mathbf{R}}}),

where :math:`E_{I}` is the :math:`I`-th eigenvalue and :math:`\nu` is the degrees of freedom for nuclei.

In this implementation, the eigenvalue gradient (the negative BO force) of a given
target state is provided with an external electronic structure
package or a customized Hamiltonian, and nuclear propagation is done by the Velocity-Verlet algorithm.

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| istate         | initial state                                  | 0(GS)   |
+----------------+------------------------------------------------+---------+
| dt             | time interval (fs)                             | 0.5     |
+----------------+------------------------------------------------+---------+
| nsteps         | Total step of nuclear propagation              | 1000    |
+----------------+------------------------------------------------+---------+

