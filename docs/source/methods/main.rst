
=======================
MQC
=======================
Mixed quantum-classical(MQC) dynamics is general method for explaining the variation of molecule including
electronic state through time propagation. This can be exactly solved by time-dependent Schrodinger equation
for all particles, but this solution requires enormous cost for numerical calculation so it is restricted for
very small system. MQC tried to describe larger system by consider nuclear as classical particle which follows
classical equation of motion.

UNI-xMD mainly targeted on MQC, and whole dynamics implemented in current version of UNI-xMD are subclass of
MQC class. In the MQC class, there are functions for update classical properties of nuclear.

=======================
BOMD
=======================
The Born-Oppenheimer molecular dynamics (BOMD) is the most simplified form of mixed quantum-classical dynamics, considering only a
*single* potential energy surface. Briefly, the molecular wave function is approximated to a product between an eigenfunction of the 
Born-Oppenheimer Hamiltonian and its corresponding nuclear wave function as follows.


.. centered:: :math:`\Psi(\underline{\underline{\mathbf{r}}},\underline{\underline{\mathbf{R}}},t) \approx \chi_{I}(\underline{\underline{\mathbf{R}}},t) \Phi_{I}(\underline{\underline{\mathbf{r}}}, \underline{\underline{\mathbf{R}}})`,


where :math:`\underline{\underline{\mathbf{r}}}` is the collective index for electrons, :math:`\underline{\underline{\mathbf{R}}}` is for nuclei,
and the subscript :math:`I` stands for an index of the eigenfunctions. The electronic wave function is paramagnetically dependent on 
:math:`\underline{\underline{\mathbf{R}}}` given at each time t. Plugging the above product into the full molecular time-dependent Schrodinger 
equation and approximating the nuclear equation further within the quantum hydrodynamics formulation, the coupled equation is derived as follows.


.. centered:: :math:`M_{\nu} \ddot{\mathbf{R}}_{\nu} = - \nabla_{\nu}E_{I}(\underline{\underline{\mathbf{R}}})`,


.. centered:: :math:`H_{BO}(\underline{\underline{\mathbf{r}}},\underline{\underline{\mathbf{R}}})\Phi_{I}(\underline{\underline{\mathbf{r}}},\underline{\underline{\mathbf{R}}})`,


where :math:`E_{I}` is the :math:`I`-th eigenvalue and :math:`\nu` is the degrees of freedom for nuclei.\n
\t In this implementation, the eigenvalue gradient (the negative BO force) of a given target state is provided with an external electronic structure
package or a customized Hamiltonizan, and nuclear degrees of freedom are propagated with the Velocity-Verlet algorithm.

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| istate         | initial state                                  | 0(GS)   |
+----------------+------------------------------------------------+---------+
| dt             | time interval (fs)                             | 0.5     |
+----------------+------------------------------------------------+---------+
| nsteps         | Total step of nuclear propagation              | 1000    |
+----------------+------------------------------------------------+---------+

=======================
Eherenfest
=======================
Ehrenfest dynamics, which is mean-field dynamics, evolves nuclei on averaging potential energy surfaces,
:math:`E(\underline{\underline{\bf R}}(t))=\sum_{i}\vert c_i \vert^2E_i`,
here, :math:`E_i` is i-th adiabatic energy.
The driving force is given by: :math:`\vec{F}=\sum_{i} \vec{F}_i + \sum_{i\neq j} c_ic_j(E_i-E_j)d_{ij}`, where :math:`d_{ij}` is non-adiabatic couping between i-th and j-th adiabatic state.

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| istate         | initial state                                  | 0(GS)   |
+----------------+------------------------------------------------+---------+
| dt             | time interval (fs)                             | 0.5     |
+----------------+------------------------------------------------+---------+
| nsteps         | Total step of nuclear propagation              | 1000    |
+----------------+------------------------------------------------+---------+
| nesteps        | Total step of electronic propagation           | 10000   |
+----------------+------------------------------------------------+---------+
| propagation    | propagation scheme                             | density |
+----------------+------------------------------------------------+---------+
| l_adjnac       | adjust nonadiabatic coupling                   | True    |
+----------------+------------------------------------------------+---------+

================================
Surface Hopping
================================

Surface hopping dynamics, often called as Tully's fewest switches surface hopping dynamics (FSSH) is basic method
for propagate of artificial wavepackets through time. It was introduced by Tully, J. C. in 1990, and many other
augmented has been introduced up to now. The basic algorism of FSSH has been implemented in the UNI_xMD with
following equations.

:math:`M_{v}R^{I}_{v}(t) =`

Nuclear equation of motion is expressed by Newtonian equation, F = ma. It is expressed fully with classical.
However, the electronic degrees of freedom are represented just as follows.

:math:`C^{(I)}_K(t) =`

BO coefficient for electronic propagator is derived from force acting on each surface and nonadiabatic coupling
vector d. Using this coefficient we can structure hopping probability express as follows.

:math:`P^{(I)}_{L{rightarrow}K}[t,t+{delta}t} =`

:math:`{H}` is Heaviside function and :math:`{rho}` represents electronic density matrix. In this algorism, hopping probability
to running state to all other states are considered(including running state) and roll a random dice to select next
running state. If coupling is strong enough to transit to other state, the probability will be increase, and the overall
trajectories will be transit to that state in stochastical behavior.

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| istate         | initial state                                  | 0(GS)   |
+----------------+------------------------------------------------+---------+
| dt             | time interval (fs)                             | 0.5     |
+----------------+------------------------------------------------+---------+
| nsteps         | Total step of nuclear propagation              | 1000    |
+----------------+------------------------------------------------+---------+
| nesteps        | Total step of electronic propagation           | 10000   |
+----------------+------------------------------------------------+---------+
| propagation    | propagation scheme                             | density |
+----------------+------------------------------------------------+---------+
| l_adjnac       | adjust nonadiabatic coupling                   | True    |
+----------------+------------------------------------------------+---------+

================================
DISH-XF
================================
Decoherence induced surface hopping based on exact factorization method is included in UNI-xMD package.
Electronic equation of motion in DISH-XF contains "decoherence term" which is derived from exact factorization, 
in addition to Eherenfest term.

Detailed description of DISH-XF method is in paper J. Phys. Chem. Lett. 2018, 9, 5, 1097-1104.

+----------------+------------------------------------------------------+---------+
| Keywords       | Work                                                 | Default |
+================+======================================================+=========+
| istate         | initial state                                        | 0(GS)   |
+----------------+------------------------------------------------------+---------+
| dt             | time interval (fs)                                   | 0.5     |
+----------------+------------------------------------------------------+---------+
| nsteps         | Total step of nuclear propagation                    | 1000    |
+----------------+------------------------------------------------------+---------+
| nesteps        | Total step of electronic propagation                 | 10000   |
+----------------+------------------------------------------------------+---------+
| propagation    | propagation scheme                                   | density |
+----------------+------------------------------------------------------+---------+
| l_adjnac       | adjust nonadiabatic coupling                         | True    |
+----------------+------------------------------------------------------+---------+
| threshold      | electronic density threshold for decoherence term    | 0.01    |
+----------------+------------------------------------------------------+---------+
| wsigma         | width of nuclear wave packet of auxiliary trajectory | 0.1     |
+----------------+------------------------------------------------------+---------+

