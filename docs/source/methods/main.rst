
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
EOMs, link module bomd.py
(tutorial, output)

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
| istep          | initial state                                  | 0(GS)   |
+----------------+------------------------------------------------+---------+
| dt             | time interval (fs)                             | 0.5     |
+----------------+------------------------------------------------+---------+
| nsteps         | Total step of nuclear propagation              | 1000    |
+----------------+------------------------------------------------+---------+
| nesteps        | Total step of electronic propagation           | 10000   |
+----------------+------------------------------------------------+---------+
| propagation    | propagation scheme                             | density |
+----------------+------------------------------------------------+---------+
| l_adjnac       | adjust non-adiabatic coupling                  | True    |
+----------------+------------------------------------------------+---------+

================================
Fewest Switch Surface Hopping
================================

Surface hopping dynamics, often called as Tully's fewest switches surface hopping dynamics(FSSH) is basic method
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



================================
DISH-XF
================================
''

Quantum momentum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
term

Phase term
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
term

Auxiliary trajectory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
explain concept, actual evaluation of QM and f

.. imports all other using toctree?
   ..toctree:
     :~~:
     molecule
     misc
     mqc/main
     bo/main
     thermostat
