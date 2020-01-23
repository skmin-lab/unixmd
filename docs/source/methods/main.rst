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
''

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
