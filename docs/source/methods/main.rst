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
