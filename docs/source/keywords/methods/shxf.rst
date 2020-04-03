
Decoherence induced surface hopping based on exact factorization (DISH-XF) :cite:`Ha2018` method is included in UNI-xMD package.
Electronic equation of motion in DISH-XF contains "decoherence term" which is derived from exact factorization,
in addition to Eherenfest term, i.e.

.. math::

    \dot C^{(I)}_K(t) =& -\frac{i}{\hbar}E^{(I)}_K(t)C^{(I)}_K(t)
    - \sum_J\sum_\nu{\bf d}^{(I)}_{KJ\nu}(t)\cdot\dot{\bf R}^{(I)}_\nu(t)C^{(I)}_J(t) \nonumber\\
    &+\sum_J\sum_\nu\frac{1}{M_\nu}\frac{\nabla_\nu|\chi|}{|\chi|}\Bigg|_{\underline{\underline{\bf R}}^{(I)}(t)}
    \cdot\left\{{\bf f}^{(I)}_{J\nu}(t)-{\bf f}^{(I)}_{K\nu}(t)\right\}|C^{(I)}_J(t)|^2 C^{(I)}_K(t)

Detailed description of DISH-XF method is in :cite:`Ha2018`

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

