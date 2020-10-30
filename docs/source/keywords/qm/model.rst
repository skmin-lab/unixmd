
QM interface for a few model sytems are provided in UNI-xMD package.

+------------+------+----+----+-----+
|            | BOMD | SH | Eh | nac |
+============+======+====+====+=====+
| SAC        | o    | o  | o  | o   |
+------------+------+----+----+-----+
| Shin Metiu | o    | o  | o  | o   |
+------------+------+----+----+-----+

simple avoided crossing (SAC) model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Shin Metiu model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1D charge transfer model system proposed by Shin and Metiu :cite:`Shin1995` is implemented.

.. math::

   \hat{H}_{BO}(r;R) = -\frac{1}{2}\frac{\partial^2}{\partial r^2}
   +\frac{1}{|\frac{L}{2}-R|}&+\frac{1}{|\frac{L}{2}+R|}\nonumber\\
   -\frac{\text{erf}\left(|R-r|/R_c\right)}{|R-r|}
   -\frac{\text{erf}\left(|r-\frac{L}{2}|/R_r\right)}{|r-\frac{L}{2}|}
   &-\frac{\text{erf}\left(|r+\frac{L}{2}|/R_r\right)}{|r+\frac{L}{2}|}

Parameters :math:`L`, :math:`R_l`, :math:`R_r`, and :math:`R_c` are set to 19.0, 3.1, 4.0,
and 5.0 in atomic units, respectively. Only lowest two BO states are calculated at a given nuclear configuration.

+----------------+------------------------------------------------+---------+
| Keywords       | Work                                           | Default |
+================+================================================+=========+
| ``qm_path``    | path for QM binary                             | ``./``  |
+----------------+------------------------------------------------+---------+

