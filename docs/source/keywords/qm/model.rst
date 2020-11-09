
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

Shin-Metiu model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1D charge transfer model system proposed by Shin and Metiu :cite:`Shin1995` is implemented.

.. math::

   \hat{H}_{BO}(r;R) = -\frac{1}{2}\frac{\partial^2}{\partial r^2}
   +\frac{1}{|\frac{L}{2}-R|}&+\frac{1}{|\frac{L}{2}+R|}\nonumber\\
   -\frac{\text{erf}\left(|R-r|/R_c\right)}{|R-r|}
   -\frac{\text{erf}\left(|r-\frac{L}{2}|/R_r\right)}{|r-\frac{L}{2}|}
   &-\frac{\text{erf}\left(|r+\frac{L}{2}|/R_l\right)}{|r+\frac{L}{2}|}

Parameters :math:`L`, :math:`R_l`, :math:`R_r`, and :math:`R_c` are set to 19.0, 4.0, 3.1,
and 5.0 in atomic units, respectively.

+----------+---------------------------------------------------+-----------+
| Keywords | Work                                              | Default   |
+==========+===================================================+===========+
| ``nx``   | the number of grid points                         | ``./``    |
+----------+---------------------------------------------------+-----------+
| ``xmin`` | lower bound of the 1D space                       | ``-20.0`` |
+----------+---------------------------------------------------+-----------+
| ``xmax`` | upper bound of the 1D space                       | ``20.0``  |
+----------+---------------------------------------------------+-----------+
| ``L``    | the distance between two fixed nuclei             | ``19.0``  |
+----------+---------------------------------------------------+-----------+
| ``Rc``   | the parameter of a moving nucleus                 | ``5.0``   |
+----------+---------------------------------------------------+-----------+
| ``Rl``   | the parameter of a fixed nucleus in the left side | ``4.0``   |
+----------+---------------------------------------------------+-----------+
| ``Rr``   | the parameter of a fixed nucleus in the right side | ``3.1``   |
+----------+---------------------------------------------------+-----------+

