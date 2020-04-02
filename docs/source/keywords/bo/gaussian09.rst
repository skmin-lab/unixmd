
Gaussian09 :cite:`Frisch2009` has been a standard program for electronic structure calculations of molecules.
The only BOMD using the DFT option is available with Gaussian09 in the current version of UNI-xMD,
because it doesn't explicitly provide with nonadiabatic coupling vectors.
Numerical calculation of the coupling elements using the CI overlap :cite:`Ryabinkin2015` is on progress, which allows the other dynamics options.

- (TD)DFT is used to provide with a potential energy and its gradient for a certain adiabatic state.

+---------+------+----+----+-----+
|         | BOMD | SH | Eh | nac |
+=========+======+====+====+=====+
| (TD)DFT | o    | x  | x  | x   |
+---------+------+----+----+-----+

(TD)DFT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with Revision A.02 version of Gaussian09 program.
   Please refer to the manual for the detailed lists for **basis_set** and **functional** variable.

+----------------+------------------------------------------------+---------------+
| Keywords       | Work                                           | Default       |
+================+================================================+===============+
| basis_set      | basis set information                          | sto-3g        |
+----------------+------------------------------------------------+---------------+
| memory         | allocatable memory in the calculations         | 500m          |
+----------------+------------------------------------------------+---------------+
| functional     | the level of DFT theory                        | BLYP          |
+----------------+------------------------------------------------+---------------+
| G09_root_path  | path for G09 root                              | /opt/gaussian |
+----------------+------------------------------------------------+---------------+
| nthreads       | number of threads in the calculations          | 1             |
+----------------+------------------------------------------------+---------------+
| version        | version of Gaussian09 program                  | Revision A.02 |
+----------------+------------------------------------------------+---------------+

