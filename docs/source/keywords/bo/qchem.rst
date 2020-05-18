
QChem is a comprehensive ab initio quantum chemistry software for accurate predictions of molecular structures, reactivities, and vibrational, electronic and NMR spectra.

- (TD)DFT is used to provide with a potential energy and its gradient for a certain adiabatic state. In QChem, analytical adiabatic energy gradients and nonadiabatic couplings are provided.

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| (TD)DFT| o    | o  | o  | o   |
+--------+------+----+----+-----+

(TD)DFT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Our interface script is generated with 5.2 version of QChem program.
   Here, you should refer to manual of QChem program if you want to see detailed
   lists for ``basis_set`` and ``functional`` variables.

+--------------------+------------------------------------------------+------------+
| Keywords           | Work                                           | Default    |
+====================+================================================+============+
| ``basis_set``      | basis set information                          | ``sto-3g`` |
+--------------------+------------------------------------------------+------------+
| ``memory``         | allocatable memory in the calculations         | ``500m``   |
+--------------------+------------------------------------------------+------------+
| ``nthreads``       | number of threads in the calculation           | ``1``      |
+--------------------+------------------------------------------------+------------+
| ``functional``     | xc functional                                  | ``blyp``   |
+--------------------+------------------------------------------------+------------+
| ``scf_max_iter``   | maximum number of SCF iterations               | ``20``     |
+--------------------+------------------------------------------------+------------+
| ``scf_rho_tol``    | density convergence for SCF iterations         | ``1E-6``   |
+--------------------+------------------------------------------------+------------+
| ``cis_max_iter``   | maximum number of CIS iterations               | ``30``     |
+--------------------+------------------------------------------------+------------+
| ``cis_en_tol``     | energy convergence for CIS iterations          | ``1E-6``   |
+--------------------+------------------------------------------------+------------+
| ``cpscf_max_iter`` | maximum number of CIS iterations               | ``30``     |
+--------------------+------------------------------------------------+------------+
| ``cpscf_en_tol``   | energy convergence for CIS iterations          | ``1E-6``   |
+--------------------+------------------------------------------------+------------+
| ``qm_path``        | path for QM binary                             | ``./``     |
+--------------------+------------------------------------------------+------------+
| ``version``        | version of Molpro program                      | ``5.2``    |
+--------------------+------------------------------------------------+------------+
