
Q-Chem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Q-Chem :cite:`qchem2015` is a comprehensive ab initio quantum chemistry software for accurate predictions of molecular structures, reactivities, and vibrational, electronic and NMR spectra.

- (TD)DFT is used to provide with a potential energy and its gradient for a certain adiabatic state. In QChem, analytical adiabatic energy gradients and nonadiabatic couplings are provided.

+--------+------+----+----+-----+
|        | BOMD | SH | Eh | nac |
+========+======+====+====+=====+
| (TD)DFT| o    | o  | o  | o   |
+--------+------+----+----+-----+

(TD)DFT
"""""""""""""""""""""""""""""""""""""

.. note:: Our interface script is generated with 5.2 version of QChem program.
   Here, you should refer to manual of QChem program if you want to see detailed
   lists for **basis_set** and **functional** variables.



+-----------------------+------------------------------------------------+--------------+
| Keywords              | Work                                           | Default      |
+=======================+================================================+==============+
| **molecule**          | molecular object                               |              |  
| (:class:`Molecule`)   |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **basis_set**         | basis set information                          | *'sto-3g'*   |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **memory**            | allocatable memory in the calculations         | *'500m'*     |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **nthreads**          | number of threads in the calculation           | *1*          |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **functional**        | xc functional                                  | *'blyp'*     |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **scf_max_iter**      | maximum number of SCF iterations               | *50*         |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **scf_rho_tol**       | density convergence for SCF iterations         | *6*          |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **cis_max_iter**      | maximum number of CIS iterations               | *30*         |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **cis_en_tol**        | energy convergence for CIS iterations          | *6*          |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **cpscf_max_iter**    | maximum number of CP iterations                | *30*         |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **cpscf_grad_tol**    | gradient convergence for CP iterations         | *6*          |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **qm_path**           | path for QChem                                 | *'./'*       |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **version**           | QChem version                                  | *'5.2'*      |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+

