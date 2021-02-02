
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

+-----------------------+------------------------------------------------+--------------+
| Keywords              | Work                                           | Default      |
+=======================+================================================+==============+
| **molecule**          | Molecular object                               |              |  
| (:class:`Molecule`)   |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **basis_set**         | Basis set information                          | *'sto-3g'*   |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **memory**            | Allocatable memory in the calculations         | *'500m'*     |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **nthreads**          | Number of threads in the calculation           | *1*          |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **functional**        | Xc functional                                  | *'blyp'*     |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **scf_max_iter**      | Maximum number of SCF iterations               | *50*         |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **scf_rho_tol**       | Density convergence for SCF iterations         | *6*          |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **cis_max_iter**      | Maximum number of CIS iterations               | *30*         |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **cis_en_tol**        | Energy convergence for CIS iterations          | *6*          |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **cpscf_max_iter**    | Maximum number of CP iterations                | *30*         |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **cpscf_grad_tol**    | Gradient convergence for CP iterations         | *6*          |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **qm_path**           | Path for Q-Chem                                | *'./'*       |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **version**           | Q-Chem version                                 | *'5.2'*      |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **basis_set** *(string)* - Default: *'sto-3g'*

  Sets the basis set to be used.
  These arguments are same with the original arguments in used in Q-Chem.
  If you want to know the detailed list for basis sets, see the manual of the Q-Chem program.

\

- **memory** *(integer)* - Default : *2000*

  The available total memory in unit of megabyte.

\

- **nthreads** *(integer)* - Default : *1*

  Number of threads in the calculation

\

- **functional** *(string)* - Default : *'blyp'*

  The exchange-correlation functional to be used.
  These arguments are same with the original arguments in used in Q-Chem.
  If you want to know the detailed list for basis sets, see the manual of the Q-Chem program.

\

- **scf_max_iter** *(integer)* - Default : *50*

  Maximum number of iteration for SCF

\

- **scf_rho_tol** *(integer)* - Default : *6*

  SCF is considered converged when the wave function error is less that 10 :sup:`-scf_rho_tol`

\

- **cis_max_iter** *(integer)* - Default : *30*

  Maximum number of iteration for CIS

\

- **cis_en_tol** *(integer)* - Default : *6*

  CIS is considered converged when error is less that 10 :sup:`-cis_en_tol`

\

- **cpscf_max_iter** *(integer)* - Default : *30*

  Maximum number of iteration for CPSCF

\

- **cpscf_grad_tol** *(integer)* - Default : *6*

  CPSCF is considered converged when gradient error is less that 10 :sup:`-cpscf_grad_tol`

\

- **qm_path** *(string)* - Default : *'./'*

  Path for Q-Chem install directory. The environment varialbes for Q-Chem are assigned by executing `qcenv.sh` in Q-Chem install directory.
  Hence, You must set **qm_path** to `'$YOUR_QCHEM_DIR'` not `'$YOUR_QCHEM_DIR/bin'`

\

- **version** *(string)* - Default : *'5.2'*

  Version of Q-Chem
