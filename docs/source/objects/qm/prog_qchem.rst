
Q-Chem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Q-Chem :cite:`qchem2015` is a comprehensive ab initio quantum chemistry software for accurate predictions of molecular structures, reactivities, and vibrational, electronic and NMR spectra.

- (TD)DFT is used to provide with a potential energy and its gradient for a certain adiabatic state. In Q-Chem, analytical adiabatic energy gradients and nonadiabatic couplings are provided.

+--------+------+--------+----+-----+
|        | BOMD | SH(XF) | Eh | nac |
+========+======+========+====+=====+
| (TD)DFT| o    | o      | o  | o   |
+--------+------+--------+----+-----+

(TD)DFT
"""""""""""""""""""""""""""""""""""""

+-----------------------+------------------------------------------------+--------------+
| Keywords              | Work                                           | Default      |
+=======================+================================================+==============+
| **molecule**          | Molecule object                                |              |  
| (:class:`Molecule`)   |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **basis_set**         | Basis set information                          | *'sto-3g'*   |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **memory**            | Allocatable memory in the calculation          | *'2000'*     |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **nthreads**          | Number of threads in the calculation           | *1*          |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **functional**        | Exchange-correlation functional                | *'blyp'*     |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **scf_max_iter**      | Maximum number of SCF iterations               | *50*         |
| *(integer)*           |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **scf_rho_tol**       | Density convergence for SCF iterations         | *8*          |
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

  This argument specifies a basis set to be used in calculation.
  If you want to know the detailed list for basis sets, see the manual of the Q-Chem.

\

- **memory** *(integer)* - Default : *2000*

  This argument determines how much memory will be allocated in calculation. The unit is MB.

\

- **nthreads** *(integer)* - Default : *1*

  This argument specifies number of threads in calculation.

\

- **functional** *(string)* - Default : *'blyp'*

  This argument specifies exchange-correlation functional to be used in calculation.
  If you want to know the detailed list for exchane-correlation functional, see the manual of the Q-Chem.

\

- **scf_max_iter** *(integer)* - Default : *50*

  This argument determines maximum number of SCF iterations.

\

- **scf_rho_tol** *(integer)* - Default : *8*

  SCF is considered converged when the density error is less than :math:`10^{-\textbf{scf_rho_tol}}`.

\

- **cis_max_iter** *(integer)* - Default : *30*

  This argument determines maximum number of CIS iterations.

\

- **cis_en_tol** *(integer)* - Default : *6*

  CIS is considered converged when the energy error is less than :math:`10^{-\textbf{cis_en_tol}}`.

\

- **cpscf_max_iter** *(integer)* - Default : *30*

  This argument determines maximum number of CPSCF iterations.

\

- **cpscf_grad_tol** *(integer)* - Default : *6*

  CPSCF is considered converged when the gradient error is less than :math:`10^{-\textbf{cpscf_grad_tol}}`.

\

- **qm_path** *(string)* - Default : *'./'*

  This argument designates a path for Q-Chem install directory. 
  To execute Q-Chem binary file, the environment varialbes for Q-Chem are assigned by executing 'qcenv.sh' in Q-Chem install directory.
  Hence, You must set **qm_path** to *'/my_disk/my_name/qchem/'* not *'/my_disk/my_name/qchem/bin/'*.

\

- **version** *(string)* - Default : *'5.2'*

  This argument determines version of Q-Chem. PyUNIxMD is currently based on 5.2 version.
