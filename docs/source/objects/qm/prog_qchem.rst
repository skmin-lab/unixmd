
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
| Parameters            | Work                                           | Default      |
+=======================+================================================+==============+
| **molecule**          | Molecule object                                |              |  
| (:class:`Molecule`)   |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **basis_set**         | Basis set information                          | *'sto-3g'*   |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **memory**            | Allocatable memory in the calculation          | *2000*       |
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
| **scf_wf_tol**        | Wave function convergence for SCF iterations   | *8*          |
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
| **root_path**         | Path for Q-Chem root directory                 | *'./'*       |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+
| **version**           | Q-Chem version                                 | *'5.2'*      |
| *(string)*            |                                                |              |
+-----------------------+------------------------------------------------+--------------+

Example input for (TD)DFT
''''''''''''''''''''''''''''''''''''

.. code-block:: python
   :linenos:

   from molecule import Molecule
   import qm, mqc

   geom = '''
   3
   example
   O  1.14  3.77  0.00  0.00  0.00  0.00
   H  2.11  3.77  0.00  0.00  0.00  0.00
   H  0.81  4.45  0.60  0.00  0.00  0.00
   '''

   mol = Molecule(geometry=geom, ndim=3, nstates=2, unit_pos='angs')

   qm = qm.qchem.DFT(molecule=mol, functional='blyp', basis_set='sto-3g', \
       root_path='/opt/qchem/')

   md = mqc.SHXF(molecule=mol, nsteps=100, nesteps=20, dt=0.5, unit_dt='au', \
       sigma=0.1, istate=1, hop_rescale='energy', hop_reject='keep', elec_object='density')

   md.run(qm=qm)

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **basis_set** *(string)* - Default: *'sto-3g'*

  This parameter specifies the basis set to be used in the calculation.
  If you want to know the detailed list for basis sets, see the manual of the Q-Chem.

\

- **memory** *(integer)* - Default : *2000*

  This parameter determines how much memory will be allocated in the calculation. The unit is MB.

\

- **nthreads** *(integer)* - Default : *1*

  This parameter specifies the number of threads in calculation.

\

- **functional** *(string)* - Default : *'blyp'*

  This parameter specifies the exchange-correlation functional to be used in the calculation.
  If you want to know the detailed list for exchange-correlation functional, see the manual of the Q-Chem.

\

- **scf_max_iter** *(integer)* - Default : *50*

  This parameter determines the maximum number of the SCF iterations.

\

- **scf_wf_tol** *(integer)* - Default : *8*

  SCF is considered converged when the wave function error is less than :math:`10^{-\textbf{scf_wf_tol}}`.

\

- **cis_max_iter** *(integer)* - Default : *30*

  This parameter determines the maximum number of CIS iterations.

\

- **cis_en_tol** *(integer)* - Default : *6*

  CIS is considered converged when the energy error is less than :math:`10^{-\textbf{cis_en_tol}}`.

\

- **cpscf_max_iter** *(integer)* - Default : *30*

  This parameter determines the maximum number of the CPSCF iterations.

\

- **cpscf_grad_tol** *(integer)* - Default : *6*

  CPSCF is considered converged when the gradient error is less than :math:`10^{-\textbf{cpscf_grad_tol}}`.

\

- **root_path** *(string)* - Default : *'./'*

  This parameter designates the path for the Q-Chem root directory. 
  To execute Q-Chem binary file, the environment variables for Q-Chem are assigned by executing 'qcenv.sh' in the Q-Chem root directory.
  Hence, You must set **root_path** to *'/my_disk/my_name/qchem/'* not *'/my_disk/my_name/qchem/bin/'*.

\

- **version** *(string)* - Default : *'5.2'*

  This parameter determines the version of Q-Chem. PyUNIxMD is currently based on version 5.2.
