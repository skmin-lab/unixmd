
TURBOMOLE
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Turbomole :cite:`Ahlrichs1989` is quantum chemical program package, initially developed
in the group of Prof. Dr. Reinhart Ahlrichs at the University of Karlsruhe and at the Forschungszentrum Karlsruhe.
(TD)DFT method is interfaced with the current version of PyUNIxMD.

- (TD)DFT provides analytical gradients, thus it can be used for Born-Oppenhiemer molecular dynamics (BOMD).

+---------+------+--------+----+-----+
|         | BOMD | SH(XF) | Eh | nac |
+=========+======+========+====+=====+
| (TD)DFT | o    | x      | x  | x   |
+---------+------+--------+----+-----+

(TD)DFT
"""""""""""""""""""""""""""""""""""""

+---------------------+---------------------------------------------+----------------+
| Parameters          | Work                                        | Default        |
+=====================+=============================================+================+
| **molecule**        | Molecule object                             |                |
| (:class:`Molecule`) |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **functional**      | Exchange-correlation functional information | *'b-lyp'*      |
| *(string)*          |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **basis_set**       | Basis set information                       | *'SV(P)'*      |
| *(string)*          |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **memory**          | Allocatable memory in the calculations      | *50*           |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **scf_max_iter**    | Maximum number of SCF iterations            | *50*           |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **scf_en_tol**      | Energy convergence for SCF iterations       | *6*            |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **cis_max_iter**    | Maximum number of CIS iterations            | *25*           |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **cis_en_tol**      | Energy convergence for CIS iterations       | *6*            |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **root_path**       | Path for Turbomole root directory           | *'./'*         |
| *(string)*          |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **nthreads**        | Number of threads in the calculations       | *1*            |
| *(integer)*         |                                             |                |
+---------------------+---------------------------------------------+----------------+
| **version**         | Version of Turbomole                        | *'6.4'*        |
| *(string)*          |                                             |                |
+---------------------+---------------------------------------------+----------------+

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

   qm = qm.turbomole.DFT(molecule=mol, functional='b-lyp', basis_set='SV(P)', \
       root_path='/opt/turbomole/')

   md = mqc.BOMD(molecule=mol, nsteps=100, nesteps=20, dt=0.5, unit_dt='au', istate=1)

   md.run(qm=qm)

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **functional** *(string)* - Default: *'b-lyp'*

  This parameter specifies the exchange-correlation functional used in the Turbomole calculation.
  The available options for this parameter are the same as the original ones of Turbomole.
  It is recommended to check the Turbomole manual for the detailed list of **functional**.

\

- **basis_set** *(string)* - Default: *'SV(P)'*

  This parameter specifies the basis set used in the Turbomole calculation.
  The available options for this parameter are the same as the original ones of Turbomole.
  It is recommended to check the Turbomole manual for the detailed list of **basis_set**.

\

- **memory** *(integer)* - Default: *50*

  This parameter determines how much memory will be allocated in the QM calculation. The unit is MB.

\

- **scf_max_iter** *(integer)* - Default: *50*

  This parameter determines the maximum number of the SCF iterations.

\

- **scf_en_tol** *(integer)* - Default: *6*

  SCF is considered converged when the energy error is less than :math:`10^{-\textbf{scf_en_tol}}`.

\

- **cis_max_iter** *(integer)* - Default: *25*

  This parameter determines the maximum number of CIS iterations.

\

- **cis_en_tol** *(integer)* - Default: *6*

  CIS is considered converged when the energy error is less than :math:`10^{-\textbf{cis_en_tol}}`.

\

- **root_path** *(string)* - Default: *'./'*

  This parameter designates the path for Turbomole root directory.
  The `$TURBODIR` environment variable determines the directory where Turbomole is installed, not the binary files themselves (For example, `$TURBODIR` is '/my_disk/my_name/TURBOMOLE/').
  Thus, **root_path** must be a *'`$TURBODIR`'*, not a *'`$TURBODIR`/define'*. 

\

- **nthreads** *(integer)* - Default: *1*

  This parameter specifies the number of threads for QM calculation.

\

- **version** *(string)* - Default: *'6.4'*

  This parameter determines the version of Turbomole. PyUNIxMD is currently based on version 7.0 of Turbomole.


