
Gaussian 09
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gaussian 09 :cite:`Frisch2009` has been a standard program for electronic structure calculations of molecules.

- The (TD)DFT implementation of Gaussian 09 supports only analytical gradients, not nonadiabatic couplings.
  Instead, nonadiabatic coupling matrix elements (NACME) are calculated by using our wavefunction overlap 
  :cite:`Ryabinkin2015` routines. Thus, it can be used for adiabatic dynamics and surface hopping dynamics.

+---------+------+--------+----+-----+
|         | BOMD | SH(XF) | Eh | nac |
+=========+======+========+====+=====+
| (TD)DFT | o    | o      | x  | x   |
+---------+------+--------+----+-----+

(TD)DFT
"""""""""""""""""""""""""""""""""""""

.. note:: Our interface is tested with version Revision A.02 of Gaussian 09.

+-----------------------+---------------------------------------------+-------------------+
| Parameters            | Work                                        | Default           |
+=======================+=============================================+===================+
| **molecule**          | Molecule object                             |                   |  
| (:class:`Molecule`)   |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **functional**        | Exchange-correlation functional information | *'BLYP'*          |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **basis_set**         | Basis set information                       | *'sto-3g'*        |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **memory**            | Allocatable memory                          | *'1gb'*           |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **guess**             | Initial guess for SCF iterations            | *'Harris'*        |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **guess_file**        | Initial guess file                          | *'./g09.chk'*     |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **root_path**         | Path for Gaussian 09 root directory         | *'./'*            |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **nthreads**          | Number of threads in the calculations       | *1*               |
| *(integer)*           |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+
| **version**           | Version of Gaussian 09                      | *'Revision A.02'* |
| *(string)*            |                                             |                   |
+-----------------------+---------------------------------------------+-------------------+

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

   qm = qm.gaussian09.DFT(molecule=mol, functional='BLYP', basis_set='STO-3g', \
       guess='Harris', root_path='/opt/gaussian/')

   md = mqc.SHXF(molecule=mol, nsteps=100, nesteps=20, dt=0.5, unit_dt='au', \
       sigma=0.1, istate=1, hop_rescale='energy', hop_reject='keep', elec_object='density')

   md.run(qm=qm)

Detailed description of parameters
'''''''''''''''''''''''''''''''''''''

- **functional** *(string)* - Default: *'BLYP'*

  This parameter specifies the exchange-correlation functional used in the Gaussian 09 calculation.
  The available options for the parameter is same as the original ones of Gaussian 09.
  It is recommended to check a Gaussian 09 manual for the detailed list of **functional**.

\

- **basis_set** *(string)* - Default: *'sto-3g'*

  This parameter specifies the basis set used in the Gaussian 09 calculation.
  The available options for the parameter is same as the original ones of Gaussian 09.
  It is recommended to check a Gaussian 09 manual for the detailed list of **basis_set**.

\

- **memory** *(string)* - Default: *'1gb'*

  This parameter determines how much memory will be allocated in the QM calculation.

\

- **guess** *(string)* - Default: *'Harris'*

  This parameter determines the initial guess method for the (TD)DFT calculation.

  + *'Harris'*: Initial guess orbitals for the DFT caculations are generated by the default method of Gaussian 09 
    (diagonalizing the Harris functional :cite:`Harris1985`).
  + *'read'*: Initial guess orbitals are read from the 'g09.chk' file which contains the orbitals calculated at the previous time step.

\

- **guess_file** *(string)* - Default: *'./g09.chk'*

  The **guess_file** determines the name of the file containing orbitals for the initial guess of orbitals for the (TD)DFT calculation at the first MD step.
  This parameter is effective only if **guess** = *'read'*.
  If the file does not exist, *'Harris'* option is applied for the initial guess for the (TD)DFT calculation at the first MD step.

\

- **root_path** *(string)* - Default: *'./'*

  This parameter designates the path for the Gaussian 09 root directory, that is, the top level directory (for example, '/my_disk/my_name/gaussian09/').

\

- **nthreads** *(integer)* - Default: *1*

  This parameter specifies the number of threads for the Gaussian 09 calculation.

\

- **version** *(string)* - Default: *'Revision A.02'*

  This parameter determines the version of Gaussian 09. PyUNIxMD is currently based on version Revision A.02 of Gaussian 09.

