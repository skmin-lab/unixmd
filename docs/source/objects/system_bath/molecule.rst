
Molecule
-------------------------------------------

At the very first stage of dynamics calculations, users need to make
a molecule object to be investigated. 

**Ex.** Making a water molecule object

.. code-block:: python

   from molecule import Molecule

   geometry = """
   3
   h2o
   O    0.00    0.00    0.00    0.00    0.00    0.00
   H    0.95   -0.55    0.00    0.00    0.00    0.00
   H   -0.95   -0.55    0.00    0.00    0.00    0.00
   """

   mol = Molecule(geometry=geometry, nstates=2)

The keywords to specify the molecule are below.

+---------------+------------------------------------------------------+-----------+
| Keywords      | Work                                                 | Default   |
+===============+======================================================+===========+
| **geometry**  | a string containing position and velocity            |           |
| *(string)*    |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **nsp**       | dimension of space                                   | *3*       |
| *(integer)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **nstates**   | number of BO states                                  | *3*       |
| *(integer)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **qmmm**      | use QMMM scheme for the calculation of large systems | *False*   |
| *(boolean)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **natoms_mm** | number of atoms in MM region                         | *None*    |
| *(integer)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **dof**       | degrees of freedom (if **model** is *False*,         | *None*    |
| *(integer)*   | the molecular DoF is given.)                         |           |
+---------------+------------------------------------------------------+-----------+
| **unit_pos**  | unit of position                                     | *'A'*     |
| *(string)*    |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **unit_vel**  | unit of velocity                                     | *'au'*    |
| *(string)*    |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **charge**    | total charge of the system                           | *0*       |
| *(integer)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **model**     | is the system a model system?                        | *False*   |
| *(boolean)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+


Detailed description of the arguments
""""""""""""""""""""""""""""""""""""""""""

- **geometry** *(string)* - Default: *None*

  **geometry** string contains information of the structure of the system. The structure of this string is the following.

  1. the number of atoms

  2. comment line

  3. specification of the molecule (atomic symbol, positions, velocities)

  When QMMM scheme is used (**qmmm** == *True*), information of MM atoms is followed by QM atoms.

\

- **nsp** *(integer)* - Default: *3*

  It specifies the dimension of space where the dynamics occurs. 

\

- **nstates** *(integer)* - Default: *3*

  It specifies the number of BO states treated in the dynamics.

\

- **qmmm** *(boolean)* - Default: *False*

  This flag determines whether QMMM scheme is used.

\

- **natoms_mm** *(integer)* - Default: *None*

  It specifies the number of atoms in the MM region when **qmmm** is *True*. 

\

- **dof** *(integer)* - Default: *None*

  It specifies the degrees of freedom of the system. This value will be set automatically if no specific value is given.

  +**model** == *False*: :math:`3 \times \textrm{(the number of atoms)}-6` (The DoF of a non-linear molecule)
  +**model** == *True*: :math:`\textrm{(the dimension)} \times \textrm{(the number of atoms)}`

  When the argument is present, it overrides the above defaults.

\

- **unit_pos** *(string)* - Default: *'A'*

  It specifies the unit of position.

  + *'A'*: Angstrom
  + *'au'*: Hartree atomic unit

\

- **unit_vel** *(string)* - Default: *'au'*

  It specifies the unit of velocity.

  + *'au'*: Hartree atomic unit
  + *'A/ps'*: Angstrom per picosecond
  + *'A/fs'*: Angstrom per femtosecond

\

- **charge** *(integer)* - Default: *0*

  It specifies the total charge of the system 

\

- **model** *(boolean)* - Default: *False*

  This flag determines whether the system is a model system or not. About model systems provided by PyUNIxMD, see the Model Systems item in the :class:`QM_calculator` section.

