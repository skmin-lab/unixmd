
Molecule
-------------------------------------------

At the very first stage of dynamics calculations, users need to make
a molecule object to be investigated. The keywords to specify the molecule are below.

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
| **unit_pos**  | unit of position (A = angstrom, au = atomic unit)    | *'A'*     |
| *(string)*    |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **unit_vel**  | unit of velocity (au = atomic unit, A/ps, A/fs)      | *'au'*    |
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

  The structure of **geometry** string is the following.

  1. the number of atoms

  2. comment line

  3. specification of the molecule (atomic symbol, positions, velocities)

  When QMMM scheme is used (**qmmm** == True), information of MM atoms is followed by QM atoms.

\

- **nsp** *(integer)* - Default: *3*

  Sets the dimension of space where the dynamics occurs. 

\

- **nstates** *(integer)* - Default: *3*

  Sets the number of BO states treated in the dynamics.

\

- **qmmm** *(boolean)* - Default: *False*

  Determines whether QMMM scheme is used.

\

- **natoms_mm** *(integer)* - Default: *None*

  Sets the number of atoms in MM region when **qmmm** is True. 

\

- **dof** *(integer)* - Default: *None*

  Sets the Degrees of freedom (if **model** is *False*, the molecular DoF is given.)

\

- **unit_pos** *(string)* - Default: *'A'*

  Defines the unit of position (A = angstrom, au = atomic unit)

\

- **unit_vel** *(string)* - Default: *'au'*

  Defines the unit of velocity (au = atomic unit, A/ps, A/fs)

\

- **charge** *(integer)* - Default: *0*

  Sets the total charge of the system 

\

- **model** *(boolean)* - Default: *False*

  Determines whether the system is a model system or not. About model systems dealt with PyUNIxMD, see the Model Systems item in the QM_calculator section.

\

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


