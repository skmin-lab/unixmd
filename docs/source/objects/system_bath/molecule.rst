
Molecule
-------------------------------------------

At the very first stage of dynamics calculations, you need to make
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

The keywords to specify a molecule are below.

+---------------+------------------------------------------------------+-----------+
| Keywords      | Work                                                 | Default   |
+===============+======================================================+===========+
| **geometry**  | A string containing atomic position and velocity     |           |
| *(string)*    |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **nsp**       | Dimension of space                                   | *3*       |
| *(integer)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **nstates**   | Number of BO states                                  | *3*       |
| *(integer)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **qmmm**      | Use QMMM scheme for the calculation of large systems | *False*   |
| *(boolean)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **natoms_mm** | Number of atoms in MM region                         | *None*    |
| *(integer)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **dof**       | Degrees of freedom                                   | *None*    |
| *(integer)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **unit_pos**  | Unit of atomic position                              | *'A'*     |
| *(string)*    |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **unit_vel**  | Unit of atomic velocity                              | *'au'*    |
| *(string)*    |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **charge**    | Total charge of the system                           | *0.0*     |
| *(double)*    |                                                      |           |
+---------------+------------------------------------------------------+-----------+
| **model**     | Is the system a model system?                        | *False*   |
| *(boolean)*   |                                                      |           |
+---------------+------------------------------------------------------+-----------+


Detailed description of the arguments
""""""""""""""""""""""""""""""""""""""""""

- **geometry** *(string)*

  **geometry** string contains information of the structure of the system. The structure of this string is the following.
  This argument does not have default, thus user must put a proper string into **geometry** having the following structure.

  1. The number of atoms

  2. Comment line

  3. Specification of the molecule (atomic symbol, positions, velocities)

  When QM/MM scheme is used (**qmmm** = *True*), information of MM atoms is followed by QM atoms.

\

- **nsp** *(integer)* - Default: *3*

  This argument specifies the dimension of space where the dynamics occurs. 

\

- **nstates** *(integer)* - Default: *3*

  This argument specifies the number of BO states treated in the dynamics.

\

- **qmmm** *(boolean)* - Default: *False*

  This argument determines whether QM/MM scheme is used.

\

- **natoms_mm** *(integer)* - Default: *None*

  This argument specifies the number of atoms in the MM region when **qmmm** is *True*. 

\

- **dof** *(integer)* - Default: *None*

  This argument specifies the degrees of freedom of the system. This value will be set automatically if no specific value is given.

  If **model** = *False*, it becomes :math:`3 \times \textrm{(the number of atoms)}-6` (The DoF of a non-linear molecule).

  If **model** = *True*, :math:`\textrm{(the dimension)} \times \textrm{(the number of atoms)}`.

  When the argument is present, it overrides the above defaults.

\

- **unit_pos** *(string)* - Default: *'A'*

  This argument specifies the unit of the atomic positions.

  + *'A'*: Angstrom
  + *'au'*: Atomic unit

\

- **unit_vel** *(string)* - Default: *'au'*

  This argument specifies the unit of the atomic velocities.

  + *'au'*: Atomic unit
  + *'A/ps'*: Angstrom per picosecond
  + *'A/fs'*: Angstrom per femtosecond

\

- **charge** *(double)* - Default: *0.0*

  This argument specifies the total charge of the system 

\

- **model** *(boolean)* - Default: *False*

  This argument determines whether the system is a model system or not. About model systems provided by PyUNIxMD, see the Model Systems item in the :class:`QM_calculator` section.

