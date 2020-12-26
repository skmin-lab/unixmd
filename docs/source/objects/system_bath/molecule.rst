
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
| **dof**       | degrees of freedom (if model is *False*,             | *None*    |
| *(integer)*   | molecular DoF is given.)                             |           |
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

The structure of **geometry** string is the following.

1. the number of atoms

2. comment line

3. specification of the molecule (atomic symbol, positions, velocities)

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


