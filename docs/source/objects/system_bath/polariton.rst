.. _Objects Polariton:

Polariton
-------------------------------------------

At the very first stage of polariton dynamics, you need to make
a polariton object to be investigated.

**Ex.** Making a water polariton object with a single photon

.. code-block:: python

   from polariton import Polariton

   geometry = """
   3
   h2o
   O    0.00    0.00    0.00    0.00    0.00    0.00
   H    0.95   -0.55    0.00    0.00    0.00    0.00
   H   -0.95   -0.55    0.00    0.00    0.00    0.00
   """

   pol = Polariton(geometry=geometry, nstates=2, nphotons=1)

The parameters to specify a polariton are below.

+-------------------+------------------------------------------------------+-----------+
| Parameters        | Work                                                 | Default   |
+===================+======================================================+===========+
| **geometry**      | A string containing atomic positions and velocities  |           |
| *(string)*        |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **ndim**          | Dimension of space                                   | *3*       |
| *(integer)*       |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **nstates**       | Number of BO states                                  | *1*       |
| *(integer)*       |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **l_qmmm**        | Use the QM/MM scheme                                 | *False*   |
| *(boolean)*       |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **natoms_mm**     | Number of atoms in the MM region                     | *None*    |
| *(integer)*       |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **ndof**          | Degrees of freedom                                   | *None*    |
| *(integer)*       |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **nphotons**      | Number of quantized photons inside the cavity        | *1*       |
| *(integer)*       |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **photon_freq**   | Frequency of the photon inside the cavity            | *0.1*     |
| *(double)*        |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **field_pol_vec** | Field polarization vector                            | *None*    |
| *(double, list)*  |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **unit_pos**      | Unit of atomic positions                             | *'angs'*  |
| *(string)*        |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **unit_vel**      | Unit of atomic velocities                            | *'au'*    |
| *(string)*        |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **unit_freq**     | Unit of photon frequency                             | *'ev'*    |
| *(string)*        |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **charge**        | Total charge of the system                           | *0.0*     |
| *(double)*        |                                                      |           |
+-------------------+------------------------------------------------------+-----------+
| **l_model**       | Is the system a model system?                        | *False*   |
| *(boolean)*       |                                                      |           |
+-------------------+------------------------------------------------------+-----------+


Detailed description of the parameters
""""""""""""""""""""""""""""""""""""""""""

- **geometry** *(string)*

  The **geometry** string contains information of the structure of the system. The structure of this string is the following.

  1. The number of atoms

  2. Comment line

  3. Specification of the molecule (atomic symbol, positions, velocities)

  This parameter does not have default, thus user must put a proper string into **geometry** having the following structure.
  When the QM/MM scheme is used (**l_qmmm** = *True*), information of the MM atoms is followed by the QM atoms.

\

- **ndim** *(integer)* - Default: *3*

  This parameter specifies the dimension of space where the dynamics occurs. 

\

- **nstates** *(integer)* - Default: *1*

  This parameter specifies the number of BO states treated in the dynamics.
  **nstates** is stored as ``polariton.nst`` in the PyUNIxMD. The number of polaritonic states
  can be calculated by using **nstates** and **nphotons**, thus it is stored as ``polariton.pst``
  which will be used throughout the code and manual.

\

- **l_qmmm** *(boolean)* - Default: *False*

  This parameter determines whether to use the QM/MM scheme.

\

- **natoms_mm** *(integer)* - Default: *None*

  This parameter specifies the number of atoms in the MM region when **l_qmmm** is *True*. 
  **natoms_mm** is stored as ``polariton.nat_mm`` in the PyUNIxMD, thus ``polariton.nat_mm`` will be used throughout the code and manual.

\

- **ndof** *(integer)* - Default: *None*

  This parameter specifies the degrees of freedom of the system. This value will be set automatically if no specific value is given.

  If **l_model** = *False*, it becomes :math:`3 \times \textrm{(the number of atoms)}-6` (The DoF of a non-linear molecule).

  If **l_model** = *True*, :math:`\textrm{(the dimension)} \times \textrm{(the number of atoms)}`.

  When the argument is present, it overrides the above defaults.

\

- **nphotons** *(integer)* - Default: *1*

  This parameter specifies the number of cavity photons.

\

- **photon_freq** *(double)* - Default: *0.1*

  This parameter specifies the frequency of cavity photons.

\

- **field_pol_vec** *(double, list)* - Default: *None*

  This parameter specifies the field polarization vector. It must be specified in the running script,
  and the number of elements of **field_pol_vec** must be **ndim**.

\

- **unit_pos** *(string)* - Default: *'angs'*

  This parameter specifies the unit of atomic positions.

  + *'angs'*: Angstrom
  + *'au'*: Atomic unit

\

- **unit_vel** *(string)* - Default: *'au'*

  This parameter specifies the unit of atomic velocities.

  + *'au'*: Atomic unit
  + *'angs/ps'*: Angstrom per picosecond
  + *'angs/fs'*: Angstrom per femtosecond

\

- **unit_freq** *(string)* - Default: *'ev'*

  This parameter specifies the unit of photon frequency.

  + *'ev'*: Electronvolt
  + *'au'*: Atomic unit

\

- **charge** *(double)* - Default: *0.0*

  This parameter specifies the total charge of the system. 

\

- **l_model** *(boolean)* - Default: *False*

  This parameter determines whether the system is a model system or not. About model systems provided by PyUNIxMD, see :ref:`Model Systems <Model Systems>`.

