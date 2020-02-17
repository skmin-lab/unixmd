
At the very first stage of dynamics calculations, users need to make a molecule object to be investigated. The keywords to specify the molecule are below.

+------------+----------------------------------------------------+-----------+
| Keywords   | Work                                               | Default   |
+============+====================================================+===========+
| geometry   | a string containing position and velocity          |           |
+------------+----------------------------------------------------+-----------+
| nsp        | dimension of space                                 | 3         |
+------------+----------------------------------------------------+-----------+
| nstates    | the number of BO states                            | 3         |
+------------+----------------------------------------------------+-----------+
| dof        | degrees of freedom (if model is False, molecular   | None      |
|            | DoF is given.)                                     |           |
+------------+----------------------------------------------------+-----------+
| unit_pos   | unit of position (A = angstrom, au = atomic unit   | A         |
|            | [bohr])                                            |           |
+------------+----------------------------------------------------+-----------+
|            | unit of velocity (au = atomic unit, A/ps = angstrom| au        |
| unit_vel   | /ps, A/fs = angstrom/fs)                           |           |
+------------+----------------------------------------------------+-----------+
| charge     | total charge of the system                         | 0         |
+------------+----------------------------------------------------+-----------+
| model      | is the system a model system?                      | False     |
+------------+----------------------------------------------------+-----------+


