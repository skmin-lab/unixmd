===========================
PyUNIxMD utility scripts
===========================

PyUNIxMD provides utility scripts for user convenience such as preparing dynamics input or post-processing MD data.
These scripts are located in `$PYUNIXMDHOME/util`.
Different from PyUNIxMD input, options of utility scripts are controlled by bash command line, not from script file itself.
Short guidance of each scripts can be also find from using :code:`-h` command.
Documents of detailed description of each scripts are as follows. 

input_gen.py
---------------------------
Python utility script for PyUNIxMD input generator.
In this script, input files of each trajectories are generated from 1) sampled geometry files and 2) running script.
The important thing is that sampled geometry files must be named in "sample_(number).xyz" and prepared running script must be read geometry from "geom.xyz" file.

+---------------------+----------------------------------------------------+-----------+
| Option              | Description                                        | Default   |
+=====================+====================================================+===========+
| **-d**, **-dir**    | Directory name of sampled files,                   | *Sampled* |
| *(string)*          | these files must be written in xyz format.         |           |
+---------------------+----------------------------------------------------+-----------+
| **-f**, **-file**   | Filename of personal running script,               | *run.py*  |
| *(string)*          | this file must be written in running script format.|           |
+---------------------+----------------------------------------------------+-----------+
| **-n**, **-ntrajs** | Total number of trajectories for sampling.         | *None*    |
| *(integer)*         |                                                    |           |
+---------------------+----------------------------------------------------+-----------+
| **-h**              | Call out help message.                             | *None*    |
|                     |                                                    |           |
+---------------------+----------------------------------------------------+-----------+

**Ex.** Making a 100 trajectory of input script, based on XYZ from `Sampled/` and running script from `run_base.py`

.. code-block:: bash

   $ python3 input_gen.py -dir Sampled -file run_base.py -ntrajs 100

motion_analysis.py
---------------------------
Python utility script for PyUNIxMD output analysis.
In this script, PyUNIxMD MOVIE.xyz output files are post-process into given geometry criterion.
Currently, three geometrical parameters can be measured: bond, angle, and dihedral angle.

1. In the bond analysis, bond length between two given atoms will be calculated from given geometry information.
2. In the angle analysis, angle between three given atoms will be calculated. here, second atom will be a vertex of angle. 
3. In the dihedral angle analysis, dihedral angle between four or six given atoms will be calculated. 
In four atom case, angle between (1,2,3),(2,3,4) plane will be calculated and dihedreal axis will be atom2-atom3. 
In six atom case, angle between (1,2,3),(4,5,6) plane will be calculated and dihedral axis will be atom3-atom4.

+------------------------+----------------------------------------------------+-----------+
| Option                 | Description                                        | Default   |
+========================+====================================================+===========+
| **-n**, **--ntrajs**   | Total number of trajectories for analysis.         | *None*    |
| *(integer)*            |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **-s**, **--nsteps**   | Total number of steps.                             | *None*    |
| *(integer)*            |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **-b**, **--bond**     | Target bond length between two atoms.              | *None*    |
| *(integer, list)*      |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **-a**, **--angle**    | Target angle between three atoms.                  | *None*    |
| *(integer, list)*      |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **-d**, **--dihedeal** | Target dihedral angle between four or six atoms.   | *None*    |
| *(integer, list)*      |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **-m**, **--mean**     | Additional option for averaging through            | *False*   |
| *(boolean)*            | total trajectories.                                |           |
+------------------------+----------------------------------------------------+-----------+
| **-h**                 | Call out help message.                             | *None*    |
|                        |                                                    |           |
+------------------------+----------------------------------------------------+-----------+

**Ex.** Analyze bond length between atom 1 and 4 in 100 trajectories, first 10 steps. Averaging through total trajectories either.

.. code-block:: bash

   $ python3 motion_analysis.py -n 100 -s 10 -a 1 4 -m

statistical_analysis.py
---------------------------
Python utility script for UNI-xMD output analysis.
In this script, PyUNIxMD output files are post-process into organized analysis data.
Currently, three statistical parameters are measured: BO population analysis, electron coherence analysis, nacme averaging.

1. BO population analysis based on the running state of each trajectories or based on density matrix of each trajectories
2. Electronic coherence analysis based on density matrix of each trajectories
3. Averaging off-diagonal non-adiabatic coupling matrix, phase is ignored with absolute value

+------------------------+----------------------------------------------------+-----------+
| Option                 | Description                                        | Default   |
+========================+====================================================+===========+
| **-n**, **-ntrajs**    | Total number of trajectories for analysis.         | *None*    |
| *(integer)*            |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **-s**, **-nsteps**    | Total number of steps.                             | *None*    |
| *(integer)*            |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **-t**, **-nstates**   | Total number of states.                            | *None*    |
| *(integer)*            |                                                    |           |
+------------------------+----------------------------------------------------+-----------+
| **-h**                 | Call out help message.                             | *None*    |
|                        |                                                    |           |
+------------------------+----------------------------------------------------+-----------+

**Ex.** Statistical analysis on 100 trajectories, first 10 steps in 2 states.

.. code-block:: bash

   $ python3 statistical_analysis.py -n 100 -s 10 -t 2

