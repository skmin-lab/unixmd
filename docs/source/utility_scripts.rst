===========================
PyUNIxMD utility scripts
===========================

PyUNIxMD provides utility scripts for user convenience such as preparing dynamics running script and post-processing MD data.
These scripts are located in '$PYUNIXMDHOME/util/'.
Different from PyUNIxMD running script, options of utility scripts are controlled by command line, not from script file itself.
Short description of each script can be also found from using :code:`-h` command.
Documents of detailed description of each script are as follows. 

input_gen.py
---------------------------
Python utility script for PyUNIxMD running script generator.
In this script, running script of each trajectory are generated from 1) sampled geometry files and 2) running script.
After running a script, 'TRAJ_(number)/' directories will be generated. Each directory contains extended XYZ coordinate file and PyUNIxMD running script.
The important thing is that sampled geometry files must be named in 'sample_(number).xyz' and prepared running script must read geometry from 'geom.xyz' file.

**Ex.** Example of 'run.py' which reads geometry from 'geom.xyz'.

.. code-block:: python

   from molecule import Molecule

   ...

   g = open('geom.xyz', 'r')
   geom = g.read()
   g.close()

   ...

   mol = Molecule(geometry=geom, ...)

+---------------------+----------------------------------------------------------------+
| Option              | Description                                                    |
+=====================+================================================================+
| **-d**, **-dir**    | Directory name of sampled files, these files must be written   |
|                     | in extended XYZ format. Default value is 'Sampled/'            |
+---------------------+----------------------------------------------------------------+
| **-f**, **-file**   | Filename of basic running script template, this file must be   |
|                     | written in running script format. Default value is 'run.py'    |
+---------------------+----------------------------------------------------------------+
| **-n**, **-ntrajs** | Total number of trajectories for sampling.                     |
|                     |                                                                |
+---------------------+----------------------------------------------------------------+
| **-h**              | Call out help message.                                         |
|                     |                                                                |
+---------------------+----------------------------------------------------------------+

**Ex.** Making 100 trajectories of running scripts, based on XYZ from 'Sampled/' and running script from 'run_base.py'.

.. code-block:: bash

   $ ls
   Sampled/ run_base.py

   $ python3 input_gen.py -dir Sampled/ -file run_base.py -ntrajs 100

After running a script, 100 input directories with name 'TRAJ_001/' to 'TRAJ_100/' will be made.
Each directory contains running script and extended XYZ file. 

.. code-block:: bash

   $ ls
   Sampled/ run_base.py TRAJ_001/ TRAJ_002/ ... TRAJ_100/

   $ ls TRAJ_001/ 
   geom.xyz run.py

statistical_analysis.py
---------------------------
Python utility script for PyUNIxMD output analysis.
In this script, PyUNIxMD output files are post-process into organized analysis data.
Currently, three statistical parameters can be measured: BO population, BO coherence, NACME averaging.

1. BO population analysis based on the running state of each trajectory or based on density matrix of each trajectory

.. math::

   P_{i}(t) = \frac{N_{i}(t)}{N_{traj}} 

.. math::

   <\rho_{ii}(t)> = \frac{\sum_{I}^{N_{traj}} \rho_{ii}^{(I)}(t)}{N_{traj}}

2. BO coherence analysis based on density matrix of each trajectory

.. math::

   <\left\vert\rho_{ij}(t)\right\vert^{2}> = \frac{\sum_{I}^{N_{traj}} \rho_{ii}^{(I)}(t)\rho_{jj}^{(I)}(t)}{N_{traj}}

3. Averaging NACME, phase is ignored with absolute value

.. math::

   <\left\vert\sigma_{ij}(t)\right\vert> = \frac{\sum_{I}^{N_{traj}} \left\vert\sigma_{ij}^{(I)}(t)\right\vert}{N_{traj}}

Here, :math:`N_{traj}` represents total trajectory number, :math:`N_i(t)` represents the number of trajectories in :math:`i` state in time :math:`t`.
After running a script, analysis results can be found in current directory.

+------------------------+---------------------------------------------------------------+
| Option                 | Description                                                   |
+========================+===============================================================+
| **-n**, **-ntrajs**    | Total number of trajectories for analysis.                    |
|                        |                                                               |
+------------------------+---------------------------------------------------------------+
| **-s**, **-nsteps**    | Total number of steps for analysis.                           |
|                        |                                                               |
+------------------------+---------------------------------------------------------------+
| **-t**, **-nstates**   | Total number of states for analysis.                          |
|                        |                                                               |
+------------------------+---------------------------------------------------------------+
| **-h**                 | Call out help message.                                        |
|                        |                                                               |
+------------------------+---------------------------------------------------------------+

**Ex.** Statistical analysis on 100 trajectories, first 10 steps in 3 states.

.. code-block:: bash

   $ ls
   TRAJ_001/ TRAJ_002/ ... TRAJ_100/

   $ python3 statistical_analysis.py -n 100 -s 10 -t 3

After running a script, 'AVG_POPRUN', 'AVG_POPRHO', 'AVG_COHRHO', 'AVG_NACME' files will be generated in running directory.

.. code-block:: bash

   $ ls
   AVG_POPRUN AVG_POPRHO AVG_COHRHO AVG_NACME TRAJ_001/ TRAJ_002/ ... TRAJ_100/

Each generated file represents BO population based on running state, BO population based on density matrix, BO coherence based on density matrix, and averaged NACME, respectively.
Format of output files are following.

- AVG_POPRUN

.. code-block:: bash

     #   Running state based averaged BO population
     <MD_step>   <population_state_0>   <population_state_1>   <population_state_2>
     <MD_step>   <population_state_0>   <population_state_1>   <population_state_2>
     <MD_step>   <population_state_0>   <population_state_1>   <population_state_2>
     ...

- AVG_POPRHO

.. code-block:: bash

     #   Density matrix based averaged BO population
     <MD_step>   <population_state_0>   <population_state_1>   <population_state_2>
     <MD_step>   <population_state_0>   <population_state_1>   <population_state_2>
     <MD_step>   <population_state_0>   <population_state_1>   <population_state_2>
     ...

- AVG_COHRHO

.. code-block:: bash

     #   Averaged electronic coherence
     <MD_step>   <coherence_state_0>   <coherence_state_1>   <coherence_state_2>
     <MD_step>   <coherence_state_0>   <coherence_state_1>   <coherence_state_2>
     <MD_step>   <coherence_state_0>   <coherence_state_1>   <coherence_state_2>
     ...

- AVG_NACME

.. code-block:: bash

     #   Averaged Non-Adiabatic Coupling Matrix Eliments: off-diagonal
     <MD_step>   <NACME_(0, 1)>   <NACME_(0, 2)>   <NACME_(1, 2)>
     <MD_step>   <NACME_(0, 1)>   <NACME_(0, 2)>   <NACME_(1, 2)>
     <MD_step>   <NACME_(0, 1)>   <NACME_(0, 2)>   <NACME_(1, 2)>
     ...

motion_analysis.py
---------------------------
Python utility script for PyUNIxMD output analysis.
In this script, PyUNIxMD 'MOVIE.xyz' output file is post-process into given geometry criterion.
Currently, three geometrical parameters can be measured: bond length, bond angle, and dihedral angle.

1. In the bond length analysis, bond length between two given atoms will be calculated from given geometry information.
2. In the bond angle analysis, angle between three given atoms will be calculated. Here, second atom will be a vertex of angle. 
3. In the dihedral angle analysis, dihedral angle between four or six given atoms will be calculated. 
   In four atom case, dihedral angle between (1,2,3),(2,3,4) plane will be calculated and dihedreal axis will be atom2-atom3.
   In six atom case, dihedral angle between (1,2,3),(4,5,6) plane will be calculated and dihedral axis will be atom3-atom4.

After running a script, analysis results will be saved in md output directory in each trajectory ('TRAJ_(number)/md/').
If averaging argument is given, averaged results can be found in current directory.

+------------------------+-------------------------------------------------------------------+
| Option                 | Description                                                       |
+========================+===================================================================+
| **-n**, **--ntrajs**   | Total number of trajectories for analysis.                        |
|                        |                                                                   |
+------------------------+-------------------------------------------------------------------+
| **-s**, **--nsteps**   | Total number of steps.                                            |
|                        |                                                                   |
+------------------------+-------------------------------------------------------------------+
| **-b**, **--bond**     | Target bond length between two atoms.                             |
|                        | Two target atom numbers must be given as arguments.               |
+------------------------+-------------------------------------------------------------------+
| **-a**, **--angle**    | Target bond angle between three atoms.                            |
|                        | Three target atom numbers must be given as arguments.             |
+------------------------+-------------------------------------------------------------------+
| **-d**, **--dihedeal** | Target dihedral angle between four or six atoms.                  |
|                        | 4 or 6 target atom numbers must be given as arguments.            |
+------------------------+-------------------------------------------------------------------+
| **-m**, **--mean**     | Calculate averaged values through total trajectories.             |
|                        | This is optional argument.                                        |
+------------------------+-------------------------------------------------------------------+
| **-h**                 | Call out help message.                                            |
|                        |                                                                   |
+------------------------+-------------------------------------------------------------------+

**Ex.** Analyze bond length between atom 1 and 4 in 100 trajectories, first 10 steps. Averaging through total trajectories either.

.. code-block:: bash

   $ ls
   TRAJ_001/ TRAJ_002/ ... TRAJ_100/

   $ python3 motion_analysis.py -n 100 -s 10 -b 1 4 -m

After running a script, 'AVG_BOND' file will be generated in running directory and 'BOND' file will be generated in each trajectory directory.
'AVG_BOND' shows averaged bond length between two input atom number through total trajectories and 
'BOND' shows bond length between two input atom number in each trajectory.

.. code-block:: bash

   $ ls
   AVG_BOND TRAJ_001/ TRAJ_002/ ... TRAJ_100/

   $ ls TRAJ_001/md/
   BOND (other MD outputs) ...

Format of output files are following. Ouput file of bond angle and dihedral angle has same template.

- BOND

.. code-block:: bash

     #   Bond length between <atom_1> and <atom_2>
     <MD_step>  <Bond_length>
     <MD_step>  <Bond_length>
     <MD_step>  <Bond_length>
     ...

- AVG_BOND

.. code-block:: bash

     #   Averaged bond length between <atom_1> and <atom_2>
     <MD_step>  <Averaged_bond_length>
     <MD_step>  <Averaged_bond_length>
     <MD_step>  <Averaged_bond_length>
     ...
 
