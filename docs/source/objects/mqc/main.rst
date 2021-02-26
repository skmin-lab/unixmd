MQC
-------------------------------------------

Mixed quantum-classical (MQC) dynamics is general method for explaining the variation of molecule including
electronic state through time propagation. This can be exactly solved by time-dependent Schrodinger equation
for all particles, but this solution requires enormous cost for numerical calculation so it is restricted for
very small system. To overcome this limit, MQC tried to describe larger system by considering nuclei as classical 
particles which follow classical equation of motion.

PyUNIxMD mainly targeted on MQC, and whole dynamics implemented in current version of PyUNIxMD are subclass of
MQC class. In the MQC class, there are functions to update classical properties of nuclei.
MQC methods implemented in PyUNIxMD are listed in the following.

.. toctree::
    :glob:
    :maxdepth: 1

    *

Far more insights about treating MQC in terms of code structure, the overall modules are controlled in fundamental
input file 'run.py'. When users select their dynamics method, they have to make md object from the subclass of
:class:`MQC` class such as :class:`SH` (:class:`mqc.SH`), and a run method (``md.run``) to run that md object. In the md object, basic dynamics
parameters such as number of steps are given as arguments. Besides, the run method includes overall dynamics condition as arguments.

PyUNIxMD saves the objects for MQC and QM in a binary formatted file ('RESTART.bin') under **input_dir** directory using the 'pickle' package at every time step.
The 'RESTART.bin' file is overwritten at every successful MD step, therefore the file contains the information of a trajectory at the last successful MD step.
You can restart the dynamics simulation by reading the 'RESTART.bin' file using 'pickle' package.

Arguments for run method are listed below. The important point is that run method is included in each
md subclass of :class:`MQC`, not :class:`MQC` itself.

+-----------------------------+-------------------------------------------------+----------+
| Keywords                    | Work                                            | Default  |
+=============================+=================================================+==========+
| **qm**                      | QM object containing on-the-fly                 |          |
| (:class:`QM_calculator`)    | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **mm**                      | MM object containing MM calculation information | *None*   |
| (:class:`MM_calculator`)    |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+
| **input_dir**               | Name of directory where outputs to be saved     | *'./'*   |
| *(string)*                  |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+
| **save_qm_log**             | Logical for saving QM calculation log           | *False*  |
| *(boolean)*                 |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+
| **save_mm_log**             | Logical for saving MM calculation log           | *False*  |
| *(boolean)*                 |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+
| **save_scr**                | Logical for saving scratch directory            | *True*   |
| *(boolean)*                 |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+
| **restart**                 | Option for controlling dynamics restarting      | *None*   |
| *(string)*                  |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+

Further information of individual MD objects is listed in each section.

**Ex.** Making an MD object with FSSH method.

.. code-block:: python

   import mqc

   md = mqc.SH(molecule=mol, nsteps=1000, dt=0.125, istate=1, propagation="density")

   md.run(qm=qm, input_dir="./TRAJ.sh", save_scr=True, save_qm_log=False)

.. note:: Making molecule and QM objects are omitted in this sample code,
   but they must be declared to use run method in advance.

**Ex.** Restarting a dynamics simulation.

.. code-block:: python
   
   import pickle
   
   # Read the binary file and load the information 
   with open('./RESTART.bin', 'rb') as f_restart:
       restart = pickle.load(f_restart)
   qm = restart["qm"]
   md = restart["md"]
   
   # The dynamics is continued from the last step stored in the 'RESTART.bin' file.
   # The output files are appended with the restarted dynamics, therefore the MD output files should exist.
   md.run(qm=qm, restart='append')
   
   # The time is set to zero and restart a new dynamics from the information in the 'RESTART.bin' file.
   # The new output files are written.
   md.run(qm=qm, restart='write')


.. raw:: html

   <h2>Detailed description of arguments</h2>

- **input_dir** *(string)* - Default: *'./'*

  This argument designates directory for dynamics output. All subdirectories ('md', 'qm_log', and 'mm_log') for output files will be generated under **input_dir**.
  If the subdirectories are already present, old subdirectories will be renamed with '_old' and new subdirectories will be made.

\

- **save_qm_log** *(boolean)* - Default: *False*

  This argument determines saving QM calculation logs. Logs will be saved in '**input_dir**/qm_log/'.
 
\

- **save_mm_log** *(boolean)* - Default: *False*

  This argument determines saving MM calculation logs. Logs will be saved in '**input_dir**/mm_log/'.
  If **MM** = *None*, this argument will be ignored.

\

- **save_scr** *(boolean)* - Default: *True*

  This argument determines saving scratch output directory in '**input_dir**/md/scr_qm/' after QM calculation, and '**input_dir**/md/scr_mm/' after MM calculation.

\

- **restart** *(string)* - Default: *None*

  If this argument is not set (*None*), all objects will be initialized for a new dynamics simulation.
  Otherwise, 'RESTART.bin' is read and a dynamics simulation is restarted from the image read from the file without initialization of objects.

  + *'write'*: The dynamics simulation is restarted at resetting the time step to zero and writing new MD output files.
  + *'append'*: The dynamics simulation is continued from the image appending the existing MD output files.

\

