MQC
-------------------------------------------

Mixed quantum-classical (MQC) dynamics is one of theoretical tools to simulate nonadiabatic processes in molecular or extended systems. This method is characterized by propagating the electrons through quantum mechanics but nuclear dynamics through classical trajectories, to overcome the computational cost of calculating the correlated systems.

PyUNIxMD provides a variety of MQC methods and they are implemented as subclasses of
MQC class. In the MQC class, the common properties and methods are defined such as functions to update classical properties of nuclei.
MQC methods implemented in PyUNIxMD are listed in the following.

.. toctree::
    :glob:
    :maxdepth: 1

    *

In the running script 'run.py', you need to specify the dynamics you want to run by making an md object from the subclasses of
:class:`MQC` class such as :class:`SH` (:class:`mqc.SH`). After setting the dynamics method, you use run method (``md.run``) to perform the dynamics.
When making an md object, basic dynamics parameters such as the number of steps are needed to be specified. Also, the run method takes various parameters such as a QM object, thermostat, etc.

PyUNIxMD saves the objects for MQC and QM in a binary formatted file ('RESTART.bin') under **input_dir** directory using the 'pickle' package at every time step.
The 'RESTART.bin' file is overwritten at every successful MD step, therefore the file contains the information of a trajectory at the last successful MD step.
You can restart the dynamics simulation by reading the 'RESTART.bin' file using 'pickle' package.

Parameters for the run method are listed below. The important point is that the run method is included in each
md subclass of :class:`MQC`, not :class:`MQC` itself.

+-----------------------------+-------------------------------------------------+----------+
| Parameters                  | Work                                            | Default  |
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

   <h2>Detailed description of parameters</h2>

- **input_dir** *(string)* - Default: *'./'*

  This parameter designates the directory for dynamics output. All subdirectories ('md', 'qm_log', and 'mm_log') for output files will be generated under **input_dir**.
  If the subdirectories are already present, old subdirectories will be renamed with '_old' and new subdirectories will be made.

\

- **save_qm_log** *(boolean)* - Default: *False*

  This parameter determines whether to save QM calculation logs. Logs will be saved in '**input_dir**/qm_log/'.
 
\

- **save_mm_log** *(boolean)* - Default: *False*

  This parameter determines whether to save MM calculation logs. Logs will be saved in '**input_dir**/mm_log/'.
  If **MM** = *None*, this parameter will be ignored.

\

- **save_scr** *(boolean)* - Default: *True*

  This parameter determines whether to save scratch output directory in '**input_dir**/md/scr_qm/' after QM calculation, and '**input_dir**/md/scr_mm/' after MM calculation.

\

- **restart** *(string)* - Default: *None*

  If this parameter is not set (*None*), all objects will be initialized for a new dynamics simulation.
  Otherwise, 'RESTART.bin' is read and a dynamics simulation is restarted from the image read from the file without initialization of objects.

  + *'write'*: The dynamics simulation is restarted at resetting the time step to zero and writing new MD output files.
  + *'append'*: The dynamics simulation is continued from the image appending the existing MD output files.

\

