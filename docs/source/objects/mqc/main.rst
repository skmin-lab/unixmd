MQC
-------------------------------------------

Mixed quantum-classical (MQC) dynamics is general method for explaining the variation of molecule including
electronic state through time propagation. This can be exactly solved by time-dependent Schrodinger equation
for all particles, but this solution requires enormous cost for numerical calculation so it is restricted for
very small system. To overcome this limit, MQC tried to describe larger system by consider nuclear as classical 
particle which follows classical equation of motion.

PyUNIxMD mainly targeted on MQC, and whole dynamics implemented in current version of PyUNIxMD are subclass of
MQC class. In the MQC class, there are functions for update classical properties of nuclear.
MQC methods implemented in UNI-xMD are listed in the following.

.. toctree::
    :glob:
    :maxdepth: 1

    *

Far more insights about treating MQC in terms of code structure, the overall modules are controlled in fundamental
input file run.py. When user select their dynamics method, they have to make md object from the subclass of
:class:`MQC` class such as :class:`SH` (:class:`mqc.SH`), and a run method (``md.run``) to run that md object. In the md object, basic dynamics
parameters such as number of steps are given as arguments. Besides, run methods includes overall dynamics condition
as arguments.

Arguments for ``run`` method are listed below. The important point is that ``run`` method is included in each
md subclasses of :class:`MQC`, not :class:`MQC` itself.

+-----------------------------+-------------------------------------------------+----------+
| Keywords                    | Work                                            | Default  |
+=============================+=================================================+==========+
| **qm**                      | QM object containing on-the-fly                 |          |
| (:class:`QM_calculator`)    | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **mm**                      | MM object containing MM                         | *None*   |
| (:class:`MM_calculator`)    | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **input_dir**               | Location of input directory                     | *'./'*   |
| *(string)*                  | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **save_qm_log**             | Logical for saving QM calculation log           | *False*  |
| *(boolean)*                 | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **save_mm_log**             | Logical for saving MM calculation log           | *False*  |
| *(boolean)*                 | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **save_scr**                | Logical for saving scratch directory            | *True*   |
| *(boolean)*                 | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **restart**                 | Option for controlling dynamics restarting      | *None*   |
| *(string)*                  | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+

Further information of each individual MD objects are listed in next section.

**Ex.** Making a MD object with FSSH method
        Making molecule and QM objects are omitted in this sample code, but they must be declared to use run method in advance.

.. code-block:: python

   import mqc

   md = mqc.SH(molecule=mol, nsteps=1000, dt=0.125, istate=1, propagation="density")

   md.run(qm=qm, input_dir="./TRAJ.sh", save_scr=True, save_qm_log=False)



Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **input_dir** *(string)* - Default: *'./'*

  Direcrtory for dynamics output. All md output, saved log, ... etc will be saved in this directory.
  If the directory already present, it will be removed and new directory will be made.

\

- **save_qm_log** *(boolean)* - Default: *False*

  Save QM calculation logs for passed timestep. Logs will be saved in **'input_dir/qm_log_dir'**.
 
\
  
- **save_mm_log** *(boolean)* - Default: *False*

  Save MM calculation logs for passed timestep. Logs will be saved in **'input_dir/mm_log_dir'**.
  If there is no MM data, this argument will be ignored.

\
   
- **save_scr** *(boolean)* - Default: *True*

  Save scratch output directory in **'input_dir'**/md/scr_qm after QM calculation, and **'input_dir'**/md/scr_mm after MM calculation.

\

- **restart** *(string)* - Default: *None*

  Writting options for restarting dynamics from halted trajectory.

  + *'write'*: Write a new output, starting from halted timestep.
  + *'append'*: Write a output continually starting from a halted timestep.

\

