
Mixed quantum-classical (MQC) dynamics is general method for explaining the variation of molecule including
electronic state through time propagation. This can be exactly solved by time-dependent Schrodinger equation
for all particles, but this solution requires enormous cost for numerical calculation so it is restricted for
very small system. To overcome this limit, MQC tried to describe larger system by consider nuclear as classical 
particle which follows classical equation of motion.

UNI-xMD mainly targeted on MQC, and whole dynamics implemented in current version of UNI-xMD are subclass of
MQC class. In the MQC class, there are functions for update classical properties of nuclear.

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
| **qm**                      | qm object containing on-the-fly                 |          |
| (:class:`QM_calculator`)    | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **mm**                      | mm object containing MM                         | *None*   |
| (:class:`MM_calculator`)    | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **input_dir**               | location of input directory                     | *'./'*   |
| *(string)*                  | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **save_qm_log**             | logical for saving QM calculation log           | *False*  |
| *(boolean)*                 | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **save_mm_log**             | logical for saving MM calculation log           | *False*  |
| *(boolean)*                 | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **save_scr**                | logical for saving scratch directory            | *True*   |
| *(boolean)*                 | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **restart**                 | option for controlling dynamics restarting      | *None*   |
| *(string)*                  | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+

Further information of each individual md objects are listed in next section.

**Ex.** Making a md object with FSSH method

.. code-block:: python

   import mqc

   md = mqc.SH(molecule=mol, nsteps=1000, dt=0.125, istate=1, propagation="density")

   md.run(molecule=mol, qm=qm, input_dir="./TRAJ.sh", save_scr=True, save_qm_log=False)

   # molecule, thermostat and qm objects must be made to instatiate md object and use run method in advance.
   # Those lines are omitted in this sample code.

Detailed description of arguments
''''''''''''''''''''''''''''''''''''

- **input_dir** *(string)* - Default: *./*

 Direcrtory for dynamics output. All md output, saved log, ... etc will be saved in this directory.
 If the directory already presence, it will be removed and new directory will be made.

\

- **save_qm_log** *(boolean)* - Default: *False*

 Save QM calculation logs for passed timestep. Logs will be saved in 'input_dir/qm_log_dir'.
 
\
  
- **save_mm_log** *(boolean)* - Default: *False*

 Save MM calculation logs for passed timestep. Logs will be saved in 'input_dir/mm_log_dir'.
 If there is no mm data, this argument will be ignored.

\
   
- **save_scr** *(boolean)* - Default: *True*

 Save scratch output directory in 'input_dir/md/scr' after qm calculation.

\

- **restart** *(string)* - Default: *None*

 Writting options for restarting dynamics from halted trajectory.

+ write: Write a new output, starting from halted timestep.
+ append: Write a output continually starting from a halted timestep.

\

