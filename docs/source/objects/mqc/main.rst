MQC
-------------------------------------------

Mixed quantum-classical (MQC) dynamics is general method for explaining the variation of molecule including
electronic state through time propagation. This can be exactly solved by time-dependent Schrodinger equation
for all particles, but this solution requires enormous cost for numerical calculation so it is restricted for
very small system. MQC tried to describe larger system by consider nuclear as classical particle which follows
classical equation of motion.

UNI-xMD mainly targeted on MQC, and whole dynamics implemented in current version of UNI-xMD are subclass of
MQC class. In the MQC class, there are functions for update classical properties of nuclear. MQC methods implemented 
in UNI-xMD are listed in the following.

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

   md = mqc.SH(molecule=mol, thermostat=bathT, nsteps=1000, dt=0.125, istate=1, propagation="density")

   md.run(qm=qm, input_dir="./TRAJ.sh", save_scr=True, save_qm_log=False)

   # molecule, thermostat and qm objects must be made to instatiate md object and use run method in advance.
   # Those lines are omitted in this sample code.

