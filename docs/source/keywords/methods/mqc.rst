
Mixed quantum-classical (MQC) dynamics is general method for explaining the variation of molecule including
electronic state through time propagation. This can be exactly solved by time-dependent Schrodinger equation
for all particles, but this solution requires enormous cost for numerical calculation so it is restricted for
very small system. MQC tried to describe larger system by consider nuclear as classical particle which follows
classical equation of motion.

UNI-xMD mainly targeted on MQC, and whole dynamics implemented in current version of UNI-xMD are subclass of
MQC class. In the MQC class, there are functions for update classical properties of nuclear.

Far more insights about treating MQC in terms of code structure, the overall modules are controlled in fundamental
input file run.py. When user select their dynamics method, they have to make md object from the subclass of
``MQC`` class such as ``SH`` (mqc.SH), and a run method (md.run) to run that md object. In the md object, basic dynamics
parameters such as number of steps are given as arguments. Besides, run methods includes overall dynamics condition
as arguments.

Arguments for ``run`` method are listed below. The important point is that ``run`` method is included in each
md subclasses of ``MQC``, not ``MQC`` itself.

+-----------------+-------------------------------------------------+-----------+
| Keywords        | Work                                            | Default   |
+=================+=================================================+===========+
| ``qm``          | qm object containing on-the-fly                 |           |
|                 | calculation information                         |           |
+-----------------+-------------------------------------------------+-----------+
| ``mm``          | mm object containing MM                         | ``None``  |
|                 | calculation information                         |           |
+-----------------+-------------------------------------------------+-----------+
| ``input_dir``   | location of input directory                     | ``./``    |
+-----------------+-------------------------------------------------+-----------+
| ``save_qm_log`` | logical for saving QM calculation log           | ``False`` |
+-----------------+-------------------------------------------------+-----------+
| ``save_mm_log`` | logical for saving MM calculation log           | ``False`` |
+-----------------+-------------------------------------------------+-----------+
| ``save_scr``    | logical for saving scratch directory            | ``True``  |
+-----------------+-------------------------------------------------+-----------+
| ``restart``     | option for controlling dynamics restarting      | ``None``  |
+-----------------+-------------------------------------------------+-----------+

Further information of each individual md objects are listed in next section.

**Ex.** Making a md object with FSSH method

.. code-block:: python

   import mqc

   md = mqc.SH(molecule=mol, thermostat=bathT, nsteps=1000, dt=0.125, istate=1, propagation="density")

   md.run(qm=qm, input_dir="./TRAJ.sh", save_scr=True, save_qm_log=False)

   # molecule, thermostat and qm objects must be made to instatiate md object and use run method in advance.
   # Those lines are omitted in this sample code.

