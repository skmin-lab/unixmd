.. _Objects MQC_QED:

MQC_QED
-------------------------------------------

For MQC dynamcis in cQED, the total wavefunction can be expanded as a superposition
of the polaritonic states instead of the electronic states. However, it is numerically
difficult to compute the derivative of the unitary matrix since the unitary matrix is
not uniquely defined. Instead, the electrons are propagated in terms of uncoupled basis (:math:`D^{(I)}_{\beta}(t)`)
via unitary transformation from the polaritonic state (:math:`C^{(I)}_{J}(t)`). Other properties such as hopping
in FSSH method are similarly decided. The unitary trasformation between the polaritonic
and uncoupled states is given by

.. math::

   D^{(I)}_{\beta}(t) = \sum_{J} U^J_{\beta} C^{(I)}_{J}(t).

PyUNIxMD provides a variety of MQC_QED methods:

.. toctree::
    :glob:
    :maxdepth: 1

    *

In your running script, you need to specify the MQC_QED method you want to use by making an object of it.
In PyUNIxMD, MQC_QED methods are provided in the form of Python classes under :class:`MQC_QED` class.
The class names are tabulated below.

+-------------------+----------------+
| MQC_QED methods   | Class names    |
+===================+================+
| BOMD              | BOMD           |
+-------------------+----------------+
| FSSH              | SH             |
+-------------------+----------------+
| DISH-XF           | SHXF           |
+-------------------+----------------+

For example, a MD object of the FSSH method can be created as follows.

.. code-block:: python

   import mqc_qed
 
   md = mqc_qed.SH(...)

The parameters for the initialization are different for each MQC_QED method. For the detailed list of these parameters, see the subsections.

.. note:: In initialization of MQC_QED class, you must define a Polariton object including information about cavity photons
          instead of a Molecule object.

All classes specifying a MQC_QED method have their own ``run`` method. The ``run`` method is used to perform the dynamics at the end of your running script.
Parameters for the run method are listed below.

Plus, The ``run`` method deals with restart options of dynamics calculations. PyUNIxMD saves the objects for MQC_QED and QM in a binary formatted file ('RESTART.bin') under **output_dir** directory using the 'pickle' package at every time step.
The 'RESTART.bin' file is overwritten at every successful MD step, therefore the file contains the information of a trajectory at the last successful MD step.
You can restart the dynamics simulation by reading the 'RESTART.bin' file using 'pickle' package.

+-----------------------------+-------------------------------------------------+----------+
| Parameters                  | Work                                            | Default  |
+=============================+=================================================+==========+
| **qed**                     | QED object containing cavity-molecule           |          |
| (:class:`QED_calculator`)   | interaction                                     |          |
+-----------------------------+-------------------------------------------------+----------+
| **qm**                      | QM object containing on-the-fly                 |          |
| (:class:`QM_calculator`)    | calculation information                         |          |
+-----------------------------+-------------------------------------------------+----------+
| **mm**                      | MM object containing MM calculation information | *None*   |
| (:class:`MM_calculator`)    |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+
| **output_dir**              | Name of directory where outputs to be saved     | *'./'*   |
| *(string)*                  |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+
| **l_save_qed_log**          | Logical for saving QED calculation log          | *False*  |
| *(boolean)*                 |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+
| **l_save_qm_log**           | Logical for saving QM calculation log           | *False*  |
| *(boolean)*                 |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+
| **l_save_mm_log**           | Logical for saving MM calculation log           | *False*  |
| *(boolean)*                 |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+
| **l_save_scr**              | Logical for saving scratch directory            | *True*   |
| *(boolean)*                 |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+
| **restart**                 | Option for controlling dynamics restarting      | *None*   |
| *(string)*                  |                                                 |          |
+-----------------------------+-------------------------------------------------+----------+

**Ex.** Running FSSH dynamics with a MD object of the FSSH method.

.. code-block:: python

   import mqc_qed

   md = mqc_qed.SH(polariton=pol, nsteps=1000, dt=0.125, istate=2, elec_object="coefficient")

   md.run(qed=qed, qm=qm, output_dir="./TRAJ.sh", l_save_scr=True, l_save_qed_log=False, l_save_qm_log=False)

.. note:: Making polariton, QM, and QED objects are omitted in this sample code,
   but they must be declared to use a run method in advance.

**Ex.** Restarting a dynamics simulation.

.. code-block:: python
   
   import pickle
   
   # Read the binary file and load the information 
   with open('./RESTART.bin', 'rb') as f_restart:
       restart = pickle.load(f_restart)
   qed = restart["qed"]
   qm = restart["qm"]
   md = restart["md"]
   
   # The dynamics is continued from the last step stored in the 'RESTART.bin' file.
   # The output files are appended with the restarted dynamics, therefore the MD output files should exist.
   md.run(qed=qed, qm=qm, restart='append')
   
   # The time is set to zero and restart a new dynamics from the information in the 'RESTART.bin' file.
   # The new output files are written.
   md.run(qed=qed, qm=qm, restart='write')


.. raw:: html

   <h2>Detailed description of parameters</h2>

- **output_dir** *(string)* - Default: *'./'*

  This parameter designates the directory for dynamics output. All subdirectories ('md/', 'qed_log/', 'qm_log/', and 'mm_log/') for output files will be generated under **output_dir**.
  If the subdirectories are already present, old subdirectories will be renamed with '_old' and new subdirectories will be made.

\

- **l_save_qed_log** *(boolean)* - Default: *False*

  This parameter determines whether to save QED calculation logs. Logs will be saved in '**output_dir**/qed_log/'.
 
\

- **l_save_qm_log** *(boolean)* - Default: *False*

  This parameter determines whether to save QM calculation logs. Logs will be saved in '**output_dir**/qm_log/'.
 
\

- **l_save_mm_log** *(boolean)* - Default: *False*

  This parameter determines whether to save MM calculation logs. Logs will be saved in '**output_dir**/mm_log/'.
  If **MM** = *None*, this parameter will be ignored.

\

- **l_save_scr** *(boolean)* - Default: *True*

  This parameter determines whether to save scratch output directory in '**output_dir**/md/scr_qed/' after QED calculation, '**output_dir**/md/scr_qm/' after QM calculation, and '**output_dir**/md/scr_mm/' after MM calculation.

\

- **restart** *(string)* - Default: *None*

  If this parameter is not set (*None*), all objects will be initialized for a new dynamics simulation.
  Otherwise, 'RESTART.bin' is read and a dynamics simulation is restarted from the image read from the file without initialization of objects.

  + *'write'*: The dynamics simulation is restarted at resetting the time step to zero and writing new MD output files.
  + *'append'*: The dynamics simulation is continued from the image appending the existing MD output files.


