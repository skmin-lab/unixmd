.. _Objects QED_calculator:

QED_calculator
-------------------------------------------

The propagation of polariton dynamics needs the information about polaritonic potential enrgy surfaces.
The object for QED calculator should be constructed in PyUNIxMD to evaluate the energies,
forces, and nonadiabatic coupling vectors with polaritonic states. PyUNIxMD interfaces with the following in the current version.

.. toctree::
   :glob:
   :maxdepth: 1

   jaynes_cummings

Only one model (i.e., Jaynes-Cummings model) is supported in the current version of PyUNIxMD.
If you want to add an interface of other models such as Tavis-Cummings model, you can refer to the Jaynes-Cummings interface.
The key method of the QED interface is ``get_data`` method which constructs model Hamiltonian and computes required information at every MD step. 
Detailed description of the ``get_data`` method is given in :ref:`QED <Module QED>`.

In the current version of PyUNIxMD, the compatibility of QED methods with QM methods is tabulated below.

+-------------------+-----------------+----------------+
| QED methods       | QM programs     | QM methods     |
+===================+=================+================+
| Jaynes_Cummings   | DFTB+           | SI-SA-REKS     |
+-------------------+-----------------+----------------+

**Ex.** Making a QED object with Jaynes-Cummings model

.. code-block:: python

   import qed

   qed = qed.Jaynes_Cummings(polariton=pol, coupling_strength=0.001, force_level="full", l_check_crossing=True)

.. note:: For polariton dynamics, you must define a Polariton object including information about cavity photons.

