QM_calculator
-------------------------------------------

The propagation of molecular dynamics needs the information about potential enrgy surfaces.
The object for QM calculator should be constructed in PyUNIxMD to evaluate the energies,
forces and nonadiabatic couplings. PyUNIxMD interfaces with the following in the current version.

.. toctree::
   :glob:
   :maxdepth: 1

   prog*
   model

Only some of methods are supported within PyUNIxMD in the current version.
If you want to add other methods such as CASPT2, please refer other interfaces.
The QM interfaces in PyUNIxMD proceed with **get_data** method in each object.
Detailed explanation of **get_data** method is given in *modules* section.

**Ex.** Making a QM object with CASSCF method of Molpro program

.. code-block:: python

   import qm

   qm = qm.molpro.CASSCF(molecule=mol, basis_set="6-31G", memory="500m", \
       active_elec=2, active_orb=2, cpscf_grad_tol=1E-7, \
       qm_path="/opt/molpro2015.1/bin/", nthreads=1, version=2015.1)


