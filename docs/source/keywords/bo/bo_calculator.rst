
The propagation of molecular dynamics needs the information about BO potential enrgy surfaces.
The object for BO calculator should be constructed in UNI-xMD to evaluate the energies,
forces and nonadiabatic couplings. UNI-xMD interfaces with followings in current version.

- Columbus, DFTB+, Gaussian09, Molpro, TeraChem, Turbomole, Model system

Only some methods supports with BO interface of UNI-xMD in current version.
If you want to add other methods such as CASPT2, see "Do It yourself" section.
The BO interfacings in UNI-xMD proceed with **get_bo** method in each object.
Detailed explanation of **get_bo** method is given in *modules* section.

**Ex.** Making a BO object with CASSCF method of Molpro program

.. code-block:: python

   import bo

   qm = bo.molpro.CASSCF(molecule=mol, basis_set="6-31G", memory="500m", \
        active_elec=2, active_orb=2, cpscf_grad_tol=1E-7, \
        qm_path="/opt/molpro2015.1/bin/", nthreads=1, version=2015.1)

