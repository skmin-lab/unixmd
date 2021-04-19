.. _Objects QM_calculator:

QM_calculator
-------------------------------------------

The propagation of molecular dynamics needs the information about potential enrgy surfaces.
The object for QM calculator should be constructed in PyUNIxMD to evaluate the energies,
forces and nonadiabatic couplings. PyUNIxMD interfaces with the following QM programs.

.. toctree::
   :glob:
   :maxdepth: 1

   prog*
   model

Not all QM methods of a QM program is available even though PyUNIxMD provides the interface for that program (for example, CASPT2 implemented in Molpro).
If there is no QM interface for some QM programs or QM methods in PyUNIxMD, you can make your own interface by refering to other interfaces.
The key method of each QM interface is ``get_data`` method which makes input files, executes calculations and extracts information at every MD step. 
Detailed description of ``get_data`` method is given in :ref:`Modules <Module QM>`.

In your running script, You need to access a specific QM interface where QM methods are provided in the form of Python classes.
The names of directories where the interface packages are and the classes of methods are tabulated below.

+-------------------+------------------------+----------------+----------------+
| QM programs       | Interface directories  | QM methods     | Class names    |
+===================+========================+================+================+
| COLUMBUS          | columbus               | (SA-)CASSCF    | CASSCF         |
|                   |                        |                |                |
|                   |                        | MRCI           | MRCI           |
+-------------------+------------------------+----------------+----------------+
| Molpro            | molpro                 | (SA-)CASSCF    | CASSCF         |
+-------------------+------------------------+----------------+----------------+
| Gaussian 09       | gaussian09             | (TD)DFT        | DFT            |
+-------------------+------------------------+----------------+----------------+
| Q-Chem            | qchem                  | (TD)DFT        | DFT            |
+-------------------+------------------------+----------------+----------------+
| TURBOMOLE         | turbomole              | (TD)DFT        | DFT            |
+-------------------+------------------------+----------------+----------------+
| TeraChem          | teraChem               | SI-SA-REKS     | SSR            |
+-------------------+------------------------+----------------+----------------+
| DFTB+             | dftbplus               | (TD)DFTB       | DFTB           |
|                   |                        |                |                |
|                   |                        | DFTB/SSR       | SSR            |
+-------------------+------------------------+----------------+----------------+

**Ex.** Making a QM object of the CASSCF method implemented in Molpro

.. code-block:: python

   import qm

   qm = qm.molpro.CASSCF(molecule=mol, basis_set="6-31G", memory="500m", \
       active_elec=2, active_orb=2, cpscf_grad_tol=1E-7, \
       qm_path="/opt/molpro2015.1/bin/", nthreads=1, version="2015.1")
