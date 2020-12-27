MM_calculator
-------------------------------------------

The propagation of molecular dynamics needs the information about potential enrgy surfaces.
The object for MM calculator should be constructed in UNI-xMD to evaluate the energies
and forces. UNI-xMD interfaces with followings in current version.

.. toctree::
    :glob:
    :maxdepth: 1

    tinker

Only one method supports with MM interface of UNI-xMD in current version.
If you want to add other programs such as Gromacs, refer "Tinker" section.
The MM interfacings in UNI-xMD proceed with ``get_data`` method in each object.
Detailed explanation of ``get_data`` method is given in *modules* section.

**Ex.** Making a MM object with Tinker program

.. code-block:: python

   import mm

   mm = mm.Tinker(molecule=mol, scheme="subtractive", periodic=False, xyz_file="./tinker.xyz", \
       key_file="./tinker.key", embedding="electrostatic", vdw="lennardjones", \
       mm_path="/opt/tinker/bin/", nthreads=1, version=8.7)

.. note:: For QM/MM calculation, you should define **scheme** variable.
   In addition, same **embedding** variable of MM object should be used when one generate a QM object.


