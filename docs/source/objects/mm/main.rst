MM_calculator
-------------------------------------------

The propagation of molecular dynamics needs the information about potential enrgy surfaces.
The object for MM calculator should be constructed in PyUNIxMD to evaluate the energies
and forces in QM/MM calculation. PyUNIxMD interfaces with followings in current version.

.. toctree::
    :glob:
    :maxdepth: 1

    prog_tinker

Only one program supports with MM interface of PyUNIxMD in current version.
If you want to add other programs such as Gromacs, refer "Tinker" section.
The MM interfacings in PyUNIxMD proceed with ``get_data`` method in each object.
Detailed explanation of ``get_data`` method is given in *modules* section.

**Ex.** Making a MM object with Tinker program

.. code-block:: python

   import mm

   mm = mm.Tinker(molecule=mol, scheme="subtractive", periodic=False, xyz_file="./tinker.xyz", \
       key_file="./tinker.key", embedding="electrostatic", vdw="lennardjones", \
       mm_path="/opt/tinker/bin/", nthreads=1, version="8.7")

.. note:: For QM/MM calculation, you must define **scheme** variable.
   In addition, same **embedding** variable of MM object should be used when one generate a QM object.


