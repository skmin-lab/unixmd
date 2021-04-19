MM_calculator
-------------------------------------------

The propagation of molecular dynamics needs the information about potential enrgy surfaces.
The object for MM calculator should be constructed in PyUNIxMD to evaluate the energies
and the forces in QM/MM calculation. PyUNIxMD interfaces with the following in the current version.

.. toctree::
    :glob:
    :maxdepth: 1

    prog_tinker

Only one MM program is supported in the current version of PyUNIxMD.
If you want to add an interface of other programs such as Gromacs, you can refer to the Tinker interface.
The key method of the MM interface is ``get_data`` method which makes input files, executes calculations and extracts information at every MD step. 
Detailed description of the ``get_data`` method is given in :ref:`Modules <Module MM>`.

**Ex.** Making an MM object with Tinker

.. code-block:: python

   import mm

   mm = mm.Tinker(molecule=mol, scheme="subtractive", l_periodic=False, xyz_file="./tinker.xyz", \
       key_file="./tinker.key", embedding="electrostatic", vdw="lennardjones", \
       mm_path="/opt/tinker/bin/", nthreads=1, version="8.7")

.. note:: For QM/MM calculation, you must define **scheme** variable.
   In addition, the same **embedding** variable of the MM object should be used when one generate a QM object.


