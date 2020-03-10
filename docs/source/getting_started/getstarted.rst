==========================
Build
==========================

UNI-xMD is python based program with a little C code for time-consuming electronic propagation part interfaced via Cython,
therefore compiliation is needed.


Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^
Python 3.5 (or newer)

Numpy

Cython https://cython.org


If you don't have numpy or Cython, you can install them using pip command.

.. code-block:: bash
   
   $ pip install --upgade numpy Cython

Compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can compile electronic propagation routine by typing following command in root directory of the program which contains setup.py file.

.. code-block:: bash

   $ python3 setup.py build_ext -b ./src/mqc/el_prop

================================
Code Flow and Program Structure
================================
Simple intro..?

Work flow
^^^^^^^^^^^^^^^^^^^^^^^^^^
Diagram 

Program structure
^^^^^^^^^^^^^^^^^^^^^^^^^^
Inheritance,...

==========================
Quick Start
==========================

Once installation is complete, the user is now able to run UNI-xMD.
A straightfoward way to perform a UNI-xMD calculation is as follows:

-  First of all, you should add the corresponding directories to your python path.

.. code-block:: bash
  
    export PYTHONPATH=:[UNIxMDHOME]/src:$PATHONPATH
 
- You should import some library in running script.

.. code-block:: python

    from molecule import Molecule
    import bo
    import mqc
    from thermostat import *
    from misc import data

- To run UNI-xMD, you should create several objects, which are ``molecule``, ``bo``, ``md`` and ``thermostat``, in your running script. The important thing is that molecule object is created in the first place. 

Define system
^^^^^^^^^^^^^^^^^^^^^^^^^^
Next, You define molecular infomation contains atom number and geometry. 

.. code-block:: python

    geom = """
    6
    c2h4
    C       -0.69284974      0.00000000      0.00000000   0.0 0.0 0.0
    C        0.66284974      0.00000000      0.00000000   0.0 0.0 0.0
    H       -1.03381877      0.72280451     -0.32280655   0.0 0.0 0.0
    H       -1.23381879     -0.72280451     -0.34280655   0.0 0.0 0.0
    H        1.03381877      0.75280451     -0.32280655   0.0 0.0 0.0
    H        1.03381877     -0.72280451     -0.31280655   0.0 0.0 0.0
    """

This infomation is used to create molecule object.

.. code-block:: python
    
    mol = Molecule(geometry=geom, nstates=2)

Here, we use two adiabatic potential energy surfaces in our tutorial.

.. note:: molecule object should be already created before creating another objects such as bo method, md method and thermostat type.

Define BO method
^^^^^^^^^^^^^^^^^^^^^^^^^^
The user determine on-the-fly calculation method to get energy, force and non-adiabatic coupling vector.

.. code-block:: python
   
    qm = bo.dftbplus.SSR(molecule=mol, scc_tol)

when on-the-fly calculation is performed, We use blablalba.

Define MD method
^^^^^^^^^^^^^^^^^^^^^^^^^^
In this part, the user determine method of MD. Ehrenfest dynamics, FSSH, BOMD and SHXF are possible in UNI-xMD.

.. code-block:: python
   
    md = mqc.SHXF(molecule=mol, nsteps=1000, dt=0.125, istate=1, propagation="coefficient", \
    wsigma=0.1, threshold=0.01)


Define Thermostat method
^^^^^^^^^^^^^^^^^^^^^^^^^^
When the user choose thermostat type, there are three types of thermostat, that () states in detail.
We use TYPE which blalblabla.

.. code-block:: python
   
    bathT = rescale1(temperature=300.0, nrescale=20)


Run MD
^^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, 

.. code-block:: python
   
    md.run(~~~)

Check output
^^^^^^^^^^^^^^^^^^^^^^^^^^
UNI-xMD provieds various output file.

- blabla

- blabla


