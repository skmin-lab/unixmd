*****************************************************************
PyUNIxMD: Python-based UNIversal eXcited state Molecular Dynamics
*****************************************************************

![Alt text](./controllers_brief.svg)
<img src="./image/logo.svg">

PyUNIxMD is an object-oriented Python program for molecular dynamics simulations involving multiple electronic states.
It is mainly for studying the nonadiabatic dynamics of excited molecules.

Requirements
============
* Python 3.6 or later
* Numpy
* Cython https://cython.org
        
You can easily install the latest Numpy and Cython via Python's pip command.

::
        
  $ pip install --upgrade numpy Cython
    
Build
=====
You can build PyUNIxMD by the following command.

:: 

  $ python3 setup.py build_ext -b ./src/build

Test
====
Without the aid of external QM programs, You can test the PyUNIxMD package with model systems.
The corresponding examples are:

* $PYUNIXMDHOME/examples/qm/Eh-Shin_Metiu

* $PYUNIXMDHOME/examples/qm/SHXF-SAC

* $PYUNIXMDHOME/examples/qm/SHXF-Shin_Metiu

$PYUNIXMDHOME is the top-level directory where this file belongs.

In each directory, you can find the running script named run.py.

Before running test jobs, add the path of the PyUNIxMD package to your Python path.

::

  $ export PYTHONPATH=$PYUNIXMDHOME/src:$PYTHONPATH

Then execute run.py as follows.

::

  $ python3 run.py >& log

Utility scripts
===============
PyUNIxMD provides other Python scripts to analyze results of dynamics calculations.
To use the scripts, you need to add the path of the scripts.

::

  $ export PYTHONPATH=$PYUNIXMDHOME/util:$PYTHONPATH

Documentation
=============
If you have Sphinx, you can locally build the manual of PyUNIxMD by the following command.

::

  $ cd docs
  $ make html

License
=======
