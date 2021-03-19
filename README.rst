*****************************************************************
PyUNIxMD: Python-based UNIversal eXcited state Molecular Dynamics
*****************************************************************
PyUNIxMD is an object-oriented Python program for molecular dynamics simulations involving multiple electronic states.
It is mainly for studying the nonadiabatic dynamics of excited molecules.

Requirements
============
* Python 3.6 or later
* Numpy
* Cython https://cython.org
        
You can easily install the latest Numpy and Cython via Python's pip command.      
command::        
        
  $ pip install --upgrade numpy Cython
    
Build
=====
You can build PyUNIxMD by the following command.
command:: 

  $ python3 setup.py build_ext -b ./src/build

Test
====
Without the aid of external QM programs, You can test the PyUNIxMD package with model systems.
The corresponding examples are:

* $PYUNIXMDHOME/examples/qm/Eh-Shin_Metiu

* $PYUNIXMDHOME/examples/qm/SHXF-SAC

* $PYUNIXMDHOME/examples/qm/SHXF-Shin_Metiu

In each directory, you can find the running script named run.py

Before running test jobs, add the path of the PyUNIxMD package to your Python path.
command::

  $ export PYTHONPATH=$PYUNIXMDHOME/src:$PYTHONPATH

Then execute run.py as follows.
command::

  $ python3 run.py >& log

Documentation
=============
If you have Sphinx, you can locally build the manual of PyUNIxMD by the following command.
command::

  $ cd docs
  $ make html

License
=======