*****************************************************************
PyUNIxMD: Python-based UNIversal eXcited state Molecular Dynamics
*****************************************************************

.. image:: image/logo.svg
      :width: 600pt
      :align: center
      
PyUNIxMD is an object-oriented Python program for molecular dynamics simulations involving multiple electronic states.
It is mainly for studying the nonadiabatic dynamics of excited molecules.

Citation
========

When you publish a work using any part of PyUNIxMD code, please cite the following publication:

* PyUNIxMD: A Python-based excited state molecular dynamics package. *J. Comp. Chem.* **2021**, DOI: `10.1002/jcc.26711 <https://doi.org/10.1002/jcc.26711>`_

* Coupled- and independent-trajectory approaches based on the exact factorization using the PyUNIxMD package. *Top. Cur. Chem.* **2022**, DOI: `10.1007/s41061-021-00361-7 <https://doi.org/10.1007/s41061-021-00361-7>`_

Requirements
============
* Python 3.6 or later
* Numpy 1.20.0
* Cython https://cython.org
* BLAS/LAPACK libraries or Math Kernel Library

You can easily install the latest Numpy and Cython via Python's pip command.

::
        
  $ pip install --upgrade numpy Cython
    
Build
=====
You can build PyUNIxMD by the following command.

:: 

  $ python3 setup.py build_ext -b ./src/build

If you want to run polariton dynamics, first build PyUNIxMD as above
and then set 'do_qed = True' in setup.py and compile using the following command:

:: 

  $ python3 setup.py build_qed_ext -b ./src/build_qed

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

