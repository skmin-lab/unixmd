===========================
PyUNIxMD Overview
===========================

Features
---------------------------
The features of PyUNIxMD are as follows.

- Conventional (non)adiabatic dynamics

  - Born-Oppenheimer molecular dynamics (BOMD)

  -  Ehrenfest dynamics

  -  Fewest switches surface hopping (FSSH) dynamics with ad hoc decoherence corrections

- Decoherence based on exact factorization

  -  Decoherence induced surface hopping based on exact factorization (DISH-XF)
  -  Coupled-trajectory mixed quantum/classical (CTMQC) method

- Numerical calculation of nonadiabatic couplings
- Accessible interface to external QM programs and built-in model Hamiltonians

  -  COLUMBUS: SA-CASSCF
  -  Molpro: SA-CASSCF
  -  Gaussian 09: TDDFT
  -  Q-Chem: TDDFT
  -  TURBOMOLE: TDDFT
  -  TeraChem: SI-SA-REKS (SSR)
  -  DFTB+: TDDFTB, DFTB/SSR
  -  Model Hamiltonians: Tully, Shin-Metiu

- QM/MM functionalities
- Utility scripts in Python

Authors
---------------------------
The current version of PyUNIxMD has been developed by Seung Kyu Min, In Seong Lee, Jong-Kwon Ha, Daeho Han, Kicheol Kim, Tae In Kim, Sung Wook Moon in the Theoretical/Computational Chemistry Group for Excited State Phenomena of Ulsan National Institute of Science and Technology (UNIST). 

..
  Acknowledgement
  ---------------------------
  This is acknowledgement.


Program Structure
---------------------------
The overall code structure is displayed in the next figure.

.. image:: ./pyunixmd_structure.png
   :width: 400pt

PyUNIxMD is an object-oriented program consisting of
several key classes closely connected with each other:

- :class:`Molecule` defines a target system. A molecule object contains information of the electronic states as well as the geometry.

- :class:`MQC` has information about molecular dynamics. Each nonadiabatic dynamics method (Ehrenfest, surface hopping, etc.) comprises its subclasses. 

- :class:`QM_calculator` interfaces several QM programs (Molpro, Gaussian 09, DFTB+, etc.) and methodologies to perform electronic structure calculations.

- :class:`MM_calculator` enables QM/MM calculations using external software such as Tinker.

- :class:`Thermostat` controls temperature of a target system.

PyUNIxMD takes advantage of the inheritance feature to organize functionalities and simplify the codes by sharing the common parameters and methods.

For detailed information of each class, see Section ?. 

