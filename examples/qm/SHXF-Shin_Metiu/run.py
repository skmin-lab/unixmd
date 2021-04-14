from molecule import Molecule
import qm, mqc
from misc import data

# For model calculation, set mass of a fictitious particle of the model.
data[f"X1"] = 1836 # au


# Random seed is fixed for test purpose to check if the same result is reproduced.
import random
random.seed(10)


# Define the target system.
geom = """
1
Shin-Metiu model
X1       -4.0     0.0
"""
mol = Molecule(geometry=geom, ndim=1, nstates=2, ndof=1, unit_pos='au', l_model=True)


# Set QM method.
qm = qm.model.Shin_Metiu(molecule=mol)


# Determine MD method.
md = mqc.SHXF(molecule=mol, nsteps=2890, nesteps=1, dt=0.5, unit_dt='au', \
     sigma=0.1, istate=1, elec_object="density")


# Execute the simulation.
md.run(qm=qm)
