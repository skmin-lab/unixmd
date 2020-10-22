from molecule import Molecule
import qm, mqc
from misc import data, fs_to_au
import math

data[f"X1"] = 2000

geom = """
1
spin-boson model
X1  -18.56599360        0.01030671
"""

mol = Molecule(geometry=geom, nsp=1, nstates=2, dof=1, unit_pos='au', unit_vel='au', model=True)

qm = qm.model.SAC(molecule=mol)

t = 0.1 / fs_to_au
md = mqc.SHXF(molecule=mol, nsteps=23000, nesteps=1, dt=t, wsigma=0.05, istate=0, propagation="density")

md.run(molecule=mol, qm=qm, input_dir=f"./", save_scr=True, save_QMlog=False)
