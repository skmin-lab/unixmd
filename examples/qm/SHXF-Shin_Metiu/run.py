from molecule import Molecule
import qm, mqc
from misc import data

data[f"X1"] = 1836 # au

geom = """
1
Shin-Metiu model
X1       -4.0     0.0
"""

mol = Molecule(geometry=geom, nsp=1, nstates=2, dof=1, unit_pos='au', model=True)

qm = qm.model.Shin_Metiu(molecule=mol)

md = mqc.SHXF(molecule=mol, nsteps=2890, nesteps=1, dt=0.5, unit_dt='au', wsigma=0.1, istate=1, propagation="density")

md.run(molecule=mol, qm=qm, input_dir=f"./", save_scr=True, save_QMlog=False)
