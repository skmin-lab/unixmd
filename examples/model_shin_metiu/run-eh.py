from molecule import Molecule
import bo
import mqc
from thermostat import *
from misc import data

data["X"] = 1836.0

geom = """
1
model
X       -4.0  0.0
"""

mol = Molecule(geometry=geom, nsp=1, nstates=2, dof=1, unit_pos='au', unit_vel='au', model=True)

qm = bo.model.Shin_Metiu(molecule=mol, qm_path="/home/islee/bin/shin-metiu/")

md = mqc.Eh(molecule=mol, nsteps=2000, dt=0.024, istate=1, propagation="density")

bathT = none()

md.run(molecule=mol, theory=qm, thermostat=bathT, input_dir="./TRAJ.eh", save_scr=True, save_QMlog=False)

