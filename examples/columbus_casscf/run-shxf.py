from molecule import Molecule
import bo
import mqc
from thermostat import *
from misc import data

geom = """
6
c2h4
C       -0.69284969      0.00000000      0.00000000   0.0 0.0 0.0
C        0.66284969      0.00000000      0.00000000   0.0 0.0 0.0
H       -1.03381870      0.72280446     -0.32280653   0.0 0.0 0.0
H       -1.23381870     -0.72280446     -0.34280653   0.0 0.0 0.0
H        1.03381870      0.75280446     -0.32280653   0.0 0.0 0.0
H        1.03381870     -0.72280446     -0.31280653   0.0 0.0 0.0
"""

mol = Molecule(geometry=geom, nstates=2)

qm = bo.columbus.CASSCF(molecule=mol, basis_set="6-31g*", memory="500", active_elec=2, active_orb=2, qm_path="/opt/Columbus7.0/Columbus", nthreads=1, version=7.0)

md = mqc.SHXF(molecule=mol, nsteps=1000, dt=0.125, istate=1, propagation="coefficient", wsigma=0.1, threshold=0.01)

bathT = rescale1(temperature=300.0, nrescale=20)

md.run(molecule=mol, theory=qm, thermostat=bathT, input_dir="./TRAJ.shxf", save_scr=True, save_QMlog=False)

