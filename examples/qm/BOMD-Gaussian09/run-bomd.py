from molecule import Molecule
import qm, mqc
from thermostat import *
from misc import data

geom = """
2
He2
He   -0.8   0.0   0.0   0.0 0.0 0.0
He    0.8   0.0   0.0   0.0 0.0 0.0
"""

mol = Molecule(geometry=geom, nstates=2)

qm = qm.gaussian09.DFT(molecule=mol, nthreads=4, memory="4gb", functional="B3LYP", basis_set="aug-cc-pVTZ", \
    g09_root_path="/opt/gaussian/")

md = mqc.BOMD(molecule=mol, nsteps=500, dt=0.125, istate=1)

bathT = rescale1(temperature=300.0, nrescale=10)

md.run(molecule=mol, qm=qm, thermostat=bathT, input_dir="./TRAJ.bomd", save_scr=True, save_QMlog=False)

