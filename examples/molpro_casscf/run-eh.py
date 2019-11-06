from molecule import Molecule
import bo
import mqc
from thermostat import *
from misc import data

geom = """
6
c2h4
C       -0.69284974      0.00000000      0.00000000   0.0 0.0 0.0
C        0.66284974      0.00000000      0.00000000   0.0 0.0 0.0
H       -1.03381877      0.72280451     -0.32280655   0.0 0.0 0.0
H       -1.23381879     -0.72280451     -0.34280655   0.0 0.0 0.0
H        1.03381877      0.75280451     -0.32280655   0.0 0.0 0.0
H        1.03381877     -0.72280451     -0.31280655   0.0 0.0 0.0
"""

mol = Molecule(geometry=geom, nstates=2)

qm = bo.molpro.CASSCF(molecule=mol, basis_set="6-31G", memory="500m", \
     active_elec=2, active_orb=2, cpscf_grad_tol=1E-7, \
     qm_path="/opt/molpro2015.1/bin/", nthreads=1, version=2015.1)

md = mqc.Eh(molecule=mol, nsteps=1000, dt=0.125, istate=1, propagation="density")

bathT = rescale1(temperature=300.0, nrescale=20)

md.run(molecule=mol, theory=qm, thermostat=bathT, input_dir="./TRAJ.eh", save_scr=True, save_QMlog=False)

