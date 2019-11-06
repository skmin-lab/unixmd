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

qm = bo.terachem.SSR(molecule=mol, ngpus=2, gpu_id="0 1", \
    functional="wPBEh", basis_set="6-31gs", scf_tol=1E-1, \
    reks22="yes", reks_scf_tol=5E-6, reks_max_scf_iter=10000, shift=1.0, \
    use_ssr_state=1, cpreks_max_tol=1E-6, cpreks_max_iter=10000, \
    qm_path="/home/filatov/terachem/bin/", version=1.92)
#    qm_path="/opt/TeraChem/bin/", version=1.93)

md = mqc.BOMD(molecule=mol, nsteps=1000, dt=0.125, istate=0)

bathT = rescale1(temperature=300.0, nrescale=20)

md.run(molecule=mol, theory=qm, thermostat=bathT, input_dir="./TRAJ.bomd", save_scr=True, save_QMlog=False)

