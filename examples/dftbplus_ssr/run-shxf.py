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

qm = bo.dftbplus.SSR(molecule=mol, scc_tol=1E-6, max_scc_iter=10000, sdftb=True, lcdftb=True, lc_method="MM", ocdftb=True, \
    sk_path="/home/islee/program/dftb/slko/ob2-1-1/base/", \
    qm_path="/home/islee/program/dftb/dftbreks-19.1/", nthreads=1, \
    shift=1.0, grad_level=1, use_ssr_state=1, periodic=False, a_axis=10., b_axis=15., c_axis=20.)

md = mqc.SHXF(molecule=mol, nsteps=1000, dt=0.125, istate=1, propagation="coefficient", \
    wsigma=0.1, threshold=0.01)

bathT = rescale1(temperature=300.0, nrescale=20)

md.run(molecule=mol, theory=qm, thermostat=bathT, input_dir="./TRAJ.shxf", save_scr=True, save_QMlog=False)

