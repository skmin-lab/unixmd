from molecule import Molecule
import qm, mqc
from thermostat import *
from misc import data

geom = """
6
dt = 0.125 fs, nsteps = 1000, istate = 2, nesteps = 10000, T = 300 K, nrescale = 20
C    2.86124018      4.31304115      4.15521402      0.00002671     -0.00010330      0.00004871
N    5.27829773      4.08368726      4.18677408     -0.00004087      0.00017187      0.00013434
H    1.88738261      3.13911452      5.54097988     -0.00100507      0.00019113     -0.00024745
H    1.76187927      5.49131646      2.81165546      0.00089679     -0.00006276     -0.00027862
H    6.54948319      4.69589988      2.84554252     -0.00034025     -0.00092241     -0.00072712
H    5.96455902      2.68641232      5.40075458      0.00068502     -0.00036011     -0.00118705
"""

mol = Molecule(geometry=geom, nstates=3, charge=+1, unit_pos="au")

qm = qm.columbus.CASSCF(molecule=mol, basis_set="6-31g*", memory="1000", \
    active_elec=6, active_orb=4, \
    qm_path="/opt/Columbus7.0/Columbus", nthreads=1, version=7.0)

md = mqc.SHXF(molecule=mol, nsteps=1000, nesteps=10000, dt=0.125, istate=2, propagation="density", wsigma=0.02)

bathT = rescale1(temperature=300.0, nrescale=20)

md.run(molecule=mol, qm=qm, thermostat=bathT, input_dir="./TRAJ.cnh4", save_scr=True, save_QMlog=False, save_MMlog=False, debug=0)

