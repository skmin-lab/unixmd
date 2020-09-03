from molecule import Molecule
import qm, mm, mqc
from thermostat import *
from misc import data

geom = """
9
cnh4 with water
C         6.79554999     7.55605723     7.51409468     0.00013271     0.00013888    -0.00009336
N         8.05410625     7.30887334     7.51001345    -0.00018449    -0.00043517    -0.00005298
H         6.39647963     7.93655869     6.56275280    -0.00105895    -0.00041136     0.00057365
H         6.11259940     7.40974363     8.36320140     0.00028111    -0.00043589     0.00071349
H         8.64921781     7.03764388     8.31439330    -0.00051486     0.00088054     0.00049195
H         8.58738296     7.66765418     6.67063625     0.00091084     0.00044118     0.00088032
O         9.05007579    10.16245887     7.03576174     0.00012189     0.00024304    -0.00005104
H         8.95211500    11.11557303     7.04483075    -0.00013367     0.00023902     0.00029662
H         8.83723871     9.88618769     7.92834906    -0.00043007    -0.00017524    -0.00031385
"""

mol = Molecule(geometry=geom, nstates=2, charge=+1., qmmm=True, natoms_mm=3)

qm = qm.dftbplus.SSR(molecule=mol, scc_max_iter=10000, ocdftb=True, \
    periodic=False, do_charge="mechanical", nthreads=1, \
    sk_path="/home/islee/program/dftb/slko/mio-1-1/", \
    qm_path="/home/islee/program/dftb/dftbreks-19.1/", \
    script_path="/opt/dftbplus-19.1/tools/dptools/")

mm = mm.Tinker(molecule=mol, scheme="subtractive", periodic=False, \
    xyz_file="/home/islee/01-software-exercise/13-pyunixmd-test/57-cnh4-dftb_ssr-tinker-no_period/tinker.xyz", \
    key_file="/home/islee/01-software-exercise/13-pyunixmd-test/57-cnh4-dftb_ssr-tinker-no_period/tinker.key", \
    do_charge="mechanical", do_vdw="lennardjones", mm_path="/home/islee/program/tinker/bin/", nthreads=1, version=8.7)

md = mqc.SH(molecule=mol, nsteps=1000, dt=0.125, istate=1, propagation="density")

bathT = rescale1(temperature=300.0, nrescale=20)

md.run(molecule=mol, qm=qm, mm=mm, thermostat=bathT, input_dir="./TRAJ-sh", save_scr=False, save_QMlog=False, save_MMlog=False, debug=0)

