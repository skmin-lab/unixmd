from molecule import Molecule
import bo
import mqc
from thermostat import *
from misc import data

geom = """
6
cnh4
 C  6.8088315554968402        7.5846355226161748        7.4985749560502279    -0.00009406545790       -0.00001175808276       -0.00008045786656
 N  8.0820844301057555        7.4049613410693684        7.4795523681544633     0.00009164789009        0.00017823604315        0.00004089607493
 H  6.2248402854142624        7.9685271260121411        6.6430441079383051    -0.00086290628068       -0.00001443847445        0.00014276449193
 H  6.2982583388968596        7.2002046125366634        8.3777424223510639     0.00119435291357       -0.00019677866340        0.00005599413463
 H  8.4920386966432595        7.1534710928315768        8.3940117651209079    -0.00077470161231       -0.00097156329491        0.00038782783683
 H  8.7629892216433749        7.7816013280546343        6.8069125442917393     0.00029744650826       -0.00115008307530       -0.00021296181301
"""

mol = Molecule(geometry=geom, nstates=3, charge=+1.)

qm = bo.dftbplus.DFTB(molecule=mol, scc_tol=1E-6, scc_max_iter=200, \
    guess="h0", sk_path="/home/islee/program/dftb/slko/mio-1-1/", \
    qm_path="/home/islee/program/dftb/dftbreks-19.1/", nthreads=1, \
    script_path="/opt/dftbplus-19.1/tools/dptools/", version=19.1)

md = mqc.SHXF(molecule=mol, nsteps=400, dt=0.125, istate=2, propagation="density")

bathT = rescale1(temperature=300.0, nrescale=20)

md.run(molecule=mol, theory=qm, thermostat=bathT, input_dir="./TRAJ-shxf", save_scr=True, save_QMlog=False, debug=0)

