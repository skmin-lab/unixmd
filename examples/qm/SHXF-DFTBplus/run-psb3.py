from molecule import Molecule
import qm, mqc
from thermostat import *
from misc import data

geom = """
14
dt = 0.24 fs, nsteps = 1250, istate = 1, nesteps = 40000, T = None
C   -5.38182941     -2.47689908      0.12466734     -0.00040736      0.00009444      0.00005174
C   -2.71860203     -2.31117238     -0.02464885      0.00013840      0.00005165     -0.00029404
H   -6.54217849     -0.97759933     -0.03750662      0.00055738     -0.00320340      0.00003499
H   -5.95771390     -4.42845123      0.39948498     -0.00010513      0.00007888      0.00080635
H   -1.67116845     -4.22645713     -0.17152631      0.00176961     -0.00048002      0.00012280
C   -1.40258790      0.02442516     -0.34294572     -0.00009273     -0.00002609      0.00022486
C    1.16881813     -0.13095589     -0.56392899     -0.00053288      0.00021602     -0.00003640
H   -2.27108326      2.06447127     -0.40248367      0.00151569     -0.00009527     -0.00056093
H    2.10433200     -2.00723433     -0.58283237      0.00134830     -0.00116368      0.00052963
C    2.52950949      2.17723791     -0.13259035      0.00074505      0.00000632     -0.00001098
H    1.53565282      3.78669175     -0.92583694     -0.00081704      0.00160399     -0.00063959
N    5.02013328      2.35360852      0.77833964     -0.00028420     -0.00012433      0.00000286
H    5.76980481      4.34807789      0.85585859      0.00118364      0.00276347      0.00054755
H    6.39633409      1.09383098      1.23614811      0.00027646     -0.00185278     -0.00010874
"""

mol = Molecule(geometry=geom, nstates=2, charge=+1., unit_pos="au")

tuning = [3.2, 3.2, 3.2]

qm = qm.dftbplus.SSR(molecule=mol, scc_tol=1E-6, scc_max_iter=10000, \
    ocdftb=False, lcdftb=True, lc_method="MatrixBased", \
    ssr22=True, state_interactions=True, tuning=tuning, save_memory=True, \
    sk_path="/home/islee/program/dftb/slko/ob2-1-1/base/", \
    install_path="/opt/dftbplus-20.1/install-openmp/", version=20.1)

sigma = 0.1

md = mqc.SHXF(molecule=mol, nsteps=1250, dt=0.24, nesteps=40000, istate=1, propagation="density", threshold=0.02, \
    wsigma=sigma, vel_rescale="nac")

md.run(molecule=mol, qm=qm, input_dir="./TRAJ.psb3", debug=0)

