from molecule import Molecule
import qm, mqc
from thermostat import NHC
import numpy as np
from qm.dftbplus.dftbpar import max_l

# This file is to perform molecular dynamics for periodic system with thermostat
# geom is geometry in a unit cell
geom = """
8

Si   -1.4988750000E-02  -1.5994220000E-02   2.8513260000E-02    -1.21470638     -0.90554557      2.57492454
Si    1.4159431200E+00   1.3921311500E+00   1.3188081300E+00     5.11077208      2.63739694     -3.28624412
Si    2.7274909000E+00   2.7225041400E+00   1.5346180000E-02     1.46824487      1.09363330      0.86635305
Si    4.0393830000E+00   4.0830272800E+00   1.3328931500E+00    -2.57181663      0.96285318     -1.60727177
Si    2.7075187300E+00  -9.5258200000E-03   2.7144786300E+00    -0.82685799     -0.82627047      0.23508769
Si    4.0447711900E+00   1.3545079200E+00   4.0816367600E+00    -1.99463869     -0.31838645      0.89824961
Si    9.7349400000E-03   2.6862434200E+00   2.6782821600E+00     0.78449756     -2.25206862     -3.08942237
Si    1.3514228700E+00   4.0683821300E+00   4.1113177300E+00    -0.75549481     -0.39161231      3.40832338
"""

mol = Molecule(geometry=geom, nstates=1, charge=0., unit_pos="A", unit_vel="A/ps")

# Define lattice vectors
l_tmp = """
   0.54270920000E+01   0.00000000000E+00   0.00000000000E+00
   0.00000000000E+00   0.54270920000E+01   0.00000000000E+00
   0.00000000000E+00   0.00000000000E+00   0.54270920000E+01
   """

# The data type of lattice must be list
lattice = list(np.array(l_tmp.split(), dtype=np.float))

# Define K-point
k_sampling = [2, 2, 2]

# p is used as the highest angular momentum for each Si atoms.
max_l["Si"] = "p"

# Define method of electronic structure calculation
qm = qm.dftbplus.DFTB(molecule=mol, scc_tol=1E-6, scc_max_iter=200, \
    ocdftb=False, lcdftb=False, lc_method="MatrixBased", \
    elec_temp=273.15, k_point=k_sampling, cell_length=lattice, periodic=True, \
    sk_path="/home/tikim/21-UNIXMD-work/example/05-Si-thermostat/", \
    install_path="/opt/dftbplus-20.1/install-mpi/", nthreads=4, version=20.1)

# Define MQC method
md = mqc.BOMD(molecule=mol, nsteps=10000, dt=1.0, istate=0, unit_dt="fs")

# Define thermostat method
bathT = NHC(temperature=273.15, coupling_strength=600.0)

md.run(molecule=mol, qm=qm, thermostat=bathT, input_dir="./TRAJ.NHC", save_QMlog=False, debug=0)

