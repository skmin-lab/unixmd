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

qm = bo.dftbplus.DFTB(molecule=mol, max_scc_iter=300, sdftb=False, \
    sk_path="/home/islee/program/dftb/slko/3ob-3-1/", \
    qm_path="/home/islee/program/dftb/dftbplus-19.1/_install-openmp/bin/", \
    mpi=False, nthreads=1, ex_symmetry="S")
#    qm_path="/home/islee/program/dftb/dftbplus-19.1/_install-mpi/bin/", \
#    mpi=True, mpi_path="/opt/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/", nthreads=1)

md = mqc.BOMD(molecule=mol, nsteps=1000, dt=0.125, istate=0)

bathT = rescale1(temperature=300.0, nrescale=20)

md.run(molecule=mol, theory=qm, thermostat=bathT, input_dir="./TRAJ.bomd", save_scr=True, save_QMlog=False)

