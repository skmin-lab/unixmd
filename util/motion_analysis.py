import argparse
import os
import numpy as np

def motion_analysis():
    """ Python utility script for UNI-xMD output analysis
        In this script, UNI-xMD MOVIE.xyz output files are post-process into given geometry criterion
    """
    parser = argparse.ArgumentParser(description="Python script for UNI-xMD motion analysis")
    parser.add_argument('-n', '-ntrajs', action='store', dest='ntrajs', type=int, \
        help="Total number of trajectories", required=True)
    parser.add_argument('-s', '-nsteps', action='store', dest='nsteps', type=int, \
        help="Total number of steps", required=True)
    parser.add_argument('-p', '-points', nargs="+", action='store', dest='points', type=int, \
        help="atom numbers for motion analysis", required=True)
    args = parser.parse_args()

    # Indexing for numbering filename
    index = len(str(args.ntrajs))
    # Indexing for job selection depend on points
    njob = len(args.points)
    # Include step 0 
    nsteps1 = args.nsteps + 1

    if (njob == 2):
        averaged_length(args.ntrajs, index, nsteps1, args.points)
    elif (njob == 3):
        averaged_angle(args.ntrajs, index, nsteps1, args.points)
    elif (njob == 4):
        averaged_dihedreal_angle(args.ntrajs, index, nsteps1, args.points)
    else:
        raise ValueError (f"Invalid number of points! {njob}")

def averaged_length(ntrajs, index, nsteps, points):
    """ Averaging length between two points
    """
    f_write = ""

    header = f"#    Averaged length between atom {points[0]} and {points[1]}"
    f_write += header

    # define empty array for summation
    avg_length = np.zeros(nsteps)
    # define variable for count trajectories except halted trajectories
    mtrajs = ntrajs

    for itraj in range(ntrajs):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "MOVIE.xyz")

        with open(path, 'r') as f:
            line = f.read()
            lines = line.split()
        try:
            # sum over lengths between two points
            # TODO: change x,y,z(1,2) into proper number, depend on istep
            avg_length += np.array([length(lines[x1, y1, z1], lines[x2, y2, z2]) for istep in range(nsteps)], dtype=np.float)
            
        except ValueError:
            # exclude halted trajectories from total trajecotry number
            mtrajs -= 1
    
    # averaging array and print
    avg_length /= mtrajs
    avg_data = "".join(["\n" + f"{istep:8d}" + "".join(f"{avg_length[istep]:15.8f}") for istep in range(nsteps)])
    f_write += avg_data

    typewriter(f_write, "AVG_LENGTH")

def length(atom1, atom2):
    """ calculating length between two points
    """
    return np.sqrt((atom1[0] - atom2[0]) ** 2 + (atom1[1] - atom2[1]) ** 2 + (atom1[2] - atom2[2]) ** 2)

def averaged_angle(ntrajs, index, nsteps, points):
    """ Averaging angle between two points
    """
    return 0

def averaged_dihedral_angle(ntrajs, index, nsteps, points):
    """ Averaging dihedral angle between two points
    """
    return 0

def typewriter(string, file_name):
    """ Function to write a string in filename
    """
    with open(file_name, "w") as f:
        f.write(string)

if (__name__ == "__main__"):
    motion_analysis()

