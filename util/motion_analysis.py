import argparse
import os
import numpy as np

import tracemalloc

def motion_analysis():
    """ Python utility script for UNI-xMD output analysis
        In this script, UNI-xMD MOVIE.xyz output files are post-process into given geometry criterion
    """
    parser = argparse.ArgumentParser(description="Python script for UNI-xMD motion analysis")
    parser.add_argument('-n', '--ntrajs', action='store', dest='ntrajs', type=int, \
        help="Total number of trajectories", required=True)
    parser.add_argument('-s', '--nsteps', action='store', dest='nsteps', type=int, \
        help="Total number of steps", required=True)
    parser.add_argument('-b', '--bond', nargs=2, action='store', dest='bond', type=int, \
        help="bond length between two atoms")
    parser.add_argument('-a', '--angle', nargs=3, action='store', dest='angles', type=int, \
        help="angle between three atoms")
    parser.add_argument('-d', '--dihedral', nargs="+", action='store', dest='dihedral', type=int, \
        help="dihedral angle between three atoms")
    parser.add_argument('-m', '--mean', action='store_true', dest='mean', \
        help="additional option for averaging motion")

    args = parser.parse_args()

    # Indexing for numbering filename
    index = len(str(args.ntrajs))
    # Include step 0 
    nsteps1 = args.nsteps + 1
    # Checking job running
    chk = 3

    bond(args.ntrajs, index, nsteps1, args.bond, args.mean)
    # bond length analysis
    try:
        bond(args.ntrajs, index, nsteps1, args.bond, args.mean)
    except TypeError:
        chk -= 1
        pass

    # angle analysis
    try:
        angle(args.ntrajs, index, nsteps1, args.angles, args.mean)
    except TypeError:
        chk -= 1
        pass

    # dih angle analysis
    try:
        if (len(args.dihedral) == 4 or len(args.dihedral) == 6):
            dihedral(args.ntrajs, index, nsteps1, args.dihedral, args.mean)
        else:
            raise ValueError (f"Invalid number of points! {len(args.dihedral)}")
    except TypeError:
        chk -= 1
        pass

    if (chk == 0):
        raise ValueError ("No analysis done -- check input arguments")
    
def bond(ntrajs, index, nsteps, points, mean):
    """ Averaging bond length between two points
    """
    if mean:
        f_write_mean = ""
        # header file for averaged trajectory analysis
        header_mean = f"#    Averaged bond length between atom {points[0]} and {points[1]}"
        f_write_mean += header_mean
        # define empty array for summation
        mean_bond = np.zeros(nsteps)

    # natom index fix: input variable #(x) atom is reading as #(x+1) in python
    points = np.array(points)
    points -= 1

    # define variable for count trajectories except halted trajectories
    mtrajs = ntrajs

    for itraj in range(ntrajs):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "MOVIE.xyz")
        # chacking line number
        nline = 0

        try:
            f_write = ""
            # header file for individual trajectory analysis
            header = f"#    bond length between atom {points[0]} and {points[1]}"
            f_write += header

            bond = []
            with open(path, 'r') as f:
                lines = f.readline()
            natoms = int(lines)

            with open(path, 'r') as f:
                while (nline < (nsteps * (2 + natoms))):
                    lines = f.readline()
                    # save xyz coordinate in every steps
                    if (nline % (2 + natoms) == points[0] + 2):
                        tmp = lines.split()
                        point1 = np.array([tmp[1], tmp[2], tmp[3]], dtype=float)
                    if (nline % (2 + natoms) == points[1] + 2):
                        tmp = lines.split()
                        point2 = np.array([tmp[1], tmp[2], tmp[3]], dtype=float)

                        # calculate bond length after both point1/2 extracted
                        bondlength = np.linalg.norm(point1 - point2)
                        bond.append(bondlength)
                    nline += 1
                    bond = np.array(bond)           
            if mean:
                # sum over bond lengths between two points
                mean_bond += bond
            data = "".join(["\n" + f"{istep:8d}" + "".join(f"{bond[istep]:15.8f}") for istep in range(nsteps)])
            f_write += data

            path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "BOND")
            typewriter(f_write, path)
            
        except ValueError:
            # exclude halted trajectories from total trajecotry number
            mtrajs -= 1
    
    if mean:
        # averaging array and print
        mean_bond /= mtrajs
        mean_data = "".join(["\n" + f"{istep:8d}" + "".join(f"{mean_bond[istep]:15.8f}") for istep in range(nsteps)])
        f_write_mean += mean_data
    
        typewriter(f_write_mean, "AVG_BOND")

def angle(ntrajs, index, nsteps, points, means):
    """ Averaging angle between two points
    """
    if mean:
        f_write_mean = ""
        # header file for averaged trajectory analysis
        header_mean = f"#    Averaged angle length between atom {points[0]} and {points[1]}"
        f_write_mean += header_mean
        # define empty array for summation
        mean_angle = np.zeros(nsteps)

    # natom index fix: input variable #(x) atom is reading as #(x+1) in python
    points = np.array(points)
    points -= 1

    # define variable for count trajectories except halted trajectories
    mtrajs = ntrajs

    for itraj in range(ntrajs):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "MOVIE.xyz")
        # chacking line number
        nline = 0

        try:
            f_write = ""
            # header file for individual trajectory analysis
            header = f"#    angle between atom {points[0]}, {points[1]}, and {points[2]}"
            f_write += header

            angle_array = []
            with open(path, 'r') as f:
                lines = f.readline()
            natoms = int(lines)

            with open(path, 'r') as f:
                while (nline < (nsteps * (2 + natoms))):
                    lines = f.readline()
                    # save xyz coordinate in every steps
                    if (nline % (2 + natoms) == points[0] + 2):
                        tmp = lines.split()
                        point1 = np.array([tmp[1], tmp[2], tmp[3]], dtype=float)
                    if (nline % (2 + natoms) == points[1] + 2):
                        tmp = lines.split()
                        point2 = np.array([tmp[1], tmp[2], tmp[3]], dtype=float)
                    if (nline % (2 + natoms) == points[2] + 2):
                        tmp = lines.split()
                        point3 = np.array([tmp[1], tmp[2], tmp[3]], dtype=float)

                        # calculate angle after point1/2/3 extracted using two unit vector
                        unit_vector1 = (point1 - point2)/np.linalg.norm(point1 - point2)
                        unit_vector2 = (point3 - point2)/np.linalg.norm(point3 - point2)
                        angle = np.linalg.norm(point1 - point2)
                        angle_array.append(angle)
                    nline += 1
                    angle_array = np.array(angle_array)           
            if mean:
                # sum over angle between three points
                mean_angle += angle_array
            data = "".join(["\n" + f"{istep:8d}" + "".join(f"{angle_array[istep]:15.8f}") for istep in range(nsteps)])
            f_write += data

            path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "ANGLE")
            typewriter(f_write, path)
            
        except ValueError:
            # exclude halted trajectories from total trajecotry number
            mtrajs -= 1
    
    if mean:
        # averaging array and print
        mean_bond /= mtrajs
        mean_data = "".join(["\n" + f"{istep:8d}" + "".join(f"{mean_bond[istep]:15.8f}") for istep in range(nsteps)])
        f_write_mean += mean_data
    
        typewriter(f_write_mean, "AVG_BOND")
    return 0

def dihedral(ntrajs, index, nsteps, points, means):
    """ Averaging dihedral angle between two points
    """
    return 0

def typewriter(string, file_name):
    """ Function to write a string in filename
    """
    with open(file_name, "w") as f:
        f.write(string)

if (__name__ == "__main__"):

    tracemalloc.start()
    motion_analysis()
    
    size, peak = tracemalloc.get_traced_memory()
    print(f"MEMORY PEAK = {peak}")
