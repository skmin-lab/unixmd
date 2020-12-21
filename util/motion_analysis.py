import argparse
import os
import numpy as np

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
        help="angle between three atoms. second atom will be the centor atom")
    parser.add_argument('-d', '--dihedral', nargs="+", action='store', dest='dihedral', type=int, \
        help="dihedral angle between four or six atoms. Angle between (1,2,3),(2,3,4) or (1,2,3)(4,5,6) plane will be calculated")
    parser.add_argument('-m', '--mean', action='store_true', dest='l_mean', \
        help="additional option for averaging motion")
    args = parser.parse_args()

    # Indexing for numbering filename
    index = len(str(args.ntrajs))
    # Include step 0 
    nsteps1 = args.nsteps + 1
    # Checking job running
    if ((args.bond == None) and (args.angles == None) and (args.dihedral == None)):
        raise ValueError ("No analysis done -- check input arguments")

    # bond length analysis
    try:
        bond(args.ntrajs, index, nsteps1, args.bond, args.l_mean)
    except TypeError:
        pass

    # angle analysis
    try:
        angle(args.ntrajs, index, nsteps1, args.angles, args.l_mean)
    except TypeError:
        pass

    # dih angle analysis
    try:
        if (len(args.dihedral) == 4 or len(args.dihedral) == 6):
            dihedral(args.ntrajs, index, nsteps1, args.dihedral, args.l_mean)
        else:
            raise ValueError (f"Invalid number of points! {len(args.dihedral)}")
    except TypeError:
        pass

    
def bond(ntrajs, index, nsteps, points, l_mean):
    """ Averaging bond length between two points
    """
    if (l_mean == True):
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
        iline = 0

        f_write = ""
        # header file for individual trajectory analysis
        header = f"#    bond length between atom {points[0]} and {points[1]}"
        f_write += header

        bond = np.array([])
        with open(path, 'r') as f:
            lines = f.readline()
        natoms = int(lines)

        with open(path, 'r') as f:
             while (True):
                 lines = f.readline()
                 if not lines:
                     break
                 # save xyz coordinates in every isteps
                 if (iline % (natoms + 2) == points[0] + 2):
                     point1 = lines.split()
                 if (iline % (natoms + 2) == points[1] + 2):
                     point2 = lines.split()
                 # every istep, calculate it
                 if (iline % (natoms + 2) == natoms + 1):
                     # calculate bond length after both point1/2 extracted
                     bondlength = np.linalg.norm(np.array([point1[1], point1[2], point1[3]], dtype=float) \
                         - np.array([point2[1], point2[2], point2[3]], dtype=float))
                     bond = np.append(bond, [bondlength])
                 iline += 1
       
        if (iline != (nsteps * (2 + natoms))):
            mtrajs -= 1

        if (l_mean == True):
            # sum over bond lengths between two points if trajectory has full steps
            if (len(bond) == nsteps):
                mean_bond += bond

        data = "".join(["\n" + f"{istep:8d}" + "".join(f"{bond[istep]:15.8f}") for istep in range(len(bond))])
        f_write += data

        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "BOND")
        typewriter(f_write, path)
    
    if (l_mean == True):
        # averaging array and print
        mean_bond /= mtrajs
        mean_data = "".join(["\n" + f"{istep:8d}" + "".join(f"{mean_bond[istep]:15.8f}") for istep in range(nsteps)])
        f_write_mean += mean_data
        typewriter(f_write_mean, "AVG_BOND")


def angle(ntrajs, index, nsteps, points, l_mean):
    """ Averaging angle between two points
    """
    if (l_mean == True):
        f_write_mean = ""
        # header file for averaged trajectory analysis
        header_mean = f"#    Averaged angle between atom {points[0]}, {points[1]}, and {points[2]}"
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
        iline = 0

        f_write = ""
        # header file for individual trajectory analysis
        header = f"#    angle between atom {points[0]}, {points[1]}, and {points[2]}"
        f_write += header

        angle_array = np.array([])
        with open(path, 'r') as f:
            lines = f.readline()
        natoms = int(lines)

        with open(path, 'r') as f:
             while (True):
                 lines = f.readline()
                 if not lines:
                     break
                 # save xyz coordinate in every steps
                 if (iline % (2 + natoms) == points[0] + 2):
                     tmp = lines.split()
                     point1 = np.array([tmp[1], tmp[2], tmp[3]], dtype=float)
                 if (iline % (2 + natoms) == points[1] + 2):
                     tmp = lines.split()
                     point2 = np.array([tmp[1], tmp[2], tmp[3]], dtype=float)
                 if (iline % (2 + natoms) == points[2] + 2):
                     tmp = lines.split()
                     point3 = np.array([tmp[1], tmp[2], tmp[3]], dtype=float)
                 # every istep, calculate it
                 if (iline % (natoms + 2) == natoms + 1):
                     # calculate angle after point1/2/3 extracted, using two unit vector
                     unit_vector1 = (point1 - point2) / np.linalg.norm(point1 - point2)
                     unit_vector2 = (point3 - point2) / np.linalg.norm(point3 - point2)
                     dot_product = np.dot(unit_vector1, unit_vector2)
                     angle = np.degrees(np.arccos(dot_product))
                     angle_array = np.append(angle_array, [angle])
                 iline += 1
       
        if (iline != (nsteps * (2 + natoms))):
            mtrajs -= 1

        if (l_mean == True):
            # sum over angle between two points if trajectory has full steps
            if (len(angle_array) == nsteps):
                mean_angle += angle_array

        data = "".join(["\n" + f"{istep:8d}" + "".join(f"{angle_array[istep]:15.8f}") for istep in range(len(angle_array))])
        f_write += data

        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "ANGLE")
        typewriter(f_write, path)
    
    if (l_mean == True):
        # averaging array and print
        mean_angle /= mtrajs
        mean_data = "".join(["\n" + f"{istep:8d}" + "".join(f"{mean_angle[istep]:15.8f}") for istep in range(nsteps)])
        f_write_mean += mean_data
        typewriter(f_write_mean, "AVG_ANGLE")


def dihedral(ntrajs, index, nsteps, points, l_mean):
    """ Averaging dihedral angle between two points

 lines for calculate dihedral angle
                 # save xyz coordinates in every isteps
                 if (iline % (natoms + 2) == points[0] + 2):
                     point1 = lines.split()
                 if (iline % (natoms + 2) == points[1] + 2):
                     point2 = lines.split()
                 if (iline % (natoms + 2) == points[2] + 2):
                     point3 = lines.split()
                 if (iline % (natoms + 2) == points[3] + 2):
                     point4 = lines.split()
                 if (iline % (natoms + 2) == points[4] + 2):
                     point5 = lines.split()
                 if (iline % (natoms + 2) == points[5] + 2):
                     point6 = lines.split()
                 # every istep, calculate it
                 if (iline % (natoms + 2) == natoms + 1):
                     vector1_1 = point1 - point2
                     vector1_2 = point3 - point2

                     vector2_1 = point4 - point5
                     vector2_2 = point6 - point5

                     #find equation of plane which contains vector 1/2
                     a1, b1, c1 = np.cross(vector1_1, vector1_2)
                     a2, b2, c2 = np.cross(vector2_1, vector2_2)

                     dih_angle = np.degree(np.arccos(np.abs(a1 * a2 + b1 * b2 + c1 * c2) /  \ 
                         (np.sqrt(a1 ** 2 + b1 ** 2 + c1 ** 2) * np.sqrt(a2 ** 2 + b2 ** 2 + c2 ** 2))))
                     dih_array = np.append(dih_array, [dih_angle])
"""
    return 0

def typewriter(string, file_name):
    """ Function to write a string in filename
    """
    with open(file_name, "w") as f:
        f.write(string)

if (__name__ == "__main__"):
    motion_analysis()
