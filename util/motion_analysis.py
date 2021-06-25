import argparse
import os, sys
import numpy as np

def motion_analysis():
    """ Utility script for PyUNIxMD output analysis
        In this script, PyUNIxMD MOVIE.xyz output files are post-process into given geometry criterion
    """
    parser = argparse.ArgumentParser(description="Utility script for PyUNIxMD motion analysis")
    parser.add_argument('-n', '--ntrajs', action='store', dest='ntrajs', type=int, \
        help="Total number of trajectories.", required=True)
    parser.add_argument('-s', '--nsteps', action='store', dest='nsteps', type=int, \
        help="Total number of steps.", required=True)
    parser.add_argument('-b', '--bond', nargs=2, action='store', dest='bond', type=int, \
        help="Bond length between two atoms.")
    parser.add_argument('-a', '--angle', nargs=3, action='store', dest='angle', type=int, \
        help="Angle between three atoms. second atom will be the centor atom.")
    parser.add_argument('-d', '--dihedral', nargs="+", action='store', dest='dihedral', type=int, \
        help="Dihedral angle between four or six atoms. In four atom case, angle between (1,2,3),(2,3,4) plane will be calculated and dihedreal axis will be atom2-atom3. \
        In six atom case, angle between (1,2,3),(4,5,6) plane will be calculated and dihedral axis will be atom3-atom4.")
    parser.add_argument('-m', '--mean', action='store_true', dest='l_mean', \
        help="Additional option for averaging motion.")
    args = parser.parse_args()

    # Indexing for numbering filename
    digit = len(str(args.ntrajs))
    # Include step 0
    nsteps1 = args.nsteps + 1
    # Checking job running
    if ((args.bond == None) and (args.angle == None) and (args.dihedral == None)):
        raise ValueError (f"( {sys._getframe().f_code.co_name} ) No analysis done -- check input arguments")

    # Bond length analysis
    if (args.bond != None):
        calculate_bond_length(args.ntrajs, digit, nsteps1, args.bond, args.l_mean)

    # Angle analysis
    if (args.angle != None):
        calculate_angle(args.ntrajs, digit, nsteps1, args.angle, args.l_mean)

    # Dihedral angle analysis
    if (args.dihedral != None):
        if (len(args.dihedral) == 4 or len(args.dihedral) == 6):
            calculate_dihedral(args.ntrajs, digit, nsteps1, args.dihedral, args.l_mean)
        else:
            raise ValueError (f"( {sys._getframe().f_code.co_name} ) Invalid length of 'args.dihedral'! {len(args.dihedral)}")

def calculate_bond_length(ntrajs, digit, nimages, atom_index, l_mean):
    """ Averaging bond length between two points
    """
    if (l_mean):
        f_write_mean = ""
        # header file for averaged trajectory analysis
        header_mean = f"#    Averaged bond length between atom {atom_index[0]} and {atom_index[1]}"
        f_write_mean += header_mean
        # define empty array for summation
        mean_bond = np.zeros(nimages)

    # define variable for count trajectories except halted trajectories
    mtrajs = ntrajs

    for itraj in range(ntrajs):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{digit}d}/md/", "MOVIE.xyz")
        # checking line number
        iline = 0

        f_write = ""
        # header file for individual trajectory analysis
        header = f"#    Bond length between atom {atom_index[0]} and {atom_index[1]}"
        f_write += header
        bond_list = []

        with open(path, 'r') as f:
            line = f.readline()
            natoms = int(line)
            count = 0
            while (True):
                line = f.readline()
                if ((not line) or (count >= nimages)):
                    break
                # save xyz coordinates in every isteps
                if (iline % (natoms + 2) == atom_index[0]):
                    atom1 = np.array(line.split()[1:4], dtype=np.float64)
                elif (iline % (natoms + 2) == atom_index[1]):
                    atom2 = np.array(line.split()[1:4], dtype=np.float64)
                # every istep, calculate it
                if (iline % (natoms + 2) == natoms):
                    # calculate bond length after both point1/2 extracted
                    bond = np.linalg.norm(atom1 - atom2)
                    bond_list += [bond]
                    count += 1
                iline += 1

        if (iline != (nimages * (2 + natoms) - 1)):
            mtrajs -= 1

        if (l_mean):
            # sum over bond lengths between two points if trajectory has full steps
            if (len(bond_list) == nimages):
                mean_bond += np.array(bond_list)

        # save data even if the trajectory is halted
        data = "".join(["\n" + f"{istep:8d}" + "".join(f"{bond_dat:15.8f}") for istep, bond_dat in enumerate(bond_list)])
        f_write += data

        path = os.path.join(f"./TRAJ_{itraj + 1:0{digit}d}/md/", "BOND")
        typewriter(f_write, path)

    if (l_mean):
        # averaging array and print
        mean_bond /= mtrajs
        mean_data = "".join(["\n" + f"{istep:8d}" + "".join(f"{mean_bond[istep]:15.8f}") for istep in range(nimages)])
        f_write_mean += mean_data
        typewriter(f_write_mean, "AVG_BOND")

def calculate_angle(ntrajs, digit, nimages, atom_index, l_mean):
    """ Averaging angle between two points
    """
    if (l_mean):
        f_write_mean = ""
        # header file for averaged trajectory analysis
        header_mean = f"#    Averaged angle between atom {atom_index[0]}, {atom_index[1]}, and {atom_index[2]}"
        f_write_mean += header_mean
        # define empty array for summation
        mean_angle = np.zeros(nimages)

    # define variable for count trajectories except halted trajectories
    mtrajs = ntrajs

    for itraj in range(ntrajs):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{digit}d}/md/", "MOVIE.xyz")
        # checking line number
        iline = 0

        f_write = ""
        # header file for individual trajectory analysis
        header = f"#    Angle between atom {atom_index[0]}, {atom_index[1]}, and {atom_index[2]}"
        f_write += header
        angle_list = []

        with open(path, 'r') as f:
            line = f.readline()
            natoms = int(line)
            count = 0
            while (True):
                line = f.readline()
                if ((not line) or (count >= nimages)):
                    break
                # save xyz coordinates in every isteps
                if (iline % (natoms + 2) == atom_index[0]):
                    atom1 = np.array(line.split()[1:4], dtype=np.float64)
                elif (iline % (natoms + 2) == atom_index[1]):
                    atom2 = np.array(line.split()[1:4], dtype=np.float64)
                elif (iline % (natoms + 2) == atom_index[2]):
                    atom3 = np.array(line.split()[1:4], dtype=np.float64)
                # every istep, calculate it
                if (iline % (natoms + 2) == natoms):
                    # calculate angle with vector calculation
                    unit_vector1 = (atom1 - atom2) / np.linalg.norm(atom1 - atom2)
                    unit_vector2 = (atom3 - atom2) / np.linalg.norm(atom3 - atom2)

                    x = np.dot(unit_vector1, unit_vector2)
                    y = np.linalg.norm(np.cross(unit_vector1, unit_vector2))
                    angle = np.degrees(np.arctan2(y, x))
                    angle_list += [angle]
                    count += 1
                iline += 1

        if (iline != (nimages * (2 + natoms) - 1)):
            mtrajs -= 1

        if (l_mean):
            # sum over angles between three points if trajectory has full steps
            if (len(angle_list) == nimages):
                mean_angle += np.array(angle_list)

        # save data even if the trajectory is halted
        data = "".join(["\n" + f"{istep:8d}" + "".join(f"{angle_dat:15.8f}") for istep, angle_dat in enumerate(angle_list)])
        f_write += data

        path = os.path.join(f"./TRAJ_{itraj + 1:0{digit}d}/md/", "ANGLE")
        typewriter(f_write, path)

    if (l_mean):
        # averaging array and print
        mean_angle /= mtrajs
        mean_data = "".join(["\n" + f"{istep:8d}" + "".join(f"{mean_angle[istep]:15.8f}") for istep in range(nimages)])
        f_write_mean += mean_data
        typewriter(f_write_mean, "AVG_ANGLE")

def calculate_dihedral(ntrajs, digit, nimages, atom_index, l_mean):
    """ Averaging dihedral angle between two points
    """
    if (l_mean):
        f_write_mean = ""
        # header file for averaged trajectory analysis
        if (len(atom_index) == 4):
            header_mean = f"#    Averaged dihedral angle between plane ({atom_index[0]}, {atom_index[1]}, {atom_index[2]}) atoms and ({atom_index[1]}, {atom_index[2]}, {atom_index[3]}) atoms"
        elif (len(atom_index) == 6):
            header_mean = f"#    Averaged dihedral angle between plane ({atom_index[0]}, {atom_index[1]}, {atom_index[2]}) atoms and ({atom_index[3]}, {atom_index[4]}, {atom_index[5]}) atoms"
        f_write_mean += header_mean
        # define empty array for summation
        mean_dihedral = np.zeros(nimages)

    # define variable for count trajectories except halted trajectories
    mtrajs = ntrajs

    for itraj in range(ntrajs):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{digit}d}/md/", "MOVIE.xyz")
        # checking line number
        iline = 0

        f_write = ""
        # header file for individual trajectory analysis
        if (len(atom_index) == 4):
            header = f"#    Dihedral angle between plane ({atom_index[0]}, {atom_index[1]}, {atom_index[2]}) atoms and ({atom_index[1]}, {atom_index[2]}, {atom_index[3]}) atoms"
        elif (len(atom_index) == 6):
            header = f"#    Dihedral angle between plane ({atom_index[0]}, {atom_index[1]}, {atom_index[2]}) atoms and ({atom_index[3]}, {atom_index[4]}, {atom_index[5]}) atoms"
        f_write += header
        dihedral_list = []

        with open(path, 'r') as f:
            line = f.readline()
            natoms = int(line)
            count = 0
            while (True):
                line = f.readline()
                if ((not line) or (count >= nimages)):
                    break
                # save xyz coordinates in every isteps
                if (iline % (natoms + 2) == atom_index[0]):
                    atom1 = np.array(line.split()[1:4], dtype=np.float64)
                elif (iline % (natoms + 2) == atom_index[1]):
                    atom2 = np.array(line.split()[1:4], dtype=np.float64)
                elif (iline % (natoms + 2) == atom_index[2]):
                    atom3 = np.array(line.split()[1:4], dtype=np.float64)
                elif (iline % (natoms + 2) == atom_index[3]):
                    atom4 = np.array(line.split()[1:4], dtype=np.float64)
                if (len(atom_index) == 6):
                    if (iline % (natoms + 2) == atom_index[4]):
                        atom5 = np.array(line.split()[1:4], dtype=np.float64)
                    elif (iline % (natoms + 2) == atom_index[5]):
                        atom6 = np.array(line.split()[1:4], dtype=np.float64)
                # every istep, calculate it
                if (iline % (natoms + 2) == natoms):
                    if (len(atom_index) == 4):
                        b1 = atom2 - atom1
                        b2 = atom3 - atom2
                        b3 = atom4 - atom3
                    elif (len(atom_index) == 6):
                        b1 = np.cross((atom3 - atom1), (atom3 - atom2))
                        b2 = atom4 - atom3
                        b3 = np.cross((atom4 - atom5), (atom4 - atom6))

                    norm1 = np.cross(b1, b2)
                    norm1 /= np.linalg.norm(norm1)
                    norm2 = np.cross(b2, b3)
                    norm2 /= np.linalg.norm(norm2)

                    m = np.cross(norm1, b2 / np.linalg.norm(b2))
                    x = np.dot(norm1, norm2)
                    y = np.dot(m, norm2)
                    dihedral_angle = np.degrees(np.arctan2(y, x))
                    dihedral_list += [dihedral_angle]
                    count += 1
                iline += 1

        if (iline != (nimages * (2 + natoms) - 1)):
            mtrajs -= 1

        if (l_mean):
            # sum over dihedral angles if trajectory has full steps
            if (len(dihedral_list) == nimages):
                mean_dihedral += np.array(dihedral_list)

        # save data even if the trajectory is halted
        data = "".join(["\n" + f"{istep:8d}" + "".join(f"{dihedral_dat:15.8f}") for istep, dihedral_dat in enumerate(dihedral_list)])
        f_write += data

        path = os.path.join(f"./TRAJ_{itraj + 1:0{digit}d}/md/", "DIHEDRAL")
        typewriter(f_write, path)

    if (l_mean):
        # averaging array and print
        mean_dihedral /= mtrajs
        mean_data = "".join(["\n" + f"{istep:8d}" + "".join(f"{mean_dihedral[istep]:15.8f}") for istep in range(nimages)])
        f_write_mean += mean_data
        typewriter(f_write_mean, "AVG_DIHEDRAL")

def typewriter(string, file_name):
    """ Function to write a string in filename
    """
    with open(file_name, "w") as f:
        f.write(string)

if (__name__ == "__main__"):
    motion_analysis()
