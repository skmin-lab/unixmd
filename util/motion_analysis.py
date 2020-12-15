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
    parser.add_argument('-l', '-length', nargs="+", action='store', dest='length', type=int, \
        help="length between two atoms")
    parser.add_argument('-a', '-angle', nargs="+", action='store', dest='angles', type=int, \
        help="angle between three atoms")
    parser.add_argument('-d', '-dihedral', nargs="+", action='store', dest='dihedral', type=int, \
        help="dihedral angle between three atoms")
    parser.add_argument('-avg', action='store_true', help="additional option for averaging motion")

    args = parser.parse_args()

    # Indexing for numbering filename
    index = len(str(args.ntrajs))
    # Include step 0 
    nsteps1 = args.nsteps + 1

    # length analysis
    try:
        if (len(args.length) == 2):
            averaged_length(args.ntrajs, index, nsteps1, args.length, args.avg)
        else:
            raise ValueError (f"Invalid number of points! {len(args.length)}")
    except TypeError:
        pass

    # angle analysis
    try:
        if (len(args.angles) == 3):
            averaged_angle(args.ntrajs, index, nsteps1, args.anglei, args.avg)
        else:
            raise ValueError (f"Invalid number of points! {len(args.length)}")
    except TypeError:
        pass

    # dih angle analysis
    try:
        if (len(args.dihedral) == 4):
            averaged_dihedral(args.ntrajs, index, nsteps1, args.dihedral, args.avg)
        else:
            raise ValueError (f"Invalid number of points! {len(args.length)}")
    except TypeError:
        pass

def averaged_length(ntrajs, index, nsteps, points, avg):
    """ Averaging length between two points
    """
    
    if avg:
        f_write_avg = ""
        header_avg = f"#    Averaged length between atom {points[0]} and {points[1]}"
        f_write_avg += header_avg
        # define empty array for summation
        avg_length = np.zeros(nsteps)

    # define variable for count trajectories except halted trajectories
    mtrajs = ntrajs

    for itraj in range(ntrajs):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "MOVIE.xyz")


        with open(path, 'r') as f:
            line = f.read()
            lines = line.split()
        
        if (itraj==0):
            natoms = int(lines[0])
       
        try:
            f_write = ""
            header = f"#    length between atom {points[0]} and {points[1]}"
            f_write += header

            # sum over lengths between two points
            # TODO: change x,y,z(1,2) into proper number, depend on istep
            length = np.array([np.linalg.norm(lines[x1, y1, z1] - lines[x2, y2, z2]) for istep in range(nsteps)], dtype=np.float)

            if avg:
                avg_length += length 

            data = "".join(["\n" + f"{istep:8d}" + "".join(f"{length[istep]:15.8f}") for istep in range(nsteps)])
            f_write += data

            path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "LENGTH")
            typewriter(f_write, path)
            
        except ValueError:
            # exclude halted trajectories from total trajecotry number
            mtrajs -= 1
    
    if avg:
        # averaging array and print
        avg_length /= mtrajs
        avg_data = "".join(["\n" + f"{istep:8d}" + "".join(f"{avg_length[istep]:15.8f}") for istep in range(nsteps)])
        f_write_avg += avg_data
    
        typewriter(f_write_avg, "AVG_LENGTH")

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

