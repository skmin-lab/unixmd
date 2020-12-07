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

    # Indexing for job selection depend on points
    njob = len(args.points)) 
    # Include step 0 
    nsteps1 = args.nsteps + 1

    if (njob == 2):
        averaged_length(args.ntrajs,nsteps1, args.points)
    elif (njob == 3):
        averaged_angle(args.ntrajs, nsteps1, args.points)
    elif (njob == 4):
        averaged_dihedreal_angle(args.ntrajs, nsteps1, args.points)
    else:
        raise ValueError (f"Invalid number of points! {njob}")

def averaged_length(ntrajs, nsteps, points):
    """ Averaging length between two points
    """
    return 0

def averaged_angle(ntrajs, nsteps, points):
    """ Averaging angle between two points
    """
    return 0

def averaged_dihedral_angle(ntrajs, nsteps, points):
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

