import argparse
import math
import os
import numpy as np

def main():
    """ Python utility script for UNI-xMD output analysis
    """
    parser = argparse.ArgumentParser(description="Python script for UNI-xMD output analysis")
    parser.add_argument('-n', action='store', dest='ntraj', type=int, \
        help="Total trajectory number", required=True)
    args = parser.parse_args()

    index = len(str(args.ntraj))
    
    path = os.path.join(f"./TRAJ_{1:0{index}d}/md/", "MDENERGY")
    f1 = open(path, 'r') 
    lines_1 = f1.readlines()
    f1.close()

    nstep = len(lines_1) - 2

    State_avg(args.ntraj, index, nstep)
    Population_avg(args.ntraj, index, nstep)
    Coherence_avg(args.ntraj, index, nstep)
    NACME_avg(args.ntraj, index, nstep)
#    Motion_avg(ntraj, index, nstep, parameter_geom)


def State_avg(ntraj, index, nstep):
    wtmp = "#    Step    Average Running State\n"
    avg_state = np.zeros(nstep + 1)

    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "SHSTATE")
        f2 = open(path, 'r')
        lines_2 = f2.readlines()
        f2.close()

        for jstep in range(nstep + 1):
            tmp = lines_2[jstep + 1].split()
            avg_state[jstep] += float(tmp[1])

    for jstep in range(nstep + 1):
        tmp = avg_state[jstep] / ntraj
        wtmp += f"{jstep:8d}{tmp:15.8f}\n"
    typewriter(wtmp, "SHSTATE_avg")


def Population_avg(ntraj, index, nstep):
    wtmp = "#    Averaged Density Matrix\n"

    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "BOPOP")
        f2 = open(path, 'r')
        lines_2 = f2.readlines()
        f2.close()
        
        if (itraj == 0):
            tmp = lines_2[1].split()
            flength = len(tmp) - 1
            avg_state = np.zeros((nstep + 1, flength))

        for jstep in range(nstep + 1):
            tmp = lines_2[jstep + 1].split()

            for klength in range(flength):
                avg_state[jstep, klength] += float(tmp[klength + 1])

    for jstep in range(nstep + 1):
        w_str = f"{jstep:8d}"
        for klength in range(flength):
            tmp = avg_state[jstep, klength] / ntraj
            w_str += f"{tmp:15.8f}"
        wtmp += w_str + "\n"
    typewriter(wtmp, "BOPOP_avg")


def Coherence_avg(ntraj, index, nstep):
    wtmp = "#    Averaged Density Matrix: coherence Re-Im\n"

    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "BOCOH")
        f2 = open(path, 'r')
        lines_2 = f2.readlines()
        f2.close()
        
        if (itraj == 0):
            tmp = lines_2[1].split()
            flength = len(tmp) - 1
            avg_state = np.zeros((nstep + 1, flength))

        for jstep in range(nstep + 1):
            tmp = lines_2[jstep + 1].split()

            for klength in range(flength):
                avg_state[jstep, klength] += float(tmp[klength + 1])

    for jstep in range(nstep + 1):
        w_str = f"{jstep:8d}"
        for klength in range(flength):
            tmp = avg_state[jstep, klength] / ntraj
            w_str += f"{tmp:15.8f}"
        wtmp += w_str + "\n"
    typewriter(wtmp, "BOCOH_avg")
            
                                      
def NACME_avg(ntraj, index, nstep):
    wtmp = "#    Averaged Non-Adiabatic Coupling Matrix Eliments: off-diagonal\n"

    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "NACME")
        f2 = open(path, 'r')
        lines_2 = f2.readlines()
        f2.close()

        if (itraj == 0):
            tmp = lines_2[1].split()
            flength = len(tmp) - 1
            avg_state = np.zeros((nstep + 1, flength))

        for jstep in range(nstep + 1):
            tmp = lines_2[jstep + 1].split()

            for klength in range(flength):
                avg_state[jstep, klength] += float(math.sqrt(float(tmp[klength + 1]) ** 2))

    for jstep in range(nstep + 1):
        w_str = f"{jstep:8d}"
        for klength in range(flength):
            tmp = avg_state[jstep, klength] / ntraj
            w_str += f"{tmp:15.8f}"
        wtmp += w_str + "\n"
    typewriter(wtmp, "NACME_avg")


def typewriter(string, file_name):
    with open(file_name, "a") as f:
        f.write(string + "\n")


if __name__ == "__main__":
    main()

