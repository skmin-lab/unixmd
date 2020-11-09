import argparse
import os
import numpy as np

def statistical_analysis():
    """ Python utility script for UNI-xMD output analysis
    """
    parser = argparse.ArgumentParser(description="Python script for UNI-xMD output analysis")
    parser.add_argument('-n', '-ntraj', action='store', dest='ntraj', type=int, \
        help="Total trajectory number", required=True)
    parser.add_argument('-s', '-nstep', action='store', dest='nstep', type=int, \
        help="Total step number for analysis", required=True)
    parser.add_argument('-t', '-nstate', action='store', dest='nstate', type=int, \
        help="Total state number for analysis", required=True)
    args = parser.parse_args()

    index = len(str(args.ntraj))
    
    # Include step 0 
    nstep1 = args.nstep + 1
    
    State_avg(args.ntraj, index, nstep1, args.nstate)
    Population_avg(args.ntraj, index, nstep1, args.nstate)
    Coherence_avg(args.ntraj, index, nstep1, args.nstate)
    NACME_avg(args.ntraj, index, nstep1, args.nstate)
#    Motion_avg(ntraj, index, nstep, parameter_geom)


def State_avg(ntraj, index, nstep, nstate):
    f_write = ""

    header = "#    Running state based averaged BO population"
    f_write += header
        
    avg_state = np.zeros((nstate, nstep))

    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "SHSTATE")

        with open(path, 'r') as f:
            # Skip header and read rest
            line = f.readline()
            line = f.read()
            lines = line.split()

        rstates = np.array(lines[1::2][:nstep], dtype=np.int)    
        avg_state += np.array([(rstates == ist) for ist in range(nstate)], dtype=np.float)

    avg_state /= ntraj
    avg_data = "".join([("\n" + f"{istep:8d}" + "".join([f"{avg_state[istate, istep]:15.8f}" \
        for istate in range(nstate)])) for istep in range(nstep)])
    f_write += avg_data

    typewriter(f_write, "BOSTATE_avg")


def Population_avg(ntraj, index, nstep, nstate):
    f_write = ""

    header = "#    Diagonal density matrix element based averaged BO population"
    f_write += header
    avg_pop = np.zeros((nstate, nstep))
    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "BOPOP")

        with open(path, 'r') as f:
            # Skip header and read rest
            line = f.readline()
            line = f.read()
            lines = line.split()

        avg_pop += np.array([lines[istate::(nstate + 1)][:nstep] for istate in range(1, nstate + 1)], dtype=np.float)

    avg_pop /= ntraj
    avg_data = "".join([("\n" + f"{istep:8d}" + "".join([f"{avg_pop[istate, istep]:15.8f}" \
        for istate in range(nstate)])) for istep in range(nstep)])
    f_write += avg_data

    typewriter(f_write, "BOPOP_avg")


def Coherence_avg(ntraj, index, nstep, nstate):
    f_write = ""

    header = "#    Averaged electronic coherence"
    f_write += header

    nstate_pair = int(nstate * (nstate - 1) / 2)
    avg_coh = np.zeros((nstate_pair, nstep))

    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "BOPOP")

        with open(path, 'r') as f:
            # Skip header and read rest
            line = f.readline()
            line = f.read()
            lines = line.split()
            lines = list(map(float, lines))

        avg_coh += np.array([np.multiply(lines[istate::(nstate + 1)][:nstep],lines[jstate::(nstate + 1)][:nstep]) \
            for istate in range(1, nstate + 1) for jstate in range(istate + 1, nstate + 1)])
 
    avg_coh /= ntraj
    avg_data = "".join([("\n" + f"{istep:8d}" + "".join([f"{avg_coh[istate, istep]:15.8f}" \
        for istate in range(nstate_pair)])) for istep in range(nstep)])
    f_write += avg_data

    typewriter(f_write, "BOCOH_avg")
            
                                      
def NACME_avg(ntraj, index, nstep, nstate):
    f_write = ""

    header = "#    Averaged Non-Adiabatic Coupling Matrix Eliments: off-diagonal"
    f_write += header

    nstate_pair = int(nstate * (nstate - 1) / 2)
    avg_nacme = np.zeros((nstate_pair, nstep))

    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "NACME")

        with open(path, 'r') as f:
            # Skip header and read rest
            line = f.readline()
            line = f.read()
            lines = line.split()

        avg_nacme += abs(np.array([lines[istate::(nstate_pair + 1)][:nstep] for istate in range(1, nstate_pair + 1)], dtype=np.float))

    avg_nacme /= ntraj
    avg_data = "".join([("\n" + f"{istep:8d}" + "".join([f"{avg_nacme[istate, istep]:15.8f}" \
        for istate in range(nstate_pair)])) for istep in range(nstep)])
    f_write += avg_data

    typewriter(f_write, "NACME_avg")


def typewriter(string, file_name):
    with open(file_name, "a") as f:
        f.write(string + "\n")


if (__name__ == "__main__"):
    statistical_analysis()

