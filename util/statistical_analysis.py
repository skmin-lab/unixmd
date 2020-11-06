import argparse
import os
import numpy as np

def statistical_analysis():
    """ Python utility script for UNI-xMD output analysis
    """
    parser = argparse.ArgumentParser(description="Python script for UNI-xMD output analysis")
    parser.add_argument('-n', action='store', dest='ntraj', type=int, \
        help="Total trajectory number", required=True)
    parser.add_argument('-s', action='store', dest='nstep', type=int, \
        help="Total step number for analysis", required=True)
    parser.add_argument('-t', action='store', dest='nstate', type=int, \
        help="Total state number for analysis", required=True)
    args = parser.parse_args()

    index = len(str(args.ntraj))
    
    State_avg(args.ntraj, index, args.nstep, args.nstate)
    Population_avg(args.ntraj, index, args.nstep, args.nstate)
    Coherence_avg(args.ntraj, index, args.nstep, args.nstate)
    NACME_avg(args.ntraj, index, args.nstep, args.nstate)
#    Motion_avg(ntraj, index, nstep, parameter_geom)


def State_avg(ntraj, index, nstep, nstate):
    f_write = ""

    header = "#    Step    Average Running State\n"
    f_write += header

    # include step 0 
    nstep += 1
        
    avg_state = np.zeros((nstep, nstate))

    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "SHSTATE")

        with open(path, 'r') as f:
            lines = f.readlines()
        
        for istep in range(nstep):
            rstate = lines[istep + 1].split()[1]
            avg_state[istep, int(rstate)] += 1

    avg_state /= ntraj
   
    for istep in range(nstep):
        avg_data = f"{istep:8d}"
        for istate in range(nstate):
            avg_data += f"{avg_state[istep, istate]:15.8f}"
        f_write += avg_data + "\n"
    typewriter(f_write.rstrip(), "BOSTATE_avg")


def Population_avg(ntraj, index, nstep, nstate):
    f_write = ""

    header = "#    Averaged Density Matrix\n"
    f_write += header
  
    # include step 0 
    nstep += 1

    avg_pop = np.zeros((nstate, nstep))

    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "BOPOP")

        with open(path, 'r') as f:
            # skip header and read rest
            line = f.readline()
            line = f.read()
            lines = line.split()

        avg_pop += np.array([lines[(istate + 1)::(nstate + 1)][:nstep] for istate in range(nstate)], dtype=np.float)

    avg_pop /= ntraj

    for istep in range(nstep):
        avg_data = f"{istep:8d}"
        for istate in range(nstate):
            avg_data += f"{avg_pop[istate, istep]:15.8f}"
        f_write += avg_data + "\n"
    typewriter(f_write.rstrip(), "BOPOP_avg")


def Coherence_avg(ntraj, index, nstep, nstate):
    f_write = ""

    header = "#    Averaged Density Matrix: coherence Re-Im\n"
    f_write += header

    # include step 0 
    nstep += 1

    nstate_pair = int(nstate * (nstate - 1) / 2)
    avg_coh = np.zeros((nstep, nstate_pair))

    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "BOPOP")

        with open(path, 'r') as f:
            lines = f.readlines()
       
        for istep in range(nstep):
            i_pop = lines[istep + 1].split()

            ntmp = 0
            for istate in range(0, nstate):
                for jstate in range(istate + 1, nstate):
                    # i_pop[0] is stepnumber
                    avg_coh[istep, ntmp] += float(i_pop[istate + 1]) * float(i_pop[jstate + 1])
                    ntmp += 1

    avg_coh /= ntraj

    for istep in range(nstep):
        avg_data = f"{istep:8d}"
        for istate in range(nstate):
            avg_data += f"{avg_coh[istep, istate]:15.8f}"
        f_write += avg_data + "\n"
    typewriter(f_write.rstrip(), "BOCOH_avg")
            
                                      
def NACME_avg(ntraj, index, nstep, nstate):
    f_write = ""

    header = "#    Averaged Non-Adiabatic Coupling Matrix Eliments: off-diagonal\n"
    f_write += header
    
    # include step 0
    nstep += 1
    
    nstate_pair = int(nstate * (nstate - 1) / 2)
    avg_nacme = np.zeros((nstate_pair, nstep))
    
    for itraj in range(ntraj):
        path = os.path.join(f"./TRAJ_{itraj + 1:0{index}d}/md/", "NACME")

        with open(path, 'r') as f:
            # skip header and read rest
            line = f.readline()
            line = f.read()
            lines = line.split()

        avg_nacme += abs(np.array([lines[(istate + 1)::(nstate_pair + 1)][:nstep] for istate in range(nstate_pair)], dtype=np.float))

    avg_nacme /= ntraj

    for istep in range(nstep):
        avg_data = f"{istep:8d}"
        for istate in range(nstate_pair):
            avg_data += f"{avg_nacme[istate, istep]:15.8f}"
        f_write += avg_data + "\n"
    typewriter(f_write.rstrip(), "NACME_avg")


def typewriter(string, file_name):
    with open(file_name, "a") as f:
        f.write(string + "\n")


if (__name__ == "__main__"):
    statistical_analysis()

