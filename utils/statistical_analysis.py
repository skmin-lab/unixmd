import math
import os

def main():
    # TODO: motion avg, error message handling?
    ntraj = int(input("Input trajectory number : "))
    index = int(math.log10(ntraj+1))+1

    path = os.path.join(f"./TRAJ_{1:0{index}d}/md/", "MDENERGY")
    f1 = open(path, 'r') 
    lines_1 = f1.readlines()
    f1.close()

    nstep = len(lines_1)-2

    State_avg(ntraj, index, nstep)
    Population_avg(ntraj, index, nstep)
    Coherence_avg(ntraj, index, nstep)
#    Motion_avg(ntraj, index, nstep, parameter_geom)


def State_avg(ntraj, index, nstep):
    w = open("SHSTATE_avg", 'w')
    w.write("#    Step    Average Running State\n")
    avg_state = [0] * (nstep + 1)

    for i in range(ntraj):
        path = os.path.join(f"./TRAJ_{i+1:0{index}d}/md/", "SHSTATE")
        f2 = open(path, 'r')
        lines_2 = f2.readlines()
        f2.close()

        for j in range(nstep + 1):
            tmp = lines_2[j+1].split()
            avg_state[j] += float(tmp[1])

    for i in range(nstep + 1):
        tmp = avg_state[j]/ntraj
        w.write(f"{i:8d}{tmp:15.8f}\n")

    w.close()


def Population_avg(ntraj, index, nstep):
    w = open("BOPOP_avg", 'w')
    w.write("#    Averaged Density Matrix\n")

    for i in range(ntraj):
        path = os.path.join(f"./TRAJ_{i+1:0{index}d}/md/", "BOPOP")
        f2 = open(path, 'r')
        lines_2 = f2.readlines()
        f2.close()
        
        if (i == 0):
            tmp = lines_2[1].split()
            flength = len(tmp) - 1
            avg_state = [[0 for j in range(flength)] for k in range(nstep+1)]

        for j in range(nstep + 1):
            tmp = lines_2[j+1].split()

            for k in range(flength):
                avg_state[j][k] += float(tmp[k+1])

    for i in range(nstep + 1):
        w_str = f"{i:8d}"
        for j in range(flength):
            tmp = avg_state[i][j]/ntraj
            w_str += f"{tmp:15.8f}"
        w.write(w_str + "\n")

    w.close()


def Coherence_avg(ntraj, index, nstep):
    w = open("BOCOH_avg", 'w')
    w.write("#    Averaged Density Matrix: coherence Re-Im\n")

    for i in range(ntraj):
        path = os.path.join(f"./TRAJ_{i+1:0{index}d}/md/", "BOCOH")
        f2 = open(path, 'r')
        lines_2 = f2.readlines()
        f2.close()
        
        if (i == 0):
            tmp = lines_2[1].split()
            flength = len(tmp) - 1
            avg_state = [[0 for j in range(flength)] for k in range(nstep+1)]

        for j in range(nstep + 1):
            tmp = lines_2[j+1].split()

            for k in range(flength):
                avg_state[j][k] += float(tmp[k+1])

    for i in range(nstep + 1):
        w_str = f"{i:8d}"
        for j in range(flength):
            tmp = avg_state[i][j]/ntraj
            w_str += f"{tmp:15.8f}"
        w.write(w_str + "\n")

    w.close()
                                      
def NACME_avg(ntraj, index, nstep):
    w = open("NACME_avg", 'w')
    w.write("#    Averaged Non-Adiabatic Coupling Matrix Eliments: off-diagonal\n")

    for i in range(ntraj):
        path = os.path.join(f"./TRAJ_{i+1:0{index}d}/md/", "NACME")
        f2 = open(path, 'r')
        lines_2 = f2.readlines()
        f2.close()

        if (i == 0):
            tmp = lines_2[1].split()
            flength = len(tmp) - 1
            avg_state = [[0 for j in range(flength)] for k in range(nstep+1)]

        for j in range(nstep + 1):
            tmp = lines_2[j+1].split()

            for k in range(flength):
                avg_state[j][k] += float(tmp[k+1])

    for i in range(nstep + 1):
        w_str = f"{i:8d}"
        for j in range(flength):
            tmp = avg_state[i][j]/ntraj
            w_str += f"{tmp:15.8f}"
        w.write(w_str + "\n")

    w.close()
         
if __name__ == "__main__":
    main()

