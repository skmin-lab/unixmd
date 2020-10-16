import math
import os
import shutil

def main():
    #TODO: error message handling?
    filename_1 = input("Input sampled filename (xyz format) : ")
    f1 = open(f'{filename_1}','r')
    lines_1 = f1.readlines()
    f1.close()
    
    natom = int(lines_1[0])
    ntraj = int(len(lines_1)/(natom+2))

    print(f"Total atom number : {natom}\n")
    print(f"Total trajectory number : {ntraj}\n")

    filename_2 = input("Input setting filename (run.py format) : ")
    f2 = open(f'{filename_2}','r')
    lines_2 = f2.readlines()
    f2.close()

    for i in range(len(lines_2)):
        tmp = lines_2[i].split()
        if (len(tmp) != 0):
            if (tmp[0] == "mol"):
                saveline1 = i

            for j in range(len(tmp)):
                dir_index = tmp[j].find("md.run")      
                if (dir_index != -1):
                    saveline2 = i

    index = int(math.log10(ntraj+1))+1

    for i in range(ntraj):
        w = open('tmp_runpi', 'w')
        w.write("from molecule import Molecule\n")
        w.write("import qm, mqc\n")
        w.write("from thermostat import *\n")
        w.write("from misc import data\n\n")
        w.write("geom = \"\"\"\n")

        for j in range(natom+2):
            nline = i * (natom+2) + j
            w.write(lines_1[nline])
        w.write("\"\"\"\n")

        for j in range(saveline2-saveline1):
            w.write(lines_2[saveline1+j])

        savedat = lines_2[saveline2].split()

        for j in range(len(savedat)):
            dir_index = savedat[j].find("input_dir")      
            if (dir_index != -1):
                savedat[j] = f"input_dir=\"./TRAJ_{i+1:0{index}d}\", " 
            w.write(savedat[j] + " ")

        os.mkdir(f'TRAJ_{i+1:0{index}d}')
        shutil.move("tmp_runpi", f"TRAJ_{i+1:0{index}d}")            
        os.rename(f"TRAJ_{i+1:0{index}d}/tmp_runpi", f"TRAJ_{i+1:0{index}d}/run.py")

        w.close()
         
if __name__ == "__main__":
    main()

