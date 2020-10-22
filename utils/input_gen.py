import os
import argparse
import shutil

def main():
    #TODO: error message handling?
    parser = argparse.ArgumentParser(description = "PYTHON script for UNI-xMD input generator")
    parser.add_argument('-s', action='store', dest='sample_XYZ', type=str, help="Directory including sampled files, sampled files must be written in xyz format", required=True)
    parser.add_argument('-r', action='store', dest='template_run_py', type=str, help="Template for setting filename, must be written in run.py format", required=True)
    args = parser.parse_args()
   
    path = os.path.join(f"./", f"{args.sample_XYZ}")
    ntraj = len([itraj for itraj in os.listdir(path) if itraj.startswith("TRAJ_")])
    str_traj = str(ntraj+1) 
    index = len(str_traj) # number for trajectory indexing
    
    print(f"Total trajectory number : {ntraj}\n")

    readfile = open(f'{args.template_run_py}','r')
    lines_runpy = readfile.readlines()
    readfile.close()

    for i in range(len(lines_runpy)):
        tmp = lines_runpy[i].split()
        if (len(tmp) != 0):
            if (tmp[0] == "mol"):
                saveline1 = i

            for j in range(len(tmp)):
                dir_index = tmp[j].find("md.run")      
                if (dir_index != -1):
                    saveline2 = i

    for itraj in range(ntraj):
        w = open('tmp_runpi', 'w')
        w.close()

        with open('tmp_runpi', 'a') as f:
            f.write("from molecule import Molecule\n")
            f.write("import qm, mqc\n")
            f.write("from thermostat import *\n")
            f.write("from misc import data\n\n")
            f.write("geom = \"\"\"\n")

        traj_path = os.path.join(path,f"TRAJ_{itraj+1:0{index}d}.xyz")
        readfile = open(f'{traj_path}','r')
        xyz = readfile.readlines()
        readfile.close()

        if (itraj == 0):
            natom = int(xyz[0])
            print(f"Total atom number : {natom}\n")

        for jatom in range(natom+2):
            with open('tmp_runpi', 'a') as f:
                f.write(xyz[jatom])

        with open('tmp_runpi', 'a') as f:
            f.write("\"\"\"\n")

        for j in range(saveline2-saveline1):
            with open('tmp_runpi', 'a') as f:
                f.write(lines_runpy[saveline1+j])

        savedat = lines_runpy[saveline2].split()

        for j in range(len(savedat)):
            dir_index = savedat[j].find("input_dir")      
            if (dir_index != -1):
#                savedat[j] = f"input_dir=\"./TRAJ_{i+1:0{index}d}\", " 
                savedat[j] = f"input_dir=\"./\", " 
            with open('tmp_runpi', 'a') as f:
                f.write(savedat[j] + " ")

        os.mkdir(f"TRAJ_{itraj+1:0{index}d}")
        shutil.move("tmp_runpi", f"TRAJ_{itraj+1:0{index}d}")            
        os.rename(f"TRAJ_{itraj+1:0{index}d}/tmp_runpi", f"TRAJ_{itraj+1:0{index}d}/run.py")
         
if __name__ == "__main__":
    main()

