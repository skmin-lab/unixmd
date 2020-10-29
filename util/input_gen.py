import os
import argparse
import shutil

def main():
    parser = argparse.ArgumentParser(description="PYTHON script for UNI-xMD input generator")
    parser.add_argument('-sample', action='store', dest='sample_XYZ', type=str, default="TRAJ", help="Directory including sampled files, sampled files must be written in xyz format")
    parser.add_argument('-runfile', action='store', dest='template_run_py', type=str, default="run.py", help="Template for setting filename, must be written in run.py format")
    parser.add_argument('-ntraj', action='store', dest='ntraj', type=int, help="Total trajectory number", required=True)
    args = parser.parse_args()
   
    str_traj = str(args.ntraj) 
    index = len(str_traj) # number for trajectory indexing
    
    print (f"Total trajectory number : {args.ntraj}\n", flush=True)

    for itraj in range(args.ntraj):
        os.mkdir(f"TRAJ_{itraj+1:0{index}d}")
        shutil.copy(f"{args.template_run_py}", f"TRAJ_{itraj+1:0{index}d}/run.py")
        shutil.copy(f"{args.sample_XYZ}/TRAJ_{itraj+1:0{index}d}.xyz", f"TRAJ_{itraj+1:0{index}d}/geom.xyz")
         
if __name__ == "__main__":
    main()

