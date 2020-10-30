import os
import argparse
import shutil

def main():
    """ Python utility script for UNI-xMD input generator
        In this script, input trajectories are generated from sampled geometry files and running script
        WARNING: Sampled geom files must be named in "sample_(number).xyz"
                 running script must be read geom from "geom.xyz" file
    """
    parser = argparse.ArgumentParser(description="Python script for UNI-xMD input generator")
    parser.add_argument('-d', action='store', dest='sample_dir', type=str, default="Sampled", \
        help="Directory including sampled files, sampled files must be written in extended xyz format")
    parser.add_argument('-f', action='store', dest='running_script', type=str, default="run.py", \
        help="Filename of personal running script, must be written in running script format. \
        The geometry section of the running script must be read from 'geom.xyz'")
    parser.add_argument('-n', action='store', dest='ntraj', type=int, \
        help="Total trajectory number", required=True)
    args = parser.parse_args()
 
    # number for trajectory indexing
    index = len(str(args.ntraj))

    print (f"Total trajectory number : {args.ntraj}\n", flush=True)

    # copy from each prepared files
    for itraj in range(args.ntraj):
        os.mkdir(f"TRAJ_{itraj + 1:0{index}d}")
        shutil.copy(f"{args.running_script}", f"TRAJ_{itraj + 1:0{index}d}/run.py")
        shutil.copy(f"{args.sample_dir}/sample_{itraj + 1:0{index}d}.xyz", "TRAJ_{itraj + 1:0{index}d}/geom.xyz")

if (__name__ == "__main__"):
    main()

