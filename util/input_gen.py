import os
import argparse
import shutil

def input_gen():
    """ Python utility script for PyUNIxMD input generator
        In this script, input trajectories are generated from sampled geometry files and running script
        WARNING: Sampled geom files must be named in "sample_(number).xyz"
                 running script must be read geom from "geom.xyz" file
    """
    parser = argparse.ArgumentParser(description="Python script for PyUNIxMD input generator", \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '-dir', action='store', dest='sample_dir', type=str, default="Sampled", \
        help="Directory name of sampled files, sampled files must be written in extended xyz format.")
    parser.add_argument('-f', '-file', action='store', dest='running_script', type=str, default="run.py", \
        help="Filename of personal running script, must be written in running script format. \
        The geometry section of the running script must be read from 'geom.xyz' file.")
    parser.add_argument('-n', '-ntrajs', action='store', dest='ntrajs', type=int, \
        help="Total number of trajectories", required=True)
    args = parser.parse_args()

    # number for trajectory indexing
    index = len(str(args.ntrajs))

    print (f"Total number of trajectories : {args.ntrajs}\n", flush=True)

    # copy from each prepared files
    for itraj in range(args.ntrajs):
        os.mkdir(f"TRAJ_{itraj + 1:0{index}d}")
        shutil.copy(f"{args.running_script}", f"TRAJ_{itraj + 1:0{index}d}/run.py")
        shutil.copy(f"{args.sample_dir}/sample_{itraj + 1:0{index}d}.xyz", f"TRAJ_{itraj + 1:0{index}d}/geom.xyz")

if (__name__ == "__main__"):
    input_gen()

