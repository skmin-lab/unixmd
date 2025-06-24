import argparse
import re

def extract_initial_conditions_from_MD():
    """ Python utility script for PyUNIxMD output analysis
        In this script, PyUNIxMD 'MOVIE.xyz' output file is post-processed to create initial conditions (positions, momenta)
    """
    parser = argparse.ArgumentParser(description="Extract geometries from PyUNIxMD 'MOVIE.xyz' output file.")
    parser.add_argument("-i", "--input", action='store', dest='input', type=str, \
        help="Path to the input file.", default='./MOVIE.xyz')
    parser.add_argument("-s", "--start", action='store', dest='start', type=int, \
        help="Start iteration number.", required=True)
    parser.add_argument("-e", "--end", action='store', dest='end', type=int, \
        help="End iteration number (inclusive).", required=True)
    parser.add_argument("-d", "--step", action='store', dest='step', type=int, \
        help="Step between iterations (default: 1)", default=1)
    args = parser.parse_args()

    geometries = parse_md_file(args.input)
    ntrajs = len(range(args.start, args.end + 1, args.step))
    index = len(str(ntrajs))

    for file_index, step in enumerate(range(args.start, args.end + 1, args.step), start=1):
        if (step in geometries):
            filename = f"sample_{file_index:0{index}d}.xyz"
            save_as_xyz(geometries[step], filename, step)
        else:
            print (f"\n\n Warning: Iteration {step} not found in file. \n\n", flush=True)

def parse_md_file(file_path):
    """ Read MOVIE.xyz and extract xyz coordinates and extras
    """ 
    with open(file_path, 'r') as f:
        lines = f.readlines()

    geometries = {}
    current_geom = []
    current_iter = None

    for line in lines:
        iter_match = re.match(r'^\s*Step:\s*(\d+)', line)
        if (iter_match):
            if ((current_iter != None) and (current_geom)):
                geometries[current_iter] = current_geom
            current_iter = int(iter_match.group(1))
            current_geom = []
        elif (re.match(r'^\s*[A-Z]', line)):
            tokens = line.split()
            atom = tokens[0]
            x, y, z = map(float, tokens[1:4])
            extras = list(map(float, tokens[-3:]))
            current_geom.append((atom, x, y, z, extras))

    if ((current_iter !=  None) and (current_geom)):
        geometries[current_iter] = current_geom

    return geometries

def save_as_xyz(geometry, filename, step):
    """ Save the conditions of samplings as xyz format
    """
    with open(filename, 'w') as f:
        f.write(f"{len(geometry)}\n")
        f.write(f"Geometry from iter {step} (x y z + extras)\n")
        for atom, x, y, z, extras in geometry:
            extra_str = " ".join(f"{val:15.8f}" for val in extras)
            f.write(f"{atom} {x:15.8f} {y:15.8f} {z:15.8f} {extra_str}\n")
    print (f"Saved {filename}", flush=True)

if __name__ == "__main__":
    extract_initial_conditions_from_MD()
