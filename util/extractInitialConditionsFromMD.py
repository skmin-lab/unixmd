import argparse
import re

def parse_md_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    geometries = {}
    current_geom = []
    current_iter = None

    for line in lines:
        iter_match = re.match(r'^\s*(?:D|MD) iter:\s*(\d+)', line)
        if iter_match:
            if current_iter is not None and current_geom:
                geometries[current_iter] = current_geom
            current_iter = int(iter_match.group(1))
            current_geom = []
        elif re.match(r'^\s*[A-Z]', line):
            tokens = line.split()
            atom = tokens[0]
            x, y, z = map(float, tokens[1:4])
            extras = list(map(float, tokens[-3:]))
            current_geom.append((atom, x, y, z, extras))

    if current_iter is not None and current_geom:
        geometries[current_iter] = current_geom

    return geometries

def save_as_xyz(geometry, filename, iter_label):
    with open(filename, 'w') as f:
        f.write(f"{len(geometry)}\n")
        f.write(f"Geometry from iter {iter_label} (x y z + extras)\n")
        for atom, x, y, z, extras in geometry:
            extra_str = " ".join(f"{val:16.8f}" for val in extras)
            f.write(f"{atom} {x:16.8f} {y:16.8f} {z:16.8f} {extra_str}\n")
    print(f"Saved {filename}")

def main():
    parser = argparse.ArgumentParser(description="Extract geometries from DFTB+ MD output file.")
    parser.add_argument("input_file", help="Path to the input file.")
    parser.add_argument("--start", type=int, required=True, help="Start iteration number.")
    parser.add_argument("--end", type=int, required=True, help="End iteration number (inclusive).")
    parser.add_argument("--step", type=int, default=1, help="Step between iterations (default: 1)")
    parser.add_argument("--renumber", type=int, help="Start renumbering output files from this number")
    parser.add_argument("--width", type=int, help="Width of the integer field in generated file name")

    args = parser.parse_args()

    geometries = parse_md_file(args.input_file)
    output_index = args.renumber if args.renumber is not None else None
    total_number = (args.end-args.start)/args.step

    for i in range(args.start, args.end + 1, args.step):
        if i in geometries:
            if output_index is not None:
                filename = f"sample_{output_index:0{args.width}d}.xyz"
                save_as_xyz(geometries[i], filename, iter_label=i)
                output_index += 1
            else:
                filename = f"sample_{i:0{args.width}d}.xyz"
                save_as_xyz(geometries[i], filename, iter_label=i)
        else:
            print(f"Warning: Iteration {i} not found in file.")

if __name__ == "__main__":
    main()
