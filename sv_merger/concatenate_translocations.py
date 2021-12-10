# Rachel Kositsky
# Merges single sample translocation tables row-wise into a matrix

import argparse
import sys

def main(args):
    """Merge single sample translocation tables row-wise into a matrix."""

    n_inputs = len(args.input_files)

    # There must be at least one input file for this program
    if n_inputs == 0:
        raise Exception("No input files provided!")

    # As you read each file in, add its data lines to the output file
    n_lines = 0
    with open(args.output_file, "w") as out_f:
        for i in range(n_inputs):
            print(f"{i+1}/{n_inputs}: reading and writing {args.input_files[i]}")
            with open(args.input_files[i], "r") as in_f:
                # Handling the header
                if i == 0:
                    # Get header from first file
                    first_header = in_f.readline()
                    out_f.write(first_header)
                else:
                    # Compare header line to first header, and don't write it out
                    header = in_f.readline()
                    if header != first_header:
                        raise Exception((f"Header in input file "
                            f"{args.input_files[i]} does not match header of "
                            f"first file {args.input_files[0]}"))
                
                # Write remaining lines to output file
                for line in in_f:
                    out_f.write(line)
                    n_lines += 1

    print(f"Merged {n_inputs} translocation tables and wrote {n_lines} rows "
        f"to {args.output_file}")


def parse_args(args=None):
    """Parse command line arguments and return constants"""
    parser = argparse.ArgumentParser(
        description="Merge multiple single sample translocation tables row-wise " \
        "into a matrix.")

    parser.add_argument("output_file",
        help="Tab-delimited output table of merged structural variants.")

    parser.add_argument("input_files", nargs="+", 
        help="List of tab-delimited input table of structural variants.")

    args = parser.parse_args(args)

    return args


if __name__ == '__main__':
    main(parse_args(sys.argv[1:]))
