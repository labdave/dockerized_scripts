# Rachel Kositsky
# Merges single sample translocation tables row-wise into a matrix

import argparse
import pandas as pd
import sys

def main(args):
    """Merge single sample translocation tables row-wise into a matrix."""

    n_inputs = len(args.input_files)

    # There must be at least one input file for this program
    if n_inputs == 0:
        raise Exception("No input files provided!")

    # Read in all inputs and check their columns match up
    inputs = [pd.read_csv(f, sep="\t", low_memory=False) for f in args.input_files]

    # Merge together all inputs row-wise
    df_merged = pd.concat(inputs, axis=0)

    # Write output
    df_merged.to_csv(args.output_file, sep="\t", index = False)

    print("Merged {0} translocation tables and wrote {1} total rows to {2}".format(
        n_inputs, len(df_merged.index), args.output_file))


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
