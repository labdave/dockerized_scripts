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
    inputs = [0]*n_inputs
    for i in range(n_inputs):
        
        inputs[i] = pd.read_csv(args.input_files[i], sep="\t")
        
        if i == 0:
            col_names = list(inputs[0].columns)
        else:
            # Check that column names match
            assert(list(inputs[i].columns) == col_names)

    # Merge together all inputs row-wise
    df_merged = inputs[0]
    for i in range(1, n_inputs):
        df_merged = pd.concat([df_merged, inputs[i]], axis=0)

    # Write output
    df_merged.to_csv(args.output_file, sep="\t", index = False)

    print((f"Merged {n_inputs} translocation tables and wrote "
           f"{len(df_merged.index)} total rows to {args.output_file}."))


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
