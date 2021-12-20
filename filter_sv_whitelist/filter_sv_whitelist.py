# Rachel Kositsky
# 2021-12-19

import argparse
import os
import pandas as pd
import sys
import time


def filter_whitelist_by_column_names(header, wl_df):
    """Return a whitelist dataframe filtered for columns matching the input table"""
    wl_df = wl_df[wl_df["column1"].isin(header)]
    wl_df = wl_df[wl_df["column2"].isin(header)]

    return wl_df


def sv_in_whitelist(line, header, whitelist_table):
    """If the line contains a structural variant that's in the whitelist table,
    return the name of the matching whitelist annotation. Otherwise return None."""
    fields = line.strip("\n").split("\t")

    for index, row in whitelist_table.iterrows():
        if ((fields[header.index(row["column1"])] == row["value1"]) and 
            (fields[header.index(row["column2"])] == row["value2"])):
            return row["name"]

    return None


def main(args):
    """Annotate whitelist translocations from structural variant table"""
    
    start_t = time.time()

    # Load in annotation resources
    whitelist_table_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
        "lymphoma_whitelist.tsv")
    whitelist_table = pd.read_csv(whitelist_table_file, sep = "\t")

    n_whitelist_found = 0

    # Input file may be too large for pandas to handle, so we handle it line by line
    with open(args.input_file, "r") as in_f, open(args.output_file, "w") as out_f:
        header = in_f.readline()
        out_f.write("whitelist_annotation\t")
        out_f.write(header)
        # Now that the header's been written over, make it into a searchable list
        header = header.strip("\n").split("\t")

        # Winnow down the whitelist table to just the lines that agree with the 
        # header
        whitelist_table = filter_whitelist_by_column_names(header, whitelist_table)

        # Iterate through the rest of the file and if there are any matches, 
        # write out the line. Prepend the annotation from the whitelist.
        for line in in_f:
            whitelist_annotation = sv_in_whitelist(line, header, whitelist_table)
            if whitelist_annotation != None:
                out_f.write(f"{whitelist_annotation}\t{line}")
                n_whitelist_found += 1
    
    print(f"Wrote {n_whitelist_found} whitelist structural variants over to {args.output_file}")

    print(f"Complete in {time.time() - start_t:.2f} seconds")


def parse_args(args=None):
    """Parse command line arguments and return constants"""
    parser = argparse.ArgumentParser(
        description="Annotate structural variant table with whitelist outputs")

    parser.add_argument("input_file",
        help="Tab-delimited input table of structural variants.")

    parser.add_argument("output_file",
        help="Tab-delimited output table of whitelisted structural variants. " \
        "Input table with 1 new column: whitelist_annotation")

    args = parser.parse_args(args)

    return args

if __name__ == '__main__':
    main(parse_args(sys.argv[1:]))

