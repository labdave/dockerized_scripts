# Rachel Kositsky
# 2021-01-28

import argparse
import pandas as pd
import pybedtools
import sys
import os


def get_expanded_bed(df, expansion_distance):
    # Return BedTool objects for breakpoint 1 and breakpoint 2 padded by the 
    # expansion distance

    return_list = [None, None]

    for idx in ["1", "2"]:
        bp_df = df[["chr" + idx, "pos" + idx]]
        bp_df.loc[:,"chrom"] = bp_df.loc[:,"chr" + idx]

        bp_df.loc[:,"start"] = bp_df.loc[:,"pos" + idx].map(lambda x: int(x - expansion_distance/2))
        bp_df.loc[:,"end"] = bp_df.loc[:,"pos" + idx].map(lambda x: int(x + expansion_distance/2))
        bp_df.loc[:,"orig_row"] = bp_df.index

        file_name = "bp{}.{}.bed".format(idx, expansion_distance)
        bp_df[["chrom", "start", "end", "orig_row"]].to_csv(
            file_name, sep="\t", header=False, index=False)

        # Read in these BED files as pybedtools objects
        return_list[int(idx)-1] = pybedtools.BedTool(file_name)

    return(return_list)

def add_collapsed_annotation_df(df, intersect_file_name, anno_col_name):
    """Helper function: returns dataframe with an added column named 
    annot_col_name with collapsed annotations from intersect_file_name."""

    anno_df = pd.read_csv(intersect_file_name, sep = "\t",
        names = ["orig_chrom", "orig_start", "orig_end", "orig_row", "anno_chrom", 
        "anno_start", "anno_end", anno_col_name, "unknown", "strand", "overlap"])
    
    if anno_df.empty:
        df[anno_col_name] = ""
        return(df)

    # Collapse annotations by orig_row
    anno_df = anno_df[["orig_row", anno_col_name]]
    anno_df.loc[:,anno_col_name] = anno_df.groupby("orig_row")[anno_col_name].transform(lambda x: ",".join(x))
    anno_df = anno_df.drop_duplicates()

    # Set up joining
    # Set index to the original row for index-on-index joining
    anno_df.index = anno_df["orig_row"]
    
    # Turn into a Series
    anno_df = anno_df[anno_col_name]

    # Join in annotations into normal dataframe
    df = df.join(other = anno_df, how = "left")

    # Fill in NAs with empty strings
    df[anno_col_name] = df[anno_col_name].fillna("")

    return(df)

def add_combined_annotations(df, annotation_bed_object, col_name, 
    breakpoint_bed_dict, breakpoint_idx=["BP1", "BP2"]):
    """Add columns for combined breakpoint 1 + 2 annotations with annotation_bed_object."""

    tmp_df = df[["chr1"]].copy()
    for idx in breakpoint_idx:
        # Do intersections with annotation file and save to a file so that we can
        # read in annotations with pandas
        # To save: use moveto instead of saveas to save time and because we're done
        # with using this file's BedTool object. Moves, doesn't copy.
        breakpoint_bed_dict[idx].intersect(
            annotation_bed_object, wo=True).moveto("{}.bed".format(idx))
        tmp_df = add_collapsed_annotation_df(tmp_df, "{}.bed".format(idx), idx)

    # Make annotation column by combining annotations from BP1 and BP2
    df[col_name] = tmp_df[breakpoint_idx].apply(lambda x: ",".join(x), axis=1)
    # Remove extra commas
    df[col_name] = df[col_name].apply(lambda x: x.strip(","))

    return(df)


def main(args):
    """Goal: Append gene names onto structural variant VCFs
    Add these columns:
    BP1_Gene,BP2_Gene,BP1_Dist,BP2_Dist"""

    df = pd.read_csv(args.input_file, sep="\t")

    # Convert a translocation table to 2 BED files and read in as pybedtools objects
    [bp1, bp2] = get_expanded_bed(df, expansion_distance = 1)
    breakpoint_dict = {"BP1": bp1, "BP2": bp2}
    
    # Load in annotation resources
    gene_regions = pybedtools.BedTool(args.gene_bed)

    df = add_combined_annotations(df, gene_regions, "BP1_gene", breakpoint_dict, ["BP1"])
    df = add_combined_annotations(df, gene_regions, "BP2_gene", breakpoint_dict, ["BP2"])

    # Save output
    df.to_csv(args.output_file, sep = "\t", index = False)


def parse_args(args=None):
    """Parse command line arguments and return constants"""
    parser = argparse.ArgumentParser(
        description="Annotate structural variant table with gene names.")

    parser.add_argument("input_file",
        help="Tab-delimited input table of structural variants. Requires " \
        "columns chr1,pos1,chr2,pos2 detailing positions of both breakpoints.")

    parser.add_argument("output_file",
        help="Tab-delimited output table of annotated structural variants. " \
        "Input table with 2 new columns: BP1_gene, BP2_gene")

    parser.add_argument("gene_bed",
        help="Annotated BED file with gene regions and names")

    args = parser.parse_args(args)

    return args

if __name__ == '__main__':
    main(parse_args(sys.argv[1:]))

