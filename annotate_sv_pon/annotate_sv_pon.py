# Rachel Kositsky
# Created: 2020-07-28
# Updated: 2020-09-13

import argparse
import pandas as pd
import pybedtools
import sys
import os


def get_pon_annotation(row):
    return(",".join([row.chr1, str(row.pos1), row.chr2, str(row.pos2), 
        str(row.n_samples)]))


def get_panel_of_normals_beds(pon_file):
    """Creates two BED files with PON annotations needed for later functions:
    pon.BP1.bed and pon.BP2.bed. """
    pon_df = pd.read_csv(pon_file, sep="\t")

    for idx in ["1", "2"]:
        out_df = pon_df[["chr" + idx, "pos" + idx]]
        out_df.loc[:,"chrom"] = out_df.loc[:,"chr" + idx]
        out_df.loc[:,"start"] = out_df.loc[:,"pos" + idx].map(lambda x: int(x-1))
        out_df.loc[:,"end"] = out_df.loc[:,"pos" + idx].map(lambda x: int(x))
        out_df.loc[:,"annot"] = pon_df.apply(get_pon_annotation, axis=1)
        
        out_df.to_csv("pon.BP{}.bed".format(idx), sep = "\t", index = False, 
            header = False, columns = ["chrom", "start", "end", "annot"])


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
    anno_df.loc[:,anno_col_name] = anno_df.groupby("orig_row")[anno_col_name].transform(lambda x: ";".join(x))
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


def intersect_pon_annotations(row):
    """Function to pass to 'apply' for intersecting panel of normal annotations"""

    bp1 = set(row.BP1.split(";"))
    bp2 = set(row.BP2.split(";"))

    # Get the coordinates which are in common
    common_annotations = bp1.intersection(bp2)

    # For ease of display, don't report coordinates but report the number of samples,
    # which is last in the PON annotation string
    common_annotations = list(map(lambda x: x.split(",")[-1], common_annotations))

    # common_annotations may be [], [""], or a list of number strings, like ["15", "4"]
    if common_annotations in [[], [""]]:
        return("")
    else:
        # If there's something found, report the max
        common_annotations = list(map(int, common_annotations))
        return(str(max(common_annotations)))


def add_panel_of_normal(df, annotation_bed_dict, breakpoint_bed_dict):
    """Add columns for combined breakpoint 1 + 2 annotations with annotation_bed_dict."""

    col_name = "PON"

    tmp_df = df[["chr1"]].copy()
    for idx in ["BP1", "BP2"]:
        # Do intersections with annotation file and save to a file so that we can
        # read in annotations with pandas
        # To save: use moveto instead of saveas to save time and because we're done
        # with using this file's BedTool object. Moves, doesn't copy.
        breakpoint_bed_dict[idx].intersect(
            annotation_bed_dict[idx], wo=True).moveto("{}.bed".format(idx))
        tmp_df = add_collapsed_annotation_df(tmp_df, "{}.bed".format(idx), idx)

    # Make annotation column by combining annotations from BP1 and BP2
    df[col_name] = tmp_df.apply(intersect_pon_annotations, axis=1)

    return(df)


def main(args):
    """Goal: Append panel of normals onto structural variant VCFs"""

    df = pd.read_csv(args.input_file, sep="\t")

    if df.empty:
        df["PON"] = []
        df.to_csv(args.output_file, sep = "\t", index = False)
        return

    # Convert a translocation table to 2 BED files and read in as pybedtools objects
    [bp1, bp2] = get_expanded_bed(df, expansion_distance = 300)
    breakpoint_dict = {"BP1": bp1, "BP2": bp2}
    
    # Load in annotation resources: convert PON table to BEDs and read them in
    get_panel_of_normals_beds(args.pon_file)
    pon_regions = {}
    pon_regions["BP1"] = pybedtools.BedTool("pon.BP1.bed")
    pon_regions["BP2"] = pybedtools.BedTool("pon.BP2.bed")

    df = add_panel_of_normal(df, pon_regions, breakpoint_dict)

    # Save output
    df.to_csv(args.output_file, sep = "\t", index = False)


def parse_args(args=None):
    """Parse command line arguments and return constants"""
    parser = argparse.ArgumentParser(
        description="Annotate structural variant table with whether there are" \
        " translocations picked up in normal samples close to the structural variant.")

    parser.add_argument("input_file",
        help="Tab-delimited input table of structural variants. Requires " \
        "columns chr1,pos1,chr2,pos2 detailing positions of both breakpoints.")

    parser.add_argument("output_file",
        help="Tab-delimited output table of annotated structural variants. " \
        "Input table with 1 new column: PON")


    parser.add_argument("pon_file", help="Panel of normals regions")

    args = parser.parse_args(args)

    return args

if __name__ == '__main__':
    main(parse_args(sys.argv[1:]))

