# Rachel Kositsky
# 2020-07-21

import argparse
import pandas as pd
import pybedtools
import re
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


def intersect_repeat_masker(row):
    """Function to pass to 'apply' for intersecting repeat masker annotations"""

    rm_1 = set(row.BP1_repeats_200bp.split(","))
    rm_2 = set(row.BP2_repeats_200bp.split(","))

    # TODO: add a step here that converts various repeat types to their more 
    # general classes, e.g. AluSc -> Alu
    # Would add an argument "level": "exact", "subfamily", or "family"

    return(",".join(rm_1.intersection(rm_2)))


def intersect_segmental_duplications(row, match_distance):
    """Function to pass to 'apply' for intersecting repeat masker annotations"""

    segdup_1 = row.BP1_segdup_200bp.split(",")
    segdup_2 = row.BP2_segdup_200bp.split(",")

    # TODO: add a step here that *matches* the segmental duplication annotations
    # e.g. chr1 and chr8 -> they match
    # Get regions which appear in both annotations
    segdup_matches = set(segdup_1).intersection(segdup_2)

    ## Add on regions in breakpoint 1 that are close to breakpoint 2
    # Filter by chromosome 2
    matches_other_chrom = filter(lambda x: (row.chr2 + ":") in x, segdup_1)
    # Filter by position 2
    matches_other_position = list(filter(
        lambda x: abs(row.pos2 - int(x.split(":")[1])) < match_distance, 
        matches_other_chrom))
    # Add on to overall matches
    segdup_matches.update(matches_other_position)

    ## Add on regions in breakpoint 2 that are close to breakpoint 1
    # Filter by chromosome 1
    matches_other_chrom = filter(lambda x: (row.chr1 + ":") in x, segdup_2)
    # Filter by position 1
    matches_other_position = list(filter(
        lambda x: abs(row.pos1 - int(x.split(":")[1])) < match_distance, 
        matches_other_chrom))
    # Add on to overall matches
    segdup_matches.update(matches_other_position)

    return(",".join(segdup_matches))


def add_collapsed_annotation_df(df, intersect_file_name, anno_col_name):
    """Helper function: returns dataframe with an added column named 
    annot_col_name with collapsed annotations from intersect_file_name."""

    anno_df = pd.read_csv(intersect_file_name, sep = "\t",
        names = ["orig_chrom", "orig_start", "orig_end", "orig_row", "anno_chrom", 
        "anno_start", "anno_end", anno_col_name, "unknown", "strand", "overlap"])
    
    # Fill in with empty strings if nothing came of the intersection
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
    anno_df = anno_df.loc[:,anno_col_name]

    # Join in annotations into normal dataframe
    df = df.join(other = anno_df, how = "left")

    # Fill in NAs with empty strings
    df[anno_col_name] = df[anno_col_name].fillna("")

    return(df)


def add_repeat_masker(df):
    """Make an additional column describing the merged repeat masker info"""
    
    for idx in ["1", "2"]:
        df = add_collapsed_annotation_df(df, "rm.bp{}.200.bed".format(idx), 
            "BP{}_repeats_200bp".format(idx))
        
    return(df)


def add_segmental_duplications(df):
    """Make an additional column describing the merged segdup info"""
    
    for idx in ["1", "2"]:
        df = add_collapsed_annotation_df(df, "segdup.bp{}.200.bed".format(idx), 
            "BP{}_segdup_200bp".format(idx))

    return(df)


def extract_polynt(rm_annotation):
    """Returns string of repeat annotations with just polynucleotides
    Example: '(A)n,(TTA)n,(T)n' -> ['(A)n', '(AC)n', '(T)n'] -> '(A)n,(T)n'"""
    rm_list = rm_annotation.split(",")
    polynt_list = list(filter(lambda x: bool(re.match(r"^\([A,G,C,T]\)n$",x)), 
        rm_list))

    return(",".join(polynt_list))


def add_polynucleotides(df):
    """Add columns for polynucleotide for BP1, BP2. Relies on repeat masker 
    columns."""
    for idx in ["1", "2"]:
        polynt_col_name = "BP{}_polynt_200bp".format(idx)
        df = add_collapsed_annotation_df(df, "rm.bp{}.200.bed".format(idx), 
            polynt_col_name)

        # Modify df column to only contain poly nucleotide repeats
        df[polynt_col_name] = list(map(extract_polynt, df[polynt_col_name]))

    return(df)


def main(args):
    """Goal: Append repeat masker and segdup onto structural variant VCFs"""

    df = pd.read_csv(args.input_file, sep="\t")

    # Convert a translocation table to 2 BED files and read in as pybedtools objects
    [bp1_200, bp2_200] = get_expanded_bed(df, expansion_distance = 200)
    
    # Load in the resources
    repeat_masker = pybedtools.BedTool(args.repeat_masker_bed)
    segmental_dups = pybedtools.BedTool(args.segmental_duplication_bed)

    # Do intersections with repeat masker and save to a file so that we can
    # read in annotations with pandas
    # To save: use moveto instead of saveas to save time and because we're done
    # with using this file's BedTool object. Moves, doesn't copy.
    bp1_200.intersect(repeat_masker, wo=True).moveto("rm.bp1.200.bed")
    bp2_200.intersect(repeat_masker, wo=True).moveto("rm.bp2.200.bed")

    # Merge in repeat masker and get matching column
    df = add_repeat_masker(df)
    df["matching_repeats"] = df.apply(intersect_repeat_masker, axis=1)

    # Add polynucleotide columns
    df = add_polynucleotides(df)

    # Merge in segmental duplications
    bp1_200.intersect(segmental_dups, wo=True).moveto("segdup.bp1.200.bed")
    bp2_200.intersect(segmental_dups, wo=True).moveto("segdup.bp2.200.bed")

    # Add segmental duplication columns and get matching column
    df = add_segmental_duplications(df)
    df["segdup_100k"] = df.apply(
        lambda row: intersect_segmental_duplications(row,100000), axis=1)
    df["segdup_1M"] = df.apply(
        lambda row: intersect_segmental_duplications(row,1000000), axis=1)
    df["segdup_10M"] = df.apply(
        lambda row: intersect_segmental_duplications(row,10000000), axis=1)
    df["segdup_100M"] = df.apply(
        lambda row: intersect_segmental_duplications(row,100000000), axis=1)

    # Save output
    df.to_csv(args.output_file, sep = "\t", index = False)



def parse_args(args=None):
    """Parse command line arguments and return constants"""
    parser = argparse.ArgumentParser(
        description="Annotate structural variant table with repeat masker and" \
        " segmental duplication columns. BP1/2 repeats and segdup columns " \
        " denote whether individual breakpoints have repeat masker/" \
        "polynucleotide/segmental duplication regions within 200bp. " \
        "matching_repeats: list of any repeat regions that match in BP1 and BP2." \
        " matching_segdup: T/F whether BP1 and BP2 have complementary segmental duplication regions."
        )

    parser.add_argument("input_file",
        help="Tab-delimited input table of structural variants. Requires " \
        "columns chr1,pos1,chr2,pos2 detailing positions of both breakpoints.")

    parser.add_argument("output_file",
        help="Tab-delimited output table of annotated structural variants. " \
        "Input table with 11 new columns: BP1_repeats_200bp,BP2_repeats_200bp," \
        "matching_repeats,BP1_polynt_200bp,BP2_polynt_200bp,BP1_segdup_200bp," \
        "BP2_segdup_200bp,segdup_100k,segdup_1M,segdup_10M,segdup_100M")

    parser.add_argument("repeat_masker_bed",
        help="BED6 file of repeat masker regions, e.g. downloaded from UCSC Table Browser")

    parser.add_argument("segmental_duplication_bed",
        help="BED6 file of segmental duplication regions, e.g. downloaded from UCSC Table Browser")

    args = parser.parse_args(args)

    return args

if __name__ == '__main__':
    main(parse_args(sys.argv[1:]))

