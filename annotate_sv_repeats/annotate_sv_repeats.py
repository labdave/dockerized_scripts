# Rachel Kositsky
# 2020-07-21

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


def intersect_repeat_masker(row):
    """Function to pass to 'apply' for intersecting repeat masker annotations"""

    rm_1 = set(row.BP1_repeats_600bp.split(","))
    rm_2 = set(row.BP2_repeats_600bp.split(","))

    # TODO: add a step here that converts various repeat types to their more 
    # general classes, e.g. AluSc -> Alu

    return(",".join(rm_1.intersection(rm_2)))


def intersect_segmental_duplications(row):
    """Function to pass to 'apply' for intersecting repeat masker annotations"""

    segdup_1 = set(row.BP1_segdup_600bp.split(","))
    segdup_2 = set(row.BP2_segdup_600bp.split(","))

    # TODO: add a step here that *matches* the segmental duplication annotations
    # e.g. chr1 and chr8 -> they match

    return(",".join(segdup_1.intersection(segdup_2)))


def merge_repeat_masker(df):
    """Make an additional column describing the merged repeat masker info"""
    
    for idx in ["1", "2"]:
        rm_df = pd.read_csv("rm.bp{}.bed".format(idx), sep = "\t",
            names = ["chrom", "start", "end", "orig_row", "rm_chrom", 
            "rm_start", "rm_end", "rm_anno", "size", "strand", "overlap"])
        
        # Collapse repeat masker annotations by orig_row
        rm_df = rm_df[["orig_row", "rm_anno"]]
        rm_df.loc[:,"rm_anno"] = rm_df.groupby("orig_row")["rm_anno"].transform(lambda x: ",".join(x))
        rm_df = rm_df.drop_duplicates()

        # Set up joining
        col_name =  "BP{}_repeats_600bp".format(idx)
        # set index to the original row for index-on-index joining later
        rm_df.index = rm_df["orig_row"]

        # Rename
        rm_df.loc[:,col_name] = rm_df.loc[:,"rm_anno"]
        
        # Turn into a Series
        rm_df = rm_df.loc[:,col_name]

        # Join in annotations into normal dataframe
        df = df.join(other = rm_df, how = "left")

        # Fill in NAs with empty strings
        df[col_name] = df[col_name].fillna("")
        
    return(df)


def merge_segmental_duplications(df):
    """Make an additional column describing the merged segdup info"""
    # TODO: change so that it's the matching segdups
    
    for idx in ["1", "2"]:
        segdup_df = pd.read_csv("segdup.bp{}.bed".format(idx), sep = "\t",
            names = ["chrom", "start", "end", "orig_row", "segdup_chrom", 
            "segdup_start", "segdup_end", "segdup_anno", "size", "strand", 
            "overlap"])
        
        # Collapse segmental duplication annotations by orig_row
        segdup_df = segdup_df[["orig_row", "segdup_anno"]]
        segdup_df.loc[:,"segdup_anno"] = segdup_df.groupby("orig_row")["segdup_anno"].transform(lambda x: ",".join(x))
        segdup_df = segdup_df.drop_duplicates()

        # Set up joining
        col_name =  "BP{}_segdup_600bp".format(idx)
        # set index to the original row for index-on-index joining later
        segdup_df.index = segdup_df["orig_row"]

        # Rename
        segdup_df.loc[:,col_name] = segdup_df.loc[:,"segdup_anno"]
        
        # Turn into a Series
        segdup_df = segdup_df.loc[:,col_name]

        # Join in annotations into normal dataframe
        df = df.join(other = segdup_df, how = "left")

        # Fill in NAs with empty strings
        df[col_name] = df[col_name].fillna("")
        
    return(df)


def main(args):
    """Goal: Append repeat masker and segdup onto structural variant VCFs"""

    df = pd.read_csv(args.input_file, sep="\t")

    # Convert a translocation table to 2 BED files and read in as pybedtools objects
    [bp1_600, bp2_600] = get_expanded_bed(df, expansion_distance = 600)
    
    # Load in the resources
    repeat_masker = pybedtools.BedTool(args.repeat_masker_bed)
    segmental_dups = pybedtools.BedTool(args.segmental_duplication_bed)

    # Do intersections with repeat masker and save to a file so that we can
    # read in annotations with pandas
    # To save: use moveto instead of saveas to save time and because we're done
    # with using this file's BedTool object. Moves, doesn't copy.
    bp1_600.intersect(repeat_masker, wo=True).moveto("rm.bp1.bed")
    bp2_600.intersect(repeat_masker, wo=True).moveto("rm.bp2.bed")    

    # Merge in repeat masker and get matching column
    df = merge_repeat_masker(df)
    df["matching_repeats"] = df.apply(intersect_repeat_masker, axis=1)

    # Merge in segmental duplications
    bp1_600.intersect(segmental_dups, wo=True).moveto("segdup.bp1.bed")
    bp2_600.intersect(segmental_dups, wo=True).moveto("segdup.bp2.bed")


    # Add segmental duplication columns and get matching column
    df = merge_segmental_duplications(df)
    df["matching_segdup"] = df.apply(intersect_segmental_duplications, axis=1)

    # Save output
    df.to_csv(args.output_file, sep = "\t", index = False)



def parse_args(args=None):
    """Parse command line arguments and return constants"""
    parser = argparse.ArgumentParser(
        description="Annotate structural variant table with repeat masker and" \
        " segmental duplication columns. BP1/2 repeats and segdup columns " \
        " denote whether individual breakpoints have repeat masker/" \
        "polynucleotide/segmental duplication regions within 600bp. " \
        "matching_repeats: list of any repeat regions that match in BP1 and BP2." \
        " matching_segdup: T/F whether BP1 and BP2 have complementary segmental duplication regions."
        )

    parser.add_argument("input_file",
        help="Tab-delimited input table of structural variants. Requires " \
        "columns chr1,pos1,chr2,pos2 detailing positions of both breakpoints.")

    parser.add_argument("output_file",
        help="Tab-delimited output table of annotated structural variants. " \
        "Input table with 8 new columns: BP1_repeats_600bp, BP2_repeats_600bp, " \
        "matching_repeats, BP1_poly_nt_100bp, BP2_poly_nt_100bp, " \
        "BP1_segdup_600bp, BP2_segdup_600bp, matching_segdup")

    parser.add_argument("repeat_masker_bed",
        help="BED6 file of repeat masker regions, e.g. downloaded from UCSC Table Browser")

    parser.add_argument("segmental_duplication_bed",
        help="BED6 file of segmental duplication regions, e.g. downloaded from UCSC Table Browser")

    args = parser.parse_args(args)

    return args

if __name__ == '__main__':
    main(parse_args(sys.argv[1:]))

