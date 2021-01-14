# Rachel Kositsky
# 2021-01-12

# Rules for merging:
# DW “chrom”, “start”, “stop” matches Delly/lumpy chr1:pos1 AND DW “partner_chrom” matches chr2
# OR: DW “chrom”, “start”, “stop” matches Delly/lumpy chr2:pos2 AND DW “partner_chrom” matches chr1

# Print out whatever is in DL, and if in DW, then print that
# add in the callers and Num_callers from DL
# Gene annotation will be in a different script

import argparse
import pandas as pd
import os

def main(args):
    """Use pandas to read in the DELLY/LUMPY (dl) and discowave (dw) tables. 
    Construct the output merged call by using the precise D/L positions as a 
    base, then supplementing them with the most supportive matching DW call."""

    # Read in D/L table and subset pertainent columns
    dl = pd.read_csv(args.DL_file, sep = "\t")

    col_names = ["dave_lab_id", "chr1", "pos1", "chr2", "pos2", "pe", "sr", 
        "pe_sr", "Num_callers", "Callers"]
    out = dl[col_names].copy()

    # Add new columns to output in anticipation of discowave merging
    out[["discowave_chrom", "discowave_start", "discowave_stop", 
        "discowave_pe", "discowave_depth", "discowave_pct_to_partner", 
        "discowave_partner_chrom", "discowave_evenness"]] = "NA"

    # Read in DW file and get padded start/stop locations
    # DW returns a window of possible locations. We'll just extend that window
    # by the padding amount.
    dw = pd.read_csv(args.DW_file, sep = "\t")
    dw["padded_start"] = dw["start"] - args.merge_distance
    dw["padded_stop"] = dw["stop"] + args.merge_distance
    # Add the number of reads supporting the partner chromosome
    dw["pe"] = (dw["all"] * dw["pct_reads_to_partner"] / 100).apply(round)

    # Now merge! Join in any discowave calls.
    for index, row in dl.iterrows():
        # Go through all of DW columns
        # set count to zero for testing if line was found in dw
        dw_count = 0

        # Get matches based on chr1 = DW window and chr2 = partner_chrom
        matching_dw_1 = dw[(dw["chrom"] == row["chr1"]) & 
            (dw["partner_chrom"] == row["chr2"]) &
            (dw["padded_start"] <= row["pos1"]) & 
            (dw["padded_stop"] >= row["pos1"]) &
            (dw["sample"] == row["dave_lab_id"])]

        # Get matches based on chr2 = DW window and chr1 = partner_chrom
        matching_dw_2 = dw[(dw["chrom"] == row["chr2"]) & 
            (dw["partner_chrom"] == row["chr1"]) &
            (dw["padded_start"] <= row["pos2"]) & 
            (dw["padded_stop"] >= row["pos2"]) &
            (dw["sample"] == row["dave_lab_id"])]

        # Put these two together
        matching_dw = matching_dw_1.append(matching_dw_2, ignore_index=True)

        # Add relevant discowave fields if there's a match
        if len(matching_dw) >= 1:

            # Have a warning if there's more than one match.
            if len(matching_dw) > 1:
                print("WARNING: more than one discowave call matches the " \
                    "DELLY/LUMPY call {0}:{1} - {2}:{3}!".format(
                    row["chr1"], row["pos1"], row["chr2"], row["pos2"]))
                print("Here are the {0} matching calls:".format(len(matching_dw)))
                print(matching_dw[["chrom", "start", "stop", "partner_chrom", 
                    "all", "pct_reads_to_partner", "pe"]])

            # Choose the window with the highest paired-end support and 
            # turn it into a Series
            matching_dw = matching_dw.loc[matching_dw["pe"].idxmax(),:]

            # Add relevant fields
            out.loc[index, "discowave_chrom"] = matching_dw["chrom"]
            out.loc[index, "discowave_start"] = matching_dw["start"]
            out.loc[index, "discowave_stop"] = matching_dw["stop"]
            out.loc[index, "discowave_partner_chrom"] = matching_dw["partner_chrom"]
            out.loc[index, "discowave_pe"] = matching_dw["pe"]
            out.loc[index, "discowave_depth"] = matching_dw["all"]
            out.loc[index, "discowave_pct_to_partner"] = matching_dw["pct_reads_to_partner"]
            out.loc[index, "discowave_evenness"] = matching_dw["evenness"]

    # Having gone through each line, now write output
    out.to_csv(args.output_file, sep = "\t", index = False)
    print("Wrote {0} calls to output file {1}".format(len(out), args.output_file))
    return


if __name__ == "__main__":
    # Read in arguments
    parser = argparse.ArgumentParser(description="DELLY/LUMPY and discowave merger")
    parser.add_argument('-dl', '--DL_file', required=True, 
        help="DELLY/LUMPY merged file")
    parser.add_argument('-dw', '--DW_file', required=True, 
        help="discowave file")
    parser.add_argument("-out", "--output_file", required=True, 
        help="Output file location for D/L + DW merged file (tab-separated)")
    parser.add_argument("-dist", "--merge_distance", default=500, type=int,
        help = "Merging distance for DELLY/LUMPY vs discowave annotations")
    
    args = parser.parse_args()

    main(args)

