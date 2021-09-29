import re
import os
import sys
import pandas as pd
import statistics

def remove_shit(bases):
    # remove end of read character
    bases = bases.replace("$", "")
    # remove start of read character
    new_bases = ""
    next = False
    for i in range(len(bases)):
        if bases[i] == "^":
            next = True
            continue
        if next:
            next = False
            continue
        new_bases += bases[i]
    bases = new_bases
    # remove insertions and deletions
    new_bases = ""
    next = 0
    for i in range(len(bases)):
        if bases[i] == "+" or bases[i] == "-":
            if bases[i+2].isdigit():
                next = 10*int(bases[i+1])+int(bases[i+2])+2
            else:
                next = int(bases[i+1])+1
            continue
        if next:
            next -= 1
            continue
        new_bases += bases[i]
    bases = new_bases

    return bases


def calc_af(base_string):
    af, ad = 0.0, 0.0
    for base in base_string:
        if base != "." and base != "," and base != "*":
            ad += 1
    af = ad/len(base_string)

    return [float("{:.3f}".format(af)), int(ad)]


def calc_std(base_string, pos_string):
    pos_list = []
    for i in range(len(base_string)):
        base = base_string[i]
        pos = pos_string[i]
        if base != "," and base != "." and base != "*":
            pos_list.append(int(pos))

    if len(pos_list) < 2:
        return 0.0
    else:
        return statistics.stdev(pos_list)


mpileup_input_file = sys.argv[1]
clean_vars_file = sys.argv[2]
depth_thresh = int(sys.argv[3])
std_thresh = int(sys.argv[4])
output_file = sys.argv[5]

chrom, pos, af, ad, ad_filter, std, std_filter = [], [], [], [], [], [], []
with open(mpileup_input_file, "r") as f:
    for line in f:
        line_arr = line.strip().split()
        good_bases = remove_shit(line_arr[4])
        af_ad = calc_af(good_bases)
        stdev = calc_std(good_bases, line_arr[6].split(","))

        chrom.append(line_arr[0])
        pos.append(int(line_arr[1]))
        af.append(af_ad[0])
        ad.append(af_ad[1])
        ad_filter.append(1 if af_ad[1] > depth_thresh else 0)
        std.append(stdev)
        std_filter.append(1 if stdev > std_thresh else 0)


mpileup_df = pd.DataFrame({"CHROM": chrom, "POS": pos, "mpileup_af": af, "mpileup_ad": ad, "mpileup_ad_pass": ad_filter, "mpileup_std": std, "mpileup_std_pass": std_filter})

cleaned_df = pd.read_csv(clean_vars_file, sep="\t")
cleaned_df = cleaned_df.merge(mpileup_df, on=["CHROM", "POS"])
cleaned_df.to_csv(output_file, index=False, sep="\t", na_rep="NA")