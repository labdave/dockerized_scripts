import re
import os
import sys
import pandas as pd
import statistics

def get_ref_alt(chrom_pos):
    if chrom_pos not in ref_alt_dict:
        return ["-", "-"]
    return ref_alt_dict[chrom_pos]

def remove_shit(bases, indel):
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

    if not indel:
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


def calc_af(dp, base_string, ref, alt):

    if int(dp) == 0:
        return [float("{:.3f}".format(af)), int(ad)]

    base_string = base_string.lower()
    dp = int(dp)

    # snv
    if len(ref) == 1 and len(alt) == 1:
        af, ad = 0.0, 0.0
        for base in base_string:
            if base == alt.lower():
                ad += 1
        af = ad/dp

        return [float("{:.3f}".format(af)), int(ad)]

    # indel
    if len(alt) != 1:
        search_string = "+"+str(len(alt)-1)+alt[1:].lower()

    if len(ref) != 1:
        search_string = "-"+str(len(ref)-1)+alt[1:].lower()

    ad = base_string.lower().count(search_string)
    af = ad/dp

    return [float("{:.3f}".format(af)), int(ad)]


def calc_std(base_string, pos_string, ref, alt):

    base_string = base_string.lower()
    pos_list = []

    # snv
    if len(ref) == 1 and len(alt) == 1:
        for i in range(len(base_string)):
            base = base_string[i]
            pos = pos_string[i]
            if base == alt.lower():
                pos_list.append(int(pos))

    # indel
    else:
        base_string = re.sub("[\\.,][-\\+][0-9]+[ACGTNacgtn*#]+", "@", base_string)
        for i in range(len(base_string)):
            base = base_string[i]
            pos = pos_string[i]
            if base == "@":
                pos_list.append(int(pos))

    if len(pos_list) < 2:
        return 0.0
    else:
        return statistics.stdev(pos_list)


def calc_alt_bases(base_string, ref, alt):
    count_list = []

    # snv
    if len(ref) == 1 and len(alt) == 1:
        count_list.append(base_string.lower().count("a"))
        count_list.append(base_string.lower().count("c"))
        count_list.append(base_string.lower().count("g"))
        count_list.append(base_string.lower().count("t"))

    else:
        count_list = [-1,-1,-1,-1]

    return ",".join([str(i) for i in count_list])

mpileup_input_file = sys.argv[1]
clean_vars_file = sys.argv[2]
depth_thresh = int(sys.argv[3])
std_thresh = int(sys.argv[4])
output_file = sys.argv[5]

ref_alt_dict = dict()
with open(clean_vars_file, "r") as f:
    for line in f:
        if "CHROM_POS_REF_ALT" in line:
            continue
        line_arr = line.strip().split()
        ref_alt_dict[line_arr[1]+"_"+line_arr[2]] = [line_arr[4], line_arr[5]]

chrom, pos, af, ad, ad_filter, std, std_filter, alt_base_list = [], [], [], [], [], [], [], []
with open(mpileup_input_file, "r") as f:
    for line in f:
        line_arr = line.strip().split()

        ref_alt = get_ref_alt(line_arr[0]+"_"+line_arr[1])
        ref = ref_alt[0]
        alt = ref_alt[1]
        insertion, deletion = False, False
        if len(alt) != 1:
            insertion = True
        if len(ref) != 1:
            deletion = True

        indel = insertion or deletion

        good_bases = remove_shit(line_arr[4], indel)
        af_ad = calc_af(line_arr[3], good_bases, ref, alt)
        stdev = calc_std(good_bases, line_arr[6].split(","), ref, alt)
        alt_bases = calc_alt_bases(good_bases, ref, alt)

        chrom.append(line_arr[0])
        pos.append(int(line_arr[1]))
        af.append(af_ad[0])
        ad.append(af_ad[1])
        ad_filter.append(1 if af_ad[1] > depth_thresh else 0)
        std.append(stdev)
        std_filter.append(1 if stdev > std_thresh else 0)
        alt_base_list.append(alt_bases)


mpileup_df = pd.DataFrame({"CHROM": chrom, "POS": pos, "mpileup_af": af, "mpileup_ad": ad, "mpileup_ad_pass": ad_filter, "mpileup_std": std, "mpileup_std_pass": std_filter, "alt_bases": alt_base_list})

cleaned_df = pd.read_csv(clean_vars_file, sep="\t")
cleaned_df = cleaned_df.merge(mpileup_df, on=["CHROM", "POS"])

cleaned_df.to_csv(output_file, index=False, sep="\t", na_rep="NA")
