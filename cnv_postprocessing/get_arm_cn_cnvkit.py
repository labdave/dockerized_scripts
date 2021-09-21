"""
Author:
Devang Thakkar

Date:
21 September 2021

Input:
This script takes in results from the CNVkit segments step.

Output:
Arm calls with weighted means of individual segments in the arm.

Usage:

python get_arm_cn_cnvkit.py <arm_bed> <cnvkit_segs> <out_file> <sample_name>
"""
import sys
import os
import pandas as pd
import pybedtools

# FORMAT: CHROM POS1 POS2 ARM
arm_bed = sys.argv[1].strip()

# FORMAT: CHROM POS1 POS2 GENES LOG2 DEPTH PROBES WEIGHT CI_LO CI_HI
seg_file = sys.argv[2].strip()

# output file
output_file = sys.argv[3].strip()

# sample name
sample = sys.argv[4].strip()

arms = []
with open(arm_bed, "r") as f:
	for line in f:
		line_arr = line.replace("chr", "").strip().split()
		arms.append(line_arr[0]+line_arr[3])

df = pd.DataFrame(arms, columns=["arm"])
df = df.set_index("arm")
df[sample] = 0.0

a = pybedtools.BedTool(arm_bed)
b = pybedtools.BedTool(seg_file)
a.intersect(b, wao=True).saveas("intersected_seg_file.cns")

data = dict()
with open("intersected_seg_file.cns", "r") as f:
	for line in f:
		line_arr = line.strip().split()

		# Get the arm and seg part
		key = "\t".join(line_arr[:4]).replace("chr", "")
		value = "\t".join(line_arr[4:])

		# If a arm has more than one segment, we can access it because 
		# it is being stored as a list; which works even if there's only one
		if key in data:
			data[key].append(value)
		else:
			data[key] = [value]

for key in data:
	key_arr = key.split()
	arm = key_arr[0]+key_arr[3]
	if len(data[key]) == 1:
		# if there's only one element, our job is rather simple
		val_arr = data[key][0].split()
		if val_arr[0] == ".":
			cnv = 0.0
		else:
			cnv = float(val_arr[4])
		df.at[arm, sample] = cnv
	else:
		# We need to get the lengths of the segments and their corresponding 
		# values. bedtools intersect -loj leaves the actual positions so we 
		# need to get the actual values by comparing both fields. Next, we 
		# calculate the weighted mean using these positions
		weighted_sum = 0.0
		length_sum = 0
		for i in range(len(data[key])):
			val_arr = data[key][i].split()
			intersect = int(val_arr[-1])
			if val_arr[0] == ".":
				cnv = 0.0
			else:
				cnv = float(val_arr[4])
			weighted_sum += (intersect*cnv)
			length_sum += intersect
		weighted_cnv = weighted_sum/length_sum
		df.at[arm, sample] = weighted_cnv

print(df)
df.to_csv(output_file, sep="\t")
