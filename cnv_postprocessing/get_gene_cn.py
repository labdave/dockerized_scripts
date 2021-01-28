"""
Author:
Devang Thakkar

Date:
27 January 2021

Input:
This script takes in results from the focal CNV postprocessing step.

Output:
Gene calls with weighted means of individual segments in the gene.

Rationale: 
What was being used earlier was the median of individual segments in the gene. 
This may seem to be okay in general but a simple scenario in which it won't work 
is when there are multiple tiny amplifications but most of the area is neutral.

Usage:

python get_gene_cn.py <gene_bed> <gene_intersected_seg> <out_file> <sample_name>
"""
import sys
import os
import pandas as pd

# FORMAT: CHROM POS1 POS2 GENE
gene_bed = sys.argv[1].strip()

# FORMAT: CHROM POS1 POS2 GENE CHROM POS1 POS2 N_SEG MEAN_L2CR CALL
seg_file = sys.argv[2].strip()

# output file
output_file = sys.argv[3].strip()

# sample name
sample = sys.argv[4].strip()

genes = []
with open(arm_bed, "r") as f:
	for line in f:
		line_arr = line.strip().split()
		genes.append(line_arr[3])

df = pd.DataFrame(genes, columns=["Gene"])
df = df.set_index("Gene")
df[sample] = 0.0

data = dict()
with open(seg_file, "r") as f:
	for line in f:
		line_arr = line.strip().split()
		
		# skip header
		if "MEAN_LOG2_COPY_RATIO" in line:
			continue

		# Get the gene and seg part
		key = "\t".join(line_arr[:4])
		value = "\t".join(line_arr[4:])

		# If a gene has more than one segment, we can access it because 
		# it is being stored as a list; which works even if there's only one
		if key in data:
			data[key].append(value)
		else:
			data[key] = [value]

for key in data:
	key_arr = key.split()
	gene = key_arr[3]
	if len(data[key]) == 1:
		# if there's only one element, our job is rather simple
		val_arr = data[key][0].split()
		if val_arr[0] == ".":
			cnv = 0.0
		else:
			cnv = float(val_arr[4])
		df.at[gene, sample] = cnv
	else:
		# We need to get the lengths of the segments and their corresponding 
		# values. bedtools intersect -loj leaves the actual positions so we 
		# need to get the actual values by comparing both fields. Next, we 
		# calculate the weighted mean using these positions
		weighted_sum = 0.0
		length_sum = 0
		for i in range(len(data[key])):
			val_arr = data[key][i].split()
			pos1 = max(int(key_arr[1]), int(val_arr[1]))
			pos2 = min(int(key_arr[2]), int(val_arr[2]))
			if val_arr[0] == ".":
				cnv = 0.0
			else:
				cnv = float(val_arr[4])
			weighted_sum += ((pos2-pos1)*cnv)
			length_sum += (pos2-pos1)
		weighted_cnv = weighted_sum/length_sum
		df.at[arm, sample] = weighted_cnv

print(df)
df.to_csv(output_file, sep="\t")