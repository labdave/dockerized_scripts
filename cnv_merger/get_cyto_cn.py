"""
Author:
Devang Thakkar

Date:
13 January 2021

Input:
This script takes in results from individual CNV pipeline results.

Output:
Cytoband calls with weighted means of individual segments in the cytoband.

Rationale: 
What was being used earlier was the median of idividual segments in the cytoband. 
This may seem to be okay in general but a simple scenario in which it won't work 
is when there are multiple tiny amplifications but most of the area is neutral.

Usage:

python get_cyto_cn.py <cyto_file> <seg_files_file> <sample_names_file> <out_file>
"""
import sys
import os
import pandas as pd

# FORMAT: CHROM POS1 POS2 CYTOBAND INFO
cyto_bed = sys.argv[1].strip()

# FORMAT: INTERSECTED_SEG_FILE_PATH
seg_files_file = sys.argv[2].strip()

# FORMAT: SAMPLE
sample_names_file = sys.argv[3].strip()

# output file
output_file = sys.argv[4].strip()

bands = []
with open(cyto_bed, "r") as f:
	for line in f:
		line_arr = line.replace("chr", "").split()
		bands.append(line_arr[0]+line_arr[3])

seg_files, sample_names = [], []
with open(seg_files_file, "r") as f:
	seg_files = [line.strip() for line in f]
with open(sample_names_file, "r") as f:
	sample_names = [line.strip() for line in f]

df = pd.DataFrame(bands, columns=["cytoBand"])
df = df.set_index("cytoBand")
for sample in sample_names:
	df[sample] = 0.0

for i in range(len(seg_files)):
	file = seg_files[i]
	sample = sample_names[i]
	data = dict()
	# FORMAT: CHROM POS1 POS2 CYTOBAND INFO CHROM POS1 POS2 N_SEG MEAN_L2CR CALL
	with open(file, "r") as f:
		for line in f:
			line_arr = line.strip().split()
			
			# Get the cytoband and seg part
			key = "\t".join(line_arr[:5])
			value = "\t".join(line_arr[5:])

			# If a cytoband has more than one segment, we can access it because 
			# it is being stored as a list; which works even if there's only one
			if key in data:
				data[key].append(value)
			else:
				data[key] = [value]
	for key in data:
		key_arr = key.split()
		cyto = key_arr[0].replace("chr", "")+key_arr[3]
		if len(data[key]) == 1:
			# if there's only one element, our job is rather simple
			val_arr = data[key][0].split()
			if val_arr[0] == ".":
				cnv = 0.0
			else:
				cnv = float(val_arr[4])
			df.at[cyto, sample] = cnv
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
			df.at[cyto, sample] = weighted_cnv

print(df)
df.to_csv(output_file, sep="\t")