"""
Author:
Devang Thakkar

Date:
13 January 2021

Input:
This script takes in results from gatk_callcopyratiosegments__dna module in the 
CNV pipeline.

Output:
Segment calls with weighted mean amplitudes at 2 without the @ headers

Rationale: 
The default normalization used by GATK does not seem adequate. Results are more 
accurate and comparable to aCGH when the weighted mean of the segments is at 2. 
This covers not only the more common case where everything is almost normal but 
also does not affect the case of aneuploidy since that would not be detected by 
the pipeline anyway. The one case where this would be an issue is when half the 
chromosomes are amplified/deleted. The calculation only considers cardinal 
autosomes when calculating the baseline. The figures that were shown as proof of 
concept of gold standard were using this modified set of values.

Usage:

python normalize.py <gatk_seg_file> <normalized_seg_file>
"""

import math
import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

good_chroms = [
	"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
	"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
	"chr18", "chr19", "chr20", "chr21", "chr22"]

cnv_sum = 0.0
length_sum = 0.0

with open(in_file, "r") as f:
	for line in f:
		# skip headers
		if "@" in line or "MEAN_LOG2_COPY_RATIO" in line:
			continue
		line_arr = line.split()
		# skip non cardinal and sex chromosomes
		chrom = line_arr[0]
		if chrom not in good_chroms:
			continue
		length = int(line_arr[2]) - int(line_arr[1])
		length_sum += length
		cnv_sum += length*float(line_arr[4])

# weighted mean (should ideally be zero)
diff = cnv_sum/length_sum
print(diff)

file_str = ""
with open(in_file, "r") as f:
	for line in f:
		# remove header lines
		if "@" in line:
			continue
		# add title header
		if "MEAN_LOG2_COPY_RATIO" in line:
			file_str += line
		# subtract deviation
		else:
			line_arr = line.split()
			value = float(line_arr[4]) - diff
			line_arr[4] = str(round(value, 6))
			if value > math.log(1.1, 2):
				line_arr[5] = "+"
			elif value < math.log(0.9, 2):
				line_arr[5] = "-"
			else:
				line_arr[5] = "0"
			file_str += ("\t".join(line_arr)+"\n")

with open(out_file, "w") as f:
	f.write(file_str)