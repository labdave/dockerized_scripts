# the point of this script is to take an intersected bed file and find 
# read pairs that intersect with two different exons of the same gene

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

read_map_dict = dict()
with open(input_file, "r") as f:
	for line in f:
		line_arr = line.strip().split()
		read_name = line_arr[3].split("/")[0]
		gene_info = "_".join(line_arr[12:16])
		if read_name in read_map_dict:
			read_map_dict[read_name].append(gene_info)
		else:
			read_map_dict[read_name] = [gene_info]

with open(output_file, "w") as f:
	for read_name in read_map_dict:

		# skip read pairs which only have one read mapped completely to an exon
		if len(read_map_dict[read_name]) == 1:
			continue

		# skip read pairs which only have both read pairs mapped to the same exon
		if len(list(set(read_map_dict[read_name]))) == 1:
			continue

		# skip read pairs which have both read pairs mapped to different genes
		if read_map_dict[read_name][0].split("_")[-1] != read_map_dict[read_name][1].split("_")[-1]:
			continue

		f.write(read_name+"\n")