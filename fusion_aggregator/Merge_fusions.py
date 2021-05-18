# this script takes in the cleaned/wl fusion files and creates 
# an aggregated matrix with fusions stacked in rows

import sys

input_files = sys.argv[1].split("-")
output_file = sys.argv[2]

# get header
with open(input_files[0], "r") as f:
	for line in f:
		if line[0] == "#":
			header = line.strip()

# concat files without header
file_str = ""
for file_ in input_files:
	with open(file_, "r") as f:
		for line in f:
			if line[0] == "#":
				continue
			file_str += line.strip("\n")+"\n"

# write to output file
with open(output_file, "w") as f:
	f.write(file_str)