# the point of this script is to go through the list of cleaned fusions
# and check if a whitelist fusion was found

import sys

input_file = sys.argv[1]
whitelist = sys.argv[2]
wl_fusions_file = sys.argv[3]

wl_coords = []
with open(whitelist, "r") as f:
	for line in f:
		line_arr = line.strip().split()
		wl_coords.append([line_arr[0], line_arr[1], line_arr[2], line_arr[3]])

wl_genes = []
with open(input_file, "r") as f:
	for line in f:
		line_arr = line.strip().split()
		bp_arr = []
		bp_arr += line_arr[4].split(",")
		bp_arr += line_arr[5].split(",")
		for bp in bp_arr:
			for wl_coord in wl_coords:
				chrom = bp.split(":")[0]
				pos = bp.split(":")[1]
				if (chrom == wl_coords[0]) and (pos > wl_coords[1]) and (pos < wl_coords[2]):
					wl_genes.append(wl_coords[4])


wl_genes = list(set(wl_genes))

with open(wl_fusions_file, "w") as f:
	for gene in wl_genes:
		f.write(gene+"\n")