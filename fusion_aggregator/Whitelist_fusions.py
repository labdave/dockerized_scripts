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
		wl_coords.append([line_arr[0], int(line_arr[1]), int(line_arr[2]), line_arr[3]])

wl_genes = []
with open(input_file, "r") as f:
	for line in f:
		# skip header
		if "breakpoint" in line:
			continue
		line_arr = line.strip().split()
		bp_arr = []
		bp_arr += line_arr[4].split(",")
		bp_arr += line_arr[5].split(",")
		for bp in bp_arr:
			for wl_coord in wl_coords:
				chrom = bp.split(":")[0]
				pos = int(bp.split(":")[1])
				if (chrom == wl_coord[0]) and (pos > wl_coord[1]) and (pos < wl_coord[2]):
					wl_genes.append(wl_coord[3])


wl_genes = list(set(wl_genes))

print(wl_genes)
with open(wl_fusions_file, "w") as f:
	for gene in wl_genes:
		f.write(gene+"\n")