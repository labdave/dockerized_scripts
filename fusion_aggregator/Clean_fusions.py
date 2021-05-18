# the point of this script is to merge fusions from the output file from arriba
# we want nearby fusions to be merged together. We can decide thereshold but 
# will use 10000 bp right now

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
merge_thresh = float(sys.argv[3])
sample_id = sys.argv[4]

i = 0
reverse = False
positions = []

with open(input_file, "r") as f:
	for line in f:
		
		# skip header
		if not i:
			i = 1
			continue
		
		# extract positions
		line_arr = line.strip().split("\t")
		chrom1 = line_arr[4].split(":")[0]
		pos1 = int(line_arr[4].split(":")[1])
		chrom2 = line_arr[5].split(":")[0]
		pos2 = int(line_arr[5].split(":")[1])

		positions.append({"chrom1": chrom1, "pos1": pos1, "chrom2": chrom2, "pos2": pos2})

delete = []
merge = []

for i in range(len(positions)):
	for j in range(i, len(positions)):
		pair1 = positions[i]
		pair2 = positions[j]
		if i == j:
			continue
		# skip lines that were already merged
		if i in delete or j in delete:
			continue
		if ((pair1["chrom1"] == pair2["chrom1"]) and 
			(abs(pair1["pos1"]-pair2["pos1"]) < merge_thresh) and 
			(pair1["chrom2"] == pair2["chrom2"]) and
			(abs(pair1["pos2"]-pair2["pos2"]) < merge_thresh)):
			# merge stuff
			merge.append([i+1, j+1])  # account for header
			delete.append(j)

print(merge)
# read file again, merge lines, write to file_str
lines = []
with open(input_file, "r") as f:
	lines = f.readlines()
	for item in merge:
		line1 = lines[item[0]]
		line2 = lines[item[1]]
		line_arr1 = line1.strip().split()
		line_arr2 = line2.strip().split()
		for i in range(len(line_arr1)):
			if line_arr2[i] == line_arr1[i] or line_arr2[i] in line_arr1[i].split(","):
				continue
			line_arr1[i] = line_arr1[i]+","+line_arr2[i]
		lines[item[0]] = "\t".join(line_arr1)
		lines[item[1]] = ""

empty = True
with open(output_file, "w") as f:
	for line in lines:
		# add to header
		if "#gene1" in line:
			line = "#sample\t"+line.lstrip("#")
			f.write(line.strip()+"\n")
		# filter out low confidence calls
		if "high" in line or "medium" in line:
			f.write(sample_id+"\t"+line.strip()+"\n")
			empty = False

	if empty:
		none_line = "\t".join(["N/A" for i in line.strip().split("\t")])
		f.write(sample_id+"\t"+none_line)