# the point of this script is to remove low confidence calls and add sample id

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
sample_id = sys.argv[3]

lines = []
with open(input_file, "r") as f:
    for line in f:
    	if "#gene1" in line:
			line = "#sample\t"+line.lstrip("#")
			lines.append(line)
    		print("Printing header")
    		print(line)
		# filter out low confidence calls
		if "high" in line or "medium" in line:
			line = sample_id+"\t"+line.strip()+"\n"
			lines.append(line)
    		print("Printing line")
    		print(line)

	if len(lines) == 1:
		print("Taking care of empty file")
		none_line = "\t".join(["N/A" for i in line.strip().split("\t")])
		f_out.write(sample_id+"\t"+none_line)

with open(output_file, "w") as f:
	print("Writing file")
	for line in lines:
		f.write(line)