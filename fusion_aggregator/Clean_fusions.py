# the point of this script is to remove low confidence calls and add sample id

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
sample_id = sys.argv[3]

empty = ""
with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
    for line in f_in:
    	if "#gene1" in line:
			line = "#sample\t"+line.lstrip("#")
			f_out.write(line.strip()+"\n")
		# filter out low confidence calls
		if "high" in line or "medium" in line:
			f_out.write(sample_id+"\t"+line.strip()+"\n")
			empty = False

	if empty:
		none_line = "\t".join(["N/A" for i in line.strip().split("\t")])
		f_out.write(sample_id+"\t"+none_line)
