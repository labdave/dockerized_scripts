# this script cats together emultiple tables with same columns

import sys

input_files = sys.argv[1].split("-")
output_file = sys.argv[2]

header, new_lines = "", ""

for input_file in input_files:
	with open(input_file, "r") as f:
		i = 0
		for line in f:
			if not header:
				header = line
				i += 1
				continue
			if not i:
				i += 1
				continue
			new_lines += line

with open(output_file, "w") as f:
	f.write(header+new_lines)