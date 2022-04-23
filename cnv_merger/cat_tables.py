# this script cats together emultiple tables with same columns

import sys

from os import listdir
from os.path import isfile, join

# declare folder where files are to be found
folder = "/data"
# obtain list of files
onlyfiles = [f for f in listdir(folder) if isfile(join(folder, f))]
# get file type
input_file_type = sys.argv[1].strip()
# filter for file type
input_files = [f for f in onlyfiles if (f.endswith(f".{input_file_type}.seg") and "arm" not in f and "gene" not in f and "cyto" not in f)]

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