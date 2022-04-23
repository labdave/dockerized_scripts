# this script merges single sample one-column tables into a matrix

import numpy as np
import pandas as pd
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
input_files = [f for f in onlyfiles if f.endswith(f".{input_file_type}.seg")]

output_file = sys.argv[2]

print(input_files)
# skip the first
template_file = input_files.pop(0)
with open(template_file, "r") as f:
	for line in f:
		header = line.strip().split()
		mode = header[0]
		break

print(mode)

# get the first
df = pd.read_csv(template_file, sep="\t", index_col=mode)

# start the big table with just the first
big_df = df
# add more if there are
for input_file in input_files:
	print(input_file)
	df1 = pd.read_csv(input_file, sep="\t", index_col=mode)
	big_df = pd.concat([big_df, df1], axis=1)

print(big_df)

big_df.to_csv(output_file, sep="\t")
