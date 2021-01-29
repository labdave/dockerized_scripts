# this script merges single sample one-column tables into a matrix

import numpy as np
import pandas as pd
import sys

input_files = sys.argv[1].split("-")
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
