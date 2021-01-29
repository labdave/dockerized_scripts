# this script merges single sample one-column tables into a matrix

import numpy as np
import pandas as pd
import sys

input_files = sys.argv[1].split("-")
output_file = sys.argv[2]

template_file = input_files.pop(0)
with open(template_file, "r") as f:
	for line in f:
		header = line.strip().split()
		mode = header[0]
		break

print(mode)

df = pd.read_csv(input_files[0], sep="\t", index_col=mode)

big_df = df
for i in range(len(input_files)):
	if i:
		input_file = input_files[i]
		df1 = pd.read_csv(input_file, sep="\t", index_col=mode)
		big_df = pd.concat([df, df1], axis=1)

print(big_df)

big_df.to_csv(output_file, sep="\t")
