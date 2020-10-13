import sys
import os
import pandas as pd

cyto_bed = sys.argv[1].strip()
seg_files = sys.argv[2].strip().split(",")
sample_names = sys.argv[3].strip().split(",")
output_file = sys.argv[4].strip()

bands = []
with open(cyto_bed, "r") as f:
	for line in f:
		line_arr = line.replace("chr", "").split()
		bands.append(line_arr[0]+line_arr[3])

df = pd.DataFrame(bands, columns=["cytoBand"])
df.drop(df.tail(1).index, inplace=True)

for i in range(len(seg_files)):
	file = seg_files[i]
	sample = sample_names[i]
	data = []
	with open(file, "r") as f:
		i = 0
		for line in f:
			# I think what is happening is to get the max amplitude of copy number for a cytoband
			if i:
				line_arr = line.strip().split()
				pos = line_arr[0].replace("chr", "")+line_arr[3]
				if pos == old:
					cr_arr.append(float(line_arr[-2]))
				else:
					cr = max(cr_arr, key=abs)
					data.append(str(cr))
					if line_arr[-1] == ".":
						cr_arr = [0]
					else:
						cr_arr = [float(line_arr[-2])]
					old = pos
			else:
				i = 1
				line_arr = line.strip().split()
				pos = line_arr[0].replace("chr", "")+line_arr[3]
				if line_arr[-1] == ".":
					cr_arr = [0]
				else:
					cr_arr = [float(line_arr[-2])]
				old = pos
	df[sample] = data


df = df.reindex(sorted(df.columns, reverse=True), axis=1)
print(df)
df.to_csv(output_file, sep="\t")