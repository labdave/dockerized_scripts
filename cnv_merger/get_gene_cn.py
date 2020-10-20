import sys
import os
import pandas as pd

gene_bed = sys.argv[1].strip()
seg_files_file = sys.argv[2].strip()
sample_names_file = sys.argv[3].strip()
output_file = sys.argv[4].strip()

genes = []
with open(gene_bed, "r") as f:
	for line in f:
		line_arr = line.strip().split()
		genes.append(line_arr[3])

df = pd.DataFrame(genes, columns=["gene"])
df.drop(df.tail(1).index, inplace=True)

seg_files, sample_names = [], []
with open(seg_files_file, "r") as f:
	seg_files = [line.strip() for line in f]
with open(sample_names_file, "r") as f:
	sample_names = [line.strip() for line in f]

for i in range(len(seg_files)):
	file = seg_files[i]
	sample = sample_names[i]
	data = []
	with open(file, "r") as f:
		i = 0
		for line in f:
			# I think what is happening is to get the max amplitude of copy number for a gene
			if i:
				line_arr = line.strip().split()
				gene = line_arr[3]
				if gene == old:
					cr_arr.append(float(line_arr[-2]))
				else:
					cr = max(cr_arr, key=abs)
					data.append(str(cr))
					if line_arr[-1] == ".":
						cr_arr = [0]
					else:
						cr_arr = [float(line_arr[-2])]
					old = gene
			else:
				i = 1
				line_arr = line.strip().split()
				gene = line_arr[3]
				if line_arr[-1] == ".":
					cr_arr = [0]
				else:
					cr_arr = [float(line_arr[-2])]
				old = gene
	df[sample] = data


df = df.reindex(sorted(df.columns, reverse=True), axis=1)
print(df)
df.to_csv(output_file, sep="\t", index=False)