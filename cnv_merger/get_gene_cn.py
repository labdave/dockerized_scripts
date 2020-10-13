import sys
import os
import pandas as pd

gene_bed = sys.argv[1]
seg_files = sys.argv[2:]

genes = []
with open(gene_bed, "r") as f:
	for line in f:
		line_arr = line.strip().split()
		genes.append(line_arr[3])

df = pd.DataFrame(genes, columns=['gene'])
df.drop(df.tail(1).index, inplace=True)

for file in seg_files:
	data = []
	sample = file_.replace('.called.seg.filt', '')
	with open(file, 'r') as f:
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
					if line_arr[-1] == '.':
						cr_arr = [0]
					else:
						cr_arr = [float(line_arr[-2])]
					old = gene
			else:
				i = 1
				line_arr = line.strip().split()
				gene = line_arr[3]
				if line_arr[-1] == '.':
					cr_arr = [0]
				else:
					cr_arr = [float(line_arr[-2])]
				old = gene
	df[sample] = data


df = df.reindex(sorted(df.columns, reverse=True), axis=1)
print(df)
df.to_csv('../all_gene_sample_matrix.tsv', sep='\t')