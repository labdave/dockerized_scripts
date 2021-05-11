# this script adds a flag column that specifies if the gene was found in
# a panel of normals blacklist

import sys
import pandas as pd

input_file = sys.argv[1]
gene_blacklist_file = sys.argv[2]

data = pd.read_csv(input_file, sep="\t")

gene_blacklist = []
with open(gene_blacklist_file, "r") as f:
	for line in f:
		gene_blacklist.append(line.strip())

data["blacklist"] = 0

for index, row in data.iterrows():
	if row["Gene"] in gene_blacklist:
		data["blacklist"][index] = 1

cols = data.columns.tolist()
cols = cols[:1] + cols[-1:] + cols[1:-1]
data = data[cols]

data.to_csv(input_file, sep="\t", index=False)