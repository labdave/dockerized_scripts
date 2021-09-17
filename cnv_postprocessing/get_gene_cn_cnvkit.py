# the point of this script is to parse the cnr file from cnvkit
# into a gene: gene_cn table

# the cnv value for a gene is the weighted average of the log2 ratios
# for each of the baits where the weights are the length of the bait
# multiplied by the weight provided by cnvkit

import pandas as pd
import sys

cnr = sys.argv[1]
gene_seg = sys.argv[2]
sample_id = sys.argv[3]

df = pd.read_csv(cnr, sep="\t", header=0, names=["chromosome", "start", "end", "gene", "depth", "log2", "weight"])

for index, row in df.iterrows():
	df.at[index, "length"] = (float(df.at[index, "end"]) - float(df.at[index, "start"]))*float(df.at[index, "weight"])
	df.at[index, "value"] = df.at[index, "length"]*df.at[index, "log2"]

df = df.drop(columns=["chromosome", "start", "end", "depth", "log2", "weight"])
df = df.groupby("gene", as_index=False, sort=False).sum()
df["cnv"] = df["value"]/df["length"]
df = df.drop(columns=["length", "value"])
df = df.rename(columns={"gene": "Gene", "cnv": sample_id})

df.to_csv(gene_seg, sep="\t", index=False, float_format='%.4f')