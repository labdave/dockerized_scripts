import sys
import os
import pandas as pd

cyto_bed = sys.argv[1].strip()
seg_files_file = sys.argv[2].strip()
sample_names_file = sys.argv[3].strip()
output_file = sys.argv[4].strip()

bands = []
with open(cyto_bed, "r") as f:
	for line in f:
		line_arr = line.replace("chr", "").split()
		bands.append(line_arr[0]+line_arr[3])

df = pd.DataFrame(bands, columns=["cytoBand"])
df = df.set_index("cytoBand")

seg_files, sample_names = [], []
with open(seg_files_file, "r") as f:
	seg_files = [line.strip() for line in f]
with open(sample_names_file, "r") as f:
	sample_names = [line.strip() for line in f]

for i in range(len(seg_files)):
	file = seg_files[i]
	sample = sample_names[i]
	data = []
	# I think what used to happen got the max amplitude of copy number for a cytoband
	# now using median instead - can use pandas
	file_df = pd.read_csv(file, sep="\t", names=["chrom", "pos1", "pos2", "cyto", "info", "chrom_", "pos1_", "pos2_", "n_seg", "cnv", "sign"])
	file_df["chrom_cyto"] = file_df["chrom"] + file_df["cyto"]
	file_df["chrom_cyto"] = file_df["chrom_cyto"].str[3:]
	file_df = file_df.groupby(["chrom_cyto"]).median().drop(["pos1", "pos2", "pos1_", "pos2_"], axis=1)
	file_df.rename(columns={"cnv":sample}, inplace=True)

	df = df.join(file_df)

df = df.reindex(sorted(df.columns, reverse=True), axis=1)
print(df)
df.to_csv(output_file, sep="\t")