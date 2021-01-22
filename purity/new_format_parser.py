import sys
import numpy as np
import pandas as pd

wl_variants = sys.argv[1]
output_file = sys.argv[2]
min_dpmax = int(sys.argv[3])

purity_arr = []
hc_count = 0
dv_count = 0
s2_count = 0

df = pd.read_csv(wl_variants, sep="\t", header=0)

for index, row in df.iterrows():
	# skip mpileup variants
	if pd.isna(row["dpMax"]):
		continue
	# add a filter for min dpMax
	if int(row["dpMax"]) < min_dpmax:
		continue
	purity_arr.append(min(2*100*row["afMax"], 100))
	if not pd.isna(row["HC_AD"]):
		hc_count += 1
	if not pd.isna(row["DV_AD"]):
		dv_count += 1
	if not pd.isna(row["S2_AD"]):
		s2_count += 1

# get mean and standard deviation of purity estimate
mean = np.mean(purity_arr)
std = np.std(purity_arr)

with open(output_file, "w") as f:
	f.write("Mean Purity\tStandard Deviation\tHC_Support\tDV_Support\tS2_Support\n")
	f.write(f"{mean}\t{std}\t{hc_count}\t{dv_count}\t{s2_count}\n")