# the point of this script is to create easy to read FARDEEP results

import os
import sys

import pandas as pd

input_file = sys.argv[1]

data = pd.read_csv(input_file, sep="\t", header=[0,1]).rename(columns={"Unnamed: 0_level_0": "Sample"})

# print(data.transpose())
data = data.T

mapping = {"B.immature..pro-b": "Pro B cell", "B.immature..pre-b": "Pre B cell", "B.immature..pro.b": "Pro B cell", "B.immature..pre.b": "Pre B cell", "B.mature..naive": "Naive B cell", "B.mature..germinal...dz": "Dark Zone GC B cell", "B.mature..germinal...lz": "Light Zone GC B cell", "B.mature..memory": "Memory B cell", "B.mature..plasma": "Plasma B cell", "T.CD8": "CD8 T cell", "T.CD4..naive": "Naive CD4 T cell", "T.CD4..memory...resting": "Resting Memory CD4 T cell", "T.CD4..memory...activated": "Activated Memory CD4 T cell", "T.follicular_helper": "Follicular Helper T cell", "T.other..regulatory": "Regulatory T cell", "T.other..gamma_delta": "Gamma Delta T cell", "T.NK..resting": "Resting NK T cell", "T.NK..activated": "Activated NK T cell", "Myeloid.macrophages..M0": "M0 Macrophage", "Myeloid.macrophages..M1": "M1 Macrophage", "Myeloid.macrophages..M2": "M2 Macrophage", "Myeloid.dendritic..resting": "Resting Dendritic cell", "Myeloid.dendritic..activated": "Activated Dendritic cell", "Myeloid.other..mast...resting": "Resting Mast cell", "Myeloid.other..mast...activated": "Activated Mast cell", "Myeloid.other..eosinophils": "Eosinophil", "Myeloid.other..neutrophils": "Neutrophil", "Myeloid.other..monocytes": "Monocyte", "Myeloid.langerhans": "Langerhans cell", "B": "B", "T":"T", "Myeloid": "Myeloid", "T.NK": "NK T Cell", "Myeloid.other..mast": "Mast cell", "T.CD4..memory": "Memory CD4 T cell", "Myeloid.dendritic": "Dendritic cell", "Myeloid.plasmacytoid": "Plasmacytoid Dendritic Cell", "T.plasmacytoid": "Plasmacytoid Dendritic Cell"}

lines = ["Sample\tDiagnosis\t#1 Cell type\t#1 Cell type percent\t#2 Cell type\t#2 Cell type percent\t#3 Cell type\t#3 Cell type percent\n"]
for index, row in data.iterrows():
	if index == ('Sample', 'Diagnosis'):
		cell_types = list(row)
		continue
	row = list(row)
	new_line = index[0]+"\t"+index[1]+"\t"

	# get #1
	first_max = max(row)
	if first_max > 0.1:
		cell_type = mapping[cell_types[row.index(first_max)]]
		new_line += cell_type+"\t"+str(int(first_max*100))+"\t"
	else:
		new_line += "N/A\tN/A\t"
	row = [-1 if x == first_max else x for x in row]
	
	# get #2
	second_max = max(row)
	if second_max > 0.1:
		cell_type = mapping[cell_types[row.index(second_max)]]
		new_line += cell_type+"\t"+str(int(second_max*100))+"\t"
	else:
		new_line += "N/A\tN/A\t"
	row = [-1 if x == second_max else x for x in row]

	# get #3
	third_max = max(row)
	if third_max > 0.1:
		cell_type = mapping[cell_types[row.index(third_max)]]
		new_line += cell_type+"\t"+str(int(third_max*100))+"\n"
	else:
		new_line += "N/A\tN/A\n"

	lines.append(new_line)
	new_line = ""

with open(sys.argv[2], "w") as f:
	for line in lines:
		f.write(line)

new_df = pd.read_csv(sys.argv[2], sep="\t")
new_df = new_df.sort_values(["Diagnosis", "#1 Cell type", "#2 Cell type", "#3 Cell type", "#1 Cell type percent"], ascending=[True, True, True, True, False])
new_df.to_csv(sys.argv[2], sep="\t", index=False)