import sys
import os
import pandas as pd

rootdir = 'intersect'
newdir = 'clean'

bands = []
with open('cytoBand.txt', 'r') as f:
	for line in f:
		line_arr = line.replace('chr', '').split()
		bands.append(line_arr[0]+line_arr[3])

df = pd.DataFrame(bands, columns=['cytoBand'])
df.drop(df.tail(1).index,inplace=True)

for subdir, dirs, files in os.walk(rootdir):
	for file_ in files:
		file = os.path.join(subdir, file_)
		data = []
		sample = file_.replace('.called.seg', '')
		with open(file, 'r') as f:
			i = 0
			for line in f:
				if i:
					line_arr = line.strip().split()
					pos = line_arr[0].replace('chr', '')+line_arr[3]
					if pos == old:
						cr_arr.append(float(line_arr[-2]))
					else:
						cr = max(cr_arr, key=abs)
						data.append(str(cr))
						if line_arr[-1] == '.':
							cr_arr = [0]
						else:
							cr_arr = [float(line_arr[-2])]
						old = pos
				else:
					i += 1
					line_arr = line.strip().split()
					pos = line_arr[0].replace('chr', '')+line_arr[3]
					if line_arr[-1] == '.':
						cr_arr = [0]
					else:
						cr_arr = [float(line_arr[-2])]
					old = pos
		df[sample] = data


df = df.reindex(sorted(df.columns, reverse=True), axis=1)
print(df)
df.to_csv('cytoband_sample_matrix.tsv', sep='\t')