import sys
import numpy as np

dv = sys.argv[1]
s2 = sys.argv[2]
hc = sys.argv[3]
output_file = sys.argv[4]

purity_arr = []
with open(dv, 'r') as f:
	for line in f:
		purity_arr.append(100*float(line.strip().split()[-1]))
with open(s2, 'r') as f:
	for line in f:
		purity_arr.append(100*float(line.strip().split()[-1]))
with open(hc, 'r') as f:
	for line in f:
		purity_arr.append(100*float(line.strip().split()[-1]))

mean = np.mean(purity_arr)
std = np.std(purity_arr)
with open(output_file, 'w') as f:
	f.write('Mean Purity\t Standard deviation\n{0}\t{1}'.format(mean, std))
