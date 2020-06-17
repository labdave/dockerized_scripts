import sys
import numpy as np

dv = sys.argv[1]
s2 = sys.argv[2]
hc = sys.argv[3]
output_file = sys.argv[4]

purity_arr = []
dv_count = 0
s2_count = 0
hc_count = 0

with open(dv, 'r') as f:
	for line in f:
		est = 2*100*float(line.strip().split()[-1])
		purity_arr.append(min(100, est))
		dv_count += 1
with open(s2, 'r') as f:
	for line in f:
		est = 2*100*float(line.strip().split()[-1])
		purity_arr.append(min(100, est))
		s2_count += 1
with open(hc, 'r') as f:
	for line in f:
		est = 2*100*float(line.strip().split()[-1])
		purity_arr.append(min(100, est))
		hc_count += 1

mean = np.mean(purity_arr)
std = np.std(purity_arr)
with open(output_file, 'w') as f:
	f.write('Mean Purity\t Standard deviation\tDeepvariant support\tStrelka2 support\tHaplotypecaller support\n')
	f.write('{0}\t{1}\t{2}\t{3}\t{4}\n\n\n\n\n'.format(mean, std, dv_count, s2_count, hc_count))

	f.write('\n\nDeepvariant calls')
	with open(dv, 'r') as g:
		for line in g:
			f.write(line)
	f.write('\n\nStrelka2 calls')
	with open(s2, 'r') as g:
		for line in g:
			f.write(line)
	f.write('\n\nHaplotypecaller calls')
	with open(hc, 'r') as g:
		for line in g:
			f.write(line)