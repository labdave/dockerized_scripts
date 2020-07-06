import sys

anno_file = sys.argv[1]
thresh = int(sys.argv[2])
out_file = sys.argv[3]

new_str = ''

with open(anno_file, 'r') as f:
	for line in f:
		line_arr = line.split()
		if int(line_arr[7]) > thresh:
			new_str += '\t'.join(line_arr[:5])+'\n'

with open(out_file, 'w') as f:
	f.write(new_str)