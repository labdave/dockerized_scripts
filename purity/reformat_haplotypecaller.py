import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

new_file = ''
with open(input_file, 'r') as f:
	for line in f:
		if '#' in line:
			new_file += line
			continue
		line_arr = line.strip().split()
		if len(line_arr[3]) > 1 or len(line_arr[4]) > 1:
			continue
		if ',' in line_arr[3] or ',' in line_arr[4]:
			continue
		info = line_arr[9].split(':')
		if info[0] == '0/0' or info[0] == '1/1':
			continue
		ref = info[1].split(',')[0]
		alt = info[1].split(',')[1]
		if float(alt) + float(ref) == 0:
			continue
		af = float(alt)/(float(alt)+float(ref))
		new_file += '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(
			line_arr[0], line_arr[1], line_arr[3], line_arr[4],
			'PASS', ref, alt, af)

with open(output_file, 'w') as f:
	f.write(new_file)
