import sys

whitelist_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

whitelist = []

with open(whitelist_file, 'r') as f:
	for line in f:
		if 'Start' in line:
			continue
		line_arr = line.strip().split()
		whitelist.append(line_arr[1]+':'+line_arr[2])

new_file = ''
with open(input_file, 'r') as f:
	for line in f:
		line_arr = line.strip().split()
		if line_arr[0]+':'+line_arr[1] in whitelist:
			new_file += line

with open(output_file, 'w') as f:
	f.write(new_file)