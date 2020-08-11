import sys
import argparse

parser = argparse.ArgumentParser(description='Filter rows/columns as needed')
parser.add_argument('-c', '--columns', required=True, help='bam file to get snapshots from')
parser.add_argument('-t', '--thresh', help='pe_sr threshold for sv calling')
parser.add_argument('-i', '--input', help='input vcf')
parser.add_argument('-o', '--output', help='output file')
parser.add_argument('-T', '--thresh_column', help='0-indexed column to threshold on')
args = parser.parse_args()

new_str = ''
i = 0

with open(args.input, 'r') as f:
	for line in f:
		line = line.replace(',', '\t').replace('\"', '')
		if not i:
			i = 1
			continue
		line_arr = line.split()
		if args.thresh_column is not None and args.thresh is not None:
			if int(line_arr[int(args.thresh_column)]) > int(args.thresh):
				if line_arr[1] in ['chr3', 'chr8', 'chr18'] or line_arr[3] in ['chr3', 'chr8', 'chr18']:
					# also filter for three callers
					print(line_arr[72])
					if int(line_arr[72]) == 3:
						new_str += '\t'.join(line_arr[:int(args.columns)])+'\n'
		else:
			new_str += '\t'.join(line_arr[:int(args.columns)])+'\n'

with open(args.output, 'w') as f:
	f.write(new_str)