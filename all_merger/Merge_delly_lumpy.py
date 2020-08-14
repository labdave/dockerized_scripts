# editing Naresh's scripts for compatilibility with CC -Devang

#!/usr/bin/python3
from __future__ import print_function
__author__ = "Naresh Prodduturi"
__email__ = "naresh.prodduturi@ddb.bio"
__status__ = "Dev"

import os
import argparse
import sys
import pwd
import time
import subprocess
import re
import shutil	
import logging



def main():
	output_file = sys.argv[1]
	output_cons_file = output_file.replace('.vcf', '_cons.vcf')
	# Filter translocations by chr3, chr8, chr18
	chr_filter = int(sys.argv[2])
	delly_file = sys.argv[3]
	lumpy_file = sys.argv[4]

	print(delly_file, file=sys.stderr)
	print(lumpy_file, file=sys.stderr)
	print(output_file, file=sys.stderr)
	gene_list=['bcl6','myc','bcl2']
	gene_start_list = [187721377,127735434,63123346]
	gene_stop_list = [187745725,127742951,63320128]
	if chr_filter:
		chr_list = ['chr3','chr8','chr18']
	else:
		chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	chr_list_all = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	
	list_head_delly = [line for line in open(delly_file, 'r')][0].strip().split("\t")
	list_head_lumpy = [line for line in open(lumpy_file, 'r')][0].strip().split("\t")
	
	str_header_delly = str.join("\t",list_head_delly[5:])
	str_header_lumpy = str.join("\t",list_head_lumpy[5:])
	
	'''Output file'''
	myfile = open(output_file, mode='wt')
	'''Output Header'''
	myfile.write("dave_lab_id\tchr1\tpos1\tchr2\tpos2\tpe\tsr\tpe_sr\tcaller\t")
	print("dave_lab_id\tchr1\tpos1\tchr2\tpos2\tpe\tsr\tpe_sr\tcaller\t"+str_header_delly+"\t"+str_header_lumpy, file=sys.stderr)
	myfile.write(str_header_delly+"\t"+str_header_lumpy+"\t"+"\n")

	print('header written', file=sys.stderr)
	'''Read Delly'''
	with open(delly_file, 'r') as f:
		i = 0
		for line in f:
			# skip header line
			if i == 0:
				i = 1
				continue
			p1 = line.strip().split("\t")
			'''columns num_discordant reads & split reads'''
			pe=p1[11]
			sr=p1[13]
			total=str(int(pe)+int(sr))
			myfile.write(str.join("\t",p1[0:5])+"\t"+pe+"\t"+sr+"\t"+total+"\tDELLY")
			# print(str.join("\t",p1[0:5])+"\t"+pe+"\t"+sr+"\t"+total+"\tDELLY", end='', file=sys.stderr)
			myfile.write("\t"+str.join("\t",p1[5:]))
			# print("\t"+str.join("\t",p1[5:]), end='', file=sys.stderr)
			myfile.write("\tNA"*len(list_head_lumpy[5:])+"\n")
			# print("\tNA"*len(list_head_lumpy[5:]), end='', file=sys.stderr)
	print('delly written', file=sys.stderr)
	'''Read Lumpy'''
	with open(lumpy_file, 'r') as f:
		i = 0
		for line in f:
			print(line, file=sys.stderr)
			# skip header line
			if i == 0:
				i = 1
				continue
			p1 = line.strip().split("\t")
			'''columns num_discordant reads & split reads'''
			pe=p1[10]
			sr=p1[11]
			total=str(int(pe)+int(sr))
			myfile.write(str.join("\t",p1[0:5])+"\t"+pe+"\t"+sr+"\t"+total+"\tLUMPY")
			# print(str.join("\t",p1[0:5])+"\t"+pe+"\t"+sr+"\t"+total+"\tLUMPY", end='', file=sys.stderr)
			myfile.write("\tNA"*len(list_head_delly[5:]))
			# print("\tNA"*len(list_head_delly[5:]), end='', file=sys.stderr)
			myfile.write("\t"+str.join("\t",p1[5:])+"\n")
			# print("\t"+str.join("\t",p1[5:]), end='', file=sys.stderr)
	print('lumpy written', file=sys.stderr)

	'''preparing input file for intersect bed'''
	dist = 500

	# devang's code for merging and fixing number of callers:
	# this code assumes the presence of only two callers,
	# in this case - delly, lumpy.
	
	'''
	========== heuristic for merging calls ==========
	-> Merge calls with both breakpoints within 100bp of each other
	-> Can only merge if called from same sample and different caller
	-> The BP positions are defined by the following hierarchy
		- The one with more reads; if equal:
		- Lumpy > Delly 
	-> Collapse rows into one for merged reads
	'''

	# order chr1 < chr2
	def switch_chr_asc(line):
		arr = line.split()
		# print(arr[1], arr[3], file=sys.stderr)
		if arr[1] == 'chrX' or arr[1] == 'chrY':
			temp_chr = arr[1]
			temp_pos = arr[2]
			arr[1] = arr[3]
			arr[2] = arr[4]
			arr[3] = temp_chr
			arr[4] = temp_pos
			return '\t'.join(arr), arr
		if arr[3] == 'chrX' or arr[3] == 'chrY':
			return '\t'.join(arr), arr
		chr1 = int(arr[1].replace('chr', ''))
		chr2 = int(arr[3].replace('chr', ''))
		if chr1 > chr2:
			temp_chr = arr[1]
			temp_pos = arr[2]
			arr[1] = arr[3]
			arr[2] = arr[4]
			arr[3] = temp_chr
			arr[4] = temp_pos
			# print(arr[1], arr[3], file=sys.stderr)
			# print('-', file=sys.stderr)
		return '\t'.join(arr), arr


	# check if two lines are valid calls from two different callers
	# works only when used from inside SV_calling due to difference in analysis IDs
	def check_proximity(line1, line2):
		line1 = line1.replace(':', ';').split(';')
		line2 = line2.replace(':', ';').split(';')
		if all([
				line1[1] == line2[1],
				line1[3] == line2[3],
				abs(int(line1[2])-int(line2[2])) < dist,
				abs(int(line1[4])-int(line2[4])) < dist,
				]):
			return True
		return False


	# parse a list of lines to get merged line
	def get_merged_line(joint_val):
		# we know passed order is [delly, lumpy]
		delly_arr = joint_val[0].split('\t')
		lumpy_arr = joint_val[1].split('\t')
		# print(delly_arr)
		# print(lumpy_arr)

		# check using hierarchy:
		delly_sr, delly_pe = delly_arr[6], delly_arr[5]
		lumpy_sr, lumpy_pe = lumpy_arr[6], lumpy_arr[5]
		delly_reads = delly_sr+delly_pe
		lumpy_reads = lumpy_sr+lumpy_pe

		# split+paired reads
		if (lumpy_reads > delly_reads):
			chosen = lumpy_arr
		elif (delly_reads > lumpy_reads):
			chosen = delly_arr
		else:
			chosen = lumpy_arr

		# create merged row:
		merged = chosen[:9]
		merged.extend(delly_arr[9:31])
		merged.extend(lumpy_arr[31:])

		# return merged line
		merged = '\t'.join(merged).strip()+'\tDELLY, LUMPY\t2\n'
		return merged			


	# read output_file and create dict
	with open(output_file, 'r') as f:
		i = True
		lumpy_dict = dict()
		delly_dict = dict()
		for line in f:
			if i:
				lines = line.strip()+'\tCallers\tNum_callers\n'
				i = False
				continue
			line, arr = switch_chr_asc(line)
			key = arr[0]+';'+arr[1]+':'+arr[2]+';'+arr[3]+':'+arr[4]
			if arr[8] == 'DELLY':
				delly_dict[key] = line
			if arr[8] == 'LUMPY':
				lumpy_dict[key] = line

	print('dict created', file=sys.stderr)
	print(len(delly_dict), file=sys.stderr)
	print(len(lumpy_dict), file=sys.stderr)
	for item in delly_dict:
		print(item, file=sys.stderr)
		print(delly_dict[item], file=sys.stderr)
		break

	''' bed file method to get mergeable rows '''
	temp_file = output_file.replace('.vcf', '.tmp')
	delly_bed_file = temp_file.replace('.tmp', '.delly.bed')
	lumpy_bed_file = temp_file.replace('.tmp', '.lumpy.bed')
	
	delly_lumpy_bed = temp_file.replace('.tmp', '.delly.lumpy.bed')
	lumpy_delly_bed = temp_file.replace('.tmp', '.lumpy.delly.bed')
	
	# get bed files
	with open(temp_file, 'w') as f:
		for item in delly_dict:
			a = item.replace(';',':').split(':')
			f.write('{0}_{1}\t{2}\t{3}\t{4}\n'.format(a[0], a[1], int(a[2])-dist/2, int(a[2])+dist/2, item))
			f.write('{0}_{1}\t{2}\t{3}\t{4}\n'.format(a[0], a[3], int(a[4])-dist/2, int(a[4])+dist/2, item))
	os.system('sort -k1,1 -k2,2n {0} -o {1}'.format(temp_file, delly_bed_file))
	with open(temp_file, 'w') as f:
		for item in lumpy_dict:
			a = item.replace(';',':').split(':')
			f.write('{0}_{1}\t{2}\t{3}\t{4}\n'.format(a[0], a[1], int(a[2])-dist/2, int(a[2])+dist/2, item))
			f.write('{0}_{1}\t{2}\t{3}\t{4}\n'.format(a[0], a[3], int(a[4])-dist/2, int(a[4])+dist/2, item))
	os.system('sort -k1,1 -k2,2n {0} -o {1}'.format(temp_file, lumpy_bed_file))

	# find intersections
	os.system('bedtools intersect -u -a {0} -b {1} > {2}'.format(delly_bed_file, lumpy_bed_file, delly_lumpy_bed))
	os.system('bedtools intersect -u -a {0} -b {1} > {2}'.format(lumpy_bed_file, delly_bed_file, lumpy_delly_bed))

	# get list of caller specific 2-caller ids
	delly_2_list, lumpy_2_list = [], []
	delly_2_dict, lumpy_2_dict = dict(), dict()

	with open(delly_lumpy_bed, 'r') as f:
		for line in f:
			delly_2_list.append(line.strip().split()[-1])
	delly_2_list = list(set(delly_2_list))
	for item in delly_2_list:
		delly_2_dict[item] = delly_dict[item]

	with open(lumpy_delly_bed, 'r') as f:
		for line in f:
			lumpy_2_list.append(line.strip().split()[-1])
	lumpy_2_list = list(set(lumpy_2_list))
	for item in lumpy_2_list:
		lumpy_2_dict[item] = lumpy_dict[item]

	print(len(delly_2_dict), file=sys.stderr)
	print(len(lumpy_2_dict), file=sys.stderr)

	''' complex procedure to get merge-able rows '''
	delly_lumpy_dict = dict()

	# two callers
	delly_remove, lumpy_remove = [], []

	# delly lumpy
	for delly_item in delly_2_dict:
		for lumpy_item in lumpy_2_dict:
			if check_proximity(delly_item, lumpy_item) and delly_item not in delly_remove and lumpy_item not in lumpy_remove:
				joint_key = delly_item+'|'+lumpy_item
				# print(joint_key)
				joint_val = [delly_dict[delly_item], lumpy_dict[lumpy_item]]
				# print(joint_val)
				delly_lumpy_dict[joint_key] = joint_val
				lines += get_merged_line(joint_val)
				delly_remove.append(delly_item)
				lumpy_remove.append(lumpy_item)
				continue
	print('two callers done', file=sys.stderr)
	for item in delly_remove:
		if item in delly_dict:
			delly_dict.pop(item)
	for item in lumpy_remove:
		if item in lumpy_dict:
			lumpy_dict.pop(item)
	# print(lines, file=sys.stderr)
	
	# print all other single call lines
	for item in delly_dict:
		# print(delly_dict[item], file=sys.stderr)
		lines += delly_dict[item].strip()+'\tDELLY\t1\n'
	for item in lumpy_dict:
		# print(item, file=sys.stderr)
		lines += lumpy_dict[item].strip()+'\tLUMPY\t1\n'

	# print output
	with open(output_cons_file, 'w') as f:
		f.write(lines)


if __name__ == "__main__":
	main()
