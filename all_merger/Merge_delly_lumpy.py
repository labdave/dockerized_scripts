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

'''preparing input file for intersect bed'''
dist = 500

# order chr1 < chr2
def switch_chr_asc(line):
	arr = line.split()
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

	# check using hierarchy:
	delly_sr, delly_pe = int(delly_arr[6]), int(delly_arr[5])
	lumpy_sr, lumpy_pe = int(lumpy_arr[6]), int(lumpy_arr[5])
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

	# if delly_arr[1] == 'chr8' or delly_arr[1] == 'chr14' or lumpy_arr[1] == 'chr14' or lumpy_arr[1] == 'chr8':
	# 	print(delly_arr)
	# 	print(lumpy_arr)
	# 	print(delly_reads)
	# 	print(lumpy_reads)
		# print(merged)
	return merged			


def main():
	output_file = sys.argv[1]
	output_cons_file = output_file.replace('.vcf', '_cons.vcf')
	# Filter translocations by chr3, chr8, chr18
	chr_filter = int(sys.argv[2])
	delly_file = sys.argv[3]
	lumpy_file = sys.argv[4]

	# print(delly_file)
	# print(lumpy_file)
	# print(output_file)
	gene_list=['bcl6','myc','bcl2']
	gene_start_list = [187721377,127735434,63123346]
	gene_stop_list = [187745725,127742951,63320128]
	if chr_filter:
		chr_list = ['chr3','chr8','chr18']
	else:
		chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	chr_list_all = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	
	with open(delly_file, 'r') as d_f, open(lumpy_file, 'r') as l_f:
		list_head_delly = d_f.readline().strip().split("\t")
		list_head_lumpy = l_f.readline().strip().split("\t")
	
	str_header_delly = str.join("\t",list_head_delly[5:])
	str_header_lumpy = str.join("\t",list_head_lumpy[5:])
	
	'''Output file'''
	myfile = open(output_file, mode='wt')
	'''Output Header'''
	header = "dave_lab_id\tchr1\tpos1\tchr2\tpos2\tpe\tsr\tpe_sr\tcaller\t"
	header += str_header_delly+"\t"+str_header_lumpy
	myfile.write(header+"\n")
	lines = header.strip()+'\tCallers\tNum_callers\n'

	print('header written')
	'''Read Delly'''
	with open(delly_file, 'r') as f:
		delly_dict = dict()
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
			string = str.join("\t",p1[0:5])+"\t"+pe+"\t"+sr+"\t"+total+"\tDELLY"
			string += "\t"+str.join("\t",p1[5:])
			string += "\tNA"*len(list_head_lumpy[5:])+"\n"
			myfile.write(string)
			string, arr = switch_chr_asc(string)
			key = arr[0]+';'+arr[1]+':'+arr[2]+';'+arr[3]+':'+arr[4]
			delly_dict[key] = string
	print('delly written')

	'''Read Lumpy'''
	with open(lumpy_file, 'r') as f:
		lumpy_dict = dict()
		i = 0
		for line in f:
			if i == 0:
				i = 1
				continue
			p1 = line.strip().split("\t")
			'''columns num_discordant reads & split reads'''
			pe=p1[10]
			sr=p1[11]
			total=str(int(pe)+int(sr))
			string = str.join("\t",p1[0:5])+"\t"+pe+"\t"+sr+"\t"+total+"\tLUMPY"
			string += "\tNA"*len(list_head_delly[5:])
			string += "\t"+str.join("\t",p1[5:])+"\n"
			myfile.write(string)
			string, arr = switch_chr_asc(string)
			key = arr[0]+';'+arr[1]+':'+arr[2]+';'+arr[3]+':'+arr[4]
			lumpy_dict[key] = string
	print('lumpy written')

	myfile.close()

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

	print('dict created')
	print(len(delly_dict))
	print(len(lumpy_dict))
	for item in delly_dict:
		print(item)
		print(delly_dict[item])
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
			f.write('{0}_{1}\t{2}\t{3}\t{4}\n'.format(a[0], a[1], int(int(a[2])-dist/2), int(int(a[2])+dist/2), item))
			f.write('{0}_{1}\t{2}\t{3}\t{4}\n'.format(a[0], a[3], int(int(a[4])-dist/2), int(int(a[4])+dist/2), item))
	os.system('sort -k1,1 -k2,2n {0} -o {1}'.format(temp_file, delly_bed_file))
	with open(temp_file, 'w') as f:
		for item in lumpy_dict:
			a = item.replace(';',':').split(':')
			f.write('{0}_{1}\t{2}\t{3}\t{4}\n'.format(a[0], a[1], int(int(a[2])-dist/2), int(int(a[2])+dist/2), item))
			f.write('{0}_{1}\t{2}\t{3}\t{4}\n'.format(a[0], a[3], int(int(a[4])-dist/2), int(int(a[4])+dist/2), item))
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

	print(len(delly_2_dict))
	print(len(lumpy_2_dict))

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
	print('two callers done')
	for item in delly_remove:
		if item in delly_dict:
			delly_dict.pop(item)
	for item in lumpy_remove:
		if item in lumpy_dict:
			lumpy_dict.pop(item)
	
	# print all other single call lines
	for item in delly_dict:
		lines += delly_dict[item].strip()+'\tDELLY\t1\n'
	for item in lumpy_dict:
		lines += lumpy_dict[item].strip()+'\tLUMPY\t1\n'

	# print output
	with open(output_cons_file, 'w') as f:
		f.write(lines)


if __name__ == "__main__":
	main()
