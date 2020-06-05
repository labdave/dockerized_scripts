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
	destruct_file = sys.argv[5]

	print(delly_file)
	print(lumpy_file)
	print(destruct_file)
	print(output_file)
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
	list_head_destruct = [line for line in open(destruct_file, 'r')][0].strip().split("\t")
	
	str_header_delly = str.join("\t",list_head_delly[5:])
	str_header_lumpy = str.join("\t",list_head_lumpy[5:])
	str_header_destruct = str.join("\t",list_head_destruct[5:])
	
	'''Output file'''
	myfile = open(output_file, mode='wt')
	'''Output Header'''
	myfile.write("dave_lab_id\tchr1\tpos1\tchr2\tpos2\tpe\tsr\tpe_sr\tcaller\t")
	# print("dave_lab_id\tchr1\tpos1\tchr2\tpos2\tpe\tsr\tpe_sr\tcaller\t"+str_header_delly+"\t"+str_header_lumpy+"\t"+str_header_destruct)
	myfile.write(str_header_delly+"\t"+str_header_lumpy+"\t"+str_header_destruct+"\n")

	print('header written')
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
			# print(str.join("\t",p1[0:5])+"\t"+pe+"\t"+sr+"\t"+total+"\tDELLY", end='')
			myfile.write("\t"+str.join("\t",p1[5:]))
			# print("\t"+str.join("\t",p1[5:]), end='')
			myfile.write("\tNA"*len(list_head_lumpy[5:]))
			# print("\tNA"*len(list_head_lumpy[5:]), end='')
			myfile.write("\tNA"*len(list_head_destruct[5:])+"\n")
			# print("\tNA"*len(list_head_destruct[5:])+"\n", end='')
	print('delly written')
	'''Read Lumpy'''
	with open(lumpy_file, 'r') as f:
		i = 0
		for line in f:
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
			# print(str.join("\t",p1[0:5])+"\t"+pe+"\t"+sr+"\t"+total+"\tLUMPY", end='')
			myfile.write("\tNA"*len(list_head_delly[5:]))
			# print("\tNA"*len(list_head_delly[5:]), end='')
			myfile.write("\t"+str.join("\t",p1[5:]))
			# print("\t"+str.join("\t",p1[5:]), end='')
			myfile.write("\tNA"*len(list_head_destruct[5:])+"\n")
			# print("\tNA"*len(list_head_destruct[5:])+"\n", end='')
	print('lumpy written')
	'''Read Destruct'''
	with open(destruct_file, 'r') as f:
		i = 0
		for line in f:
			# skip header line
			if i == 0:
				i = 1
				continue
			p1 = line.strip().split("\t")
			pe=p1[16]
			sr=p1[9]
			total=str(int(pe)+int(sr))
			myfile.write(str.join("\t",p1[0:5])+"\t"+pe+"\t"+sr+"\t"+total+"\tDESTRUCT")
			# print(str.join("\t",p1[0:5])+"\t"+pe+"\t"+sr+"\t"+total+"\tDESTRUCT", end='')
			myfile.write("\tNA"*len(list_head_delly[5:]))
			# print("\tNA"*len(list_head_delly[5:]), end='')
			myfile.write("\tNA"*len(list_head_lumpy[5:]))
			# print("\tNA"*len(list_head_lumpy[5:]), end='')
			myfile.write("\t"+str.join("\t",p1[5:])+"\n")
			# print("\t"+str.join("\t",p1[5:])+"\n", end='')
	print('Destruct written')
	myfile.close()

	
if __name__ == "__main__":
	main()
