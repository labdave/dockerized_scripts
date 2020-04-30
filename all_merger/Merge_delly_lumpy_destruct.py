# editing Naresh's scripts for compatilibility with CC -Devang

#!/usr/bin/python3
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
	davelab_ids = sys.argv[2].split('?')
	# Filter translocations by chr3, chr8, chr18
	chr_filter = int(sys.argv[3])
	delly_file = sys.argv[4]
	lumpy_file = sys.argv[5]
	destruct_file = sys.argv[6]

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
	
	list_head_delly = [line for line in open(delly_file, 'r')][0].split("\t")
	list_head_lumpy = [line for line in open(lumpy_file, 'r')][0].split("\t")
	list_head_destruct = [line for line in open(destruct_file, 'r')][0].split("\t")
	
	str_header_delly = str.join("\t",list_head_delly[5:])
	str_header_lumpy = str.join("\t",list_head_lumpy[5:])
	str_header_destruct = str.join("\t",list_head_destruct[5:])
	
	'''Output file'''
	myfile = open(output_file, mode='wt')
	'''Output Header'''
	myfile.write("dave_lab_id\tchr1\tpos1\tchr2\tpos2\tpe\tsr\tpe_sr\tcaller\t")
	myfile.write(str_header_delly+"\t"+str_header_lumpy+"\t"+str_header_destruct+"\n")

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
			myfile.write("\t"+str.join("\t",p1[5:]))
			myfile.write("\tNA"*len(list_head_lumpy[5:]))
			myfile.write("\tNA"*len(list_head_destruct[5:])+"\n")
	
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
			myfile.write("\tNA"*len(list_head_delly[5:]))
			myfile.write("\t"+str.join("\t",p1[5:]))
			myfile.write("\tNA"*len(list_head_destruct[5:])+"\n")
	
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
			myfile.write("\tNA"*len(list_head_delly[5:]))
			myfile.write("\tNA"*len(list_head_lumpy[5:]))
			myfile.write("\t"+str.join("\t",p1[5:])+"\n")
	
	myfile.close()

	'''preparing input file for intersect bed'''
	distance_num = 100
	fobj = open(output_file)
	tmp_bed = output_file+'.tmp.bed'
	tmp1_bed = output_file+'.tmp1.bed'
	tmp1_cp_bed = output_file+'.tmp1.cp.bed'
	tmp_all_bed = output_file+'tmp.all.bed'
	myfile = open(tmp_bed, mode='wt')
	header = fobj.readline()
	linenum = 0
	for i in fobj:
		linenum = linenum+1
		i = i.strip()
		arr = i.split("\t")
		'''merging the sample id with chromosome to do intersect bed at sample level'''
		myfile.write(arr[0]+"__"+arr[1]+"\t"+str(int(arr[2])-distance_num)+"\t"+str(int(arr[2])+distance_num)+"\t"+str(linenum)+"\t"+str(1)+"\t"+arr[8]+"\n")
		myfile.write(arr[0]+"__"+arr[3]+"\t"+str(int(arr[4])-distance_num)+"\t"+str(int(arr[4])+distance_num)+"\t"+str(linenum)+"\t"+str(2)+"\t"+arr[8]+"\n")
	fobj.close()
	myfile.close()
	
	'''Intersect bed'''
	cmd = 'sort -T ' + dir_path + '  -k2,2n -k3,3n  ' + tmp_bed + ' > ' + tmp1_bed
	os.system(cmd)
	
	cmd = 'cp ' + tmp1_bed + ' ' + tmp1_cp_bed
	os.system(cmd)
	
	cmd = 'intersectBed -a ' + tmp1_bed + ' -b ' + tmp1_cp_bed + ' -wao > '+ tmp_all_bed
	os.system(cmd)
	
	'''Reading the intersect bed output file'''
	merged_bed={}
	fobj = open(tmp_all_bed)
	for i in fobj:
		i = i.strip()
		arr = i.split("\t")
		'''separating sample & chr of first & second chr'''
		lst1 = arr[0].split("__")
		lst2 = arr[6].split("__")
		'''Overlap > 0, for same sample, and ignoring same caller'''
		if int(arr[12]) > 0 and lst1[0]==lst2[0] and arr[5]!=arr[11]:
			'''Creating key sample_caller_linenumber_BP(1 or 2)'''
			key = lst1[0]+' '+arr[5]+' '+arr[3]+' '+arr[4]
			'''Other caller information'''
			val = arr[11]
			if key in merged_bed :
				merged_bed[key]=merged_bed[key]+','+val
			else:
				merged_bed[key]=val	
	fobj.close()
	
	'''Creating new output file for other caller information'''
	header = [line for line in open(delly_file, 'r')][0].strip()

	myfile = open(output_cons_file, mode='wt')
	myfile.write(header + "\tNum_Callers\tCallers\n")
	linenum = 0
	with open(delly_file, 'r') as f:
		for i in f:
			i = i.strip()
			arr = i.split("\t")
			linenum = linenum + 1
			'''checking if both BPs in the dict'''
			key1 = arr[0]+' '+arr[12]+' '+str(linenum)+' '+"1"
			key2 = arr[0]+' '+arr[12]+' '+str(linenum)+' '+"2"
			val = "NA\tNA"
			if key1 in merged_bed and key2 in merged_bed:
				tmp1 = arr[12] + ',' + merged_bed[key1]
				tmp1_list = tmp1.split(',')
				tmp2 = arr[12] + ',' + merged_bed[key2]
				tmp2_list = tmp2.split(',')
				tmp_list = list(set(tmp1_list) & set(tmp2_list))
				tmp = str.join(',', tmp_list)
				if len(tmp_list) > 1:
					val = tmp + "\t" + str(len(tmp_list))
			myfile.write(i + "\t" + val + "\n")	
	myfile.close()
	
if __name__ == "__main__":
	main()