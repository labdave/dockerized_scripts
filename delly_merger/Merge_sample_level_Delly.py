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
	input_files = []
	output_file = sys.argv[1]
	temp_file = output_file+'.tmp'
	davelab_ids = sys.argv[2].split('?')
	chr_switch = int(sys.argv[3])
	# Filter translocations by chr3, chr8, chr18
	chr_filter = int(sys.argv[4])
	for i in range(len(sys.argv)-5):
		input_files.append(sys.argv[i+5])

	gene_list=['bcl6','myc','bcl2']
	gene_start_list = [187721377,127735434,63123346]
	gene_stop_list = [187745725,127742951,63320128]
	if chr_filter:
		chr_list = ['chr3','chr8','chr18']
	else:
		chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	chr_list_all = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
		
	'''Temp output file'''
	myfile = open(temp_file, mode='wt')
	
	'''Header'''
	myfile.write("dave_lab_id\t")
	myfile.write("Delly_CHR1\tDelly_POS1\tDelly_CHR2\tDelly_POS2\tDelly_REF\tDelly_ALT\tDelly_Distance1\tDelly_Distance2\tDelly_Filter\tDelly_PRECISE_status\t")
	myfile.write("Delly_PE_NReads\tDelly_PE_MAPQ\tDelly_SR_NReads\tDelly_SR_MAPQ\tDelly_GT\tDelly_GL\tDelly_GQ\tDelly_FT\tDelly_RCL\tDelly_RC\tDelly_RCR\tDelly_CN\tDelly_DR\tDelly_DV\tDelly_RR\tDelly_RV\n")
	
	i = -1
	for file_ in input_files:
		i += 1

		dict_clean = {}
		with open(file_, 'r') as f:
			for line in f:
				line=line.strip()
				''' Selecting rows with Translocations'''
				if "SVTYPE=BND" in line:
					'''Extracting fields'''
					p = line.split("\t")
					chr1 = p[0]
					pos1 = int(p[1])
					id = p[2]
					ref = p[3]
					alt = p[4]
					qual = p[5]
					filter = p[6]	
					info = p[7]    
					format = p[8]
					sample = p[9]
					gt_list = sample.replace(":","\t")
					precise = '0'
					if info.startswith("PRECISE"):
						precise = '1'
					chr2 = re.search(';CHR2=(.+?);', info).group(1)
					pos2 = int(re.search(';END=(.+?);', info).group(1))
					pe = re.search(';PE=(.+?);', info).group(1)
					try:
						sr = re.search(';SR=(.+?);', info).group(1)
						srmq = re.search(';SRMAPQ=(.+?);', info).group(1)
					except:
						sr = "0"
						srmq = "0"
					pemq = re.search(';MAPQ=(.+?);', info).group(1)
					'''Ignoring translocation events between same primary and alt chrs'''
					chr1_list = chr1.split("_")
					chr2_list = chr2.split("_")
					if chr1 != chr2 and ( chr1 in chr_list or chr2 in chr_list) and ( chr1_list[0] != chr2_list[0]) and (chr1_list[0] in chr_list_all and chr2_list[0] in chr_list_all):
						'''Ignore the Distance, it is going to be recalculated in the later steps'''
						dist1 = -1
						dist2 = -1
						if chr_filter == 1:
							if chr1 in chr_list:
								indx1 = chr_list.index(chr1)
								indx_pos1 = pos1
								if indx_pos1 < gene_start_list[indx1]:
									dist1 = indx_pos1-gene_start_list[indx1]
								elif indx_pos1 > gene_stop_list[indx1]:
									dist1 = indx_pos1-gene_stop_list[indx1]
								else:
									dist1 = 0
							if chr2 in chr_list:
								indx2 = chr_list.index(chr2)
								indx_pos2 = pos2
								if indx_pos2 < gene_start_list[indx2]:
									dist2 = indx_pos2-gene_start_list[indx2]
								elif indx_pos2 > gene_stop_list[indx2]:
									dist2 = indx_pos2-gene_stop_list[indx2]
								else:
									dist2 = 0
		
						myfile.write(davelab_ids[i]+"\t")
						'''Switch breakpoints if required, First BP MYC, BCL2, BCL6'''
						if chr_switch == 1:
							if chr1 in chr_list :
								myfile.write(chr1+"\t"+str(pos1)+"\t"+chr2+"\t"+str(pos2))
							else:
								myfile.write(chr2+"\t"+str(pos2)+"\t"+chr1+"\t"+str(pos1))
						else:
							myfile.write(chr1+"\t"+str(pos1)+"\t"+chr2+"\t"+str(pos2))
						myfile.write("\t"+ref+"\t"+alt+"\t"+str(dist1)+"\t"+str(dist2)+"\t"+filter+"\t"+precise+"\t")
						myfile.write(pe+"\t"+pemq+"\t"+sr+"\t"+srmq+"\t"+gt_list+"\n")

	myfile.close()

	''' Remove duplicates '''
	pos1 = -1
	pos2 = -1
	i = 0
	myfile = open(output_file, mode='wt')
	with open(temp_file, 'r') as f:
		for line in f:
			# print header line
			if i == 0:
				i += 1
				myfile.write(line)
				continue
			# print first line
			if i == 1:
				pos1 = line_arr[2]
				pos2 = line_arr[4]
				myfile.write(line)
				continue
			# every other line
			line_arr = line.strip().split()
			if line_arr[2] != pos2 or line_arr[4] != pos1:
				myfile.write(line)
			pos1 = line_arr[2]
			pos2 = line_arr[4]
	myfile.close()

	
if __name__ == "__main__":
	main()