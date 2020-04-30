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


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


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

	print(input_files)

	gene_list=['bcl6','myc','bcl2']
	gene_start_list = [187721377,127735434,63123346]
	gene_stop_list = [187745725,127742951,63320128]
	if chr_filter:
		chr_list = ['chr3','chr8','chr18']
	else:
		chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	chr_list_all = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	
	'''Output  file'''
	myfile = open(temp_file, mode='wt')
	
	'''Header'''
	myfile.write("dave_lab_id\t")
	myfile.write("Destruct_chromosome_1\tDestruct_position_1\tDestruct_chromosome_2\tDestruct_position_2\tDestruct_strand_1\tDestruct_strand_2\tDestruct_prediction_id\tDestruct_homology\tDestruct_num_split\tDestruct_mate_score\tDestruct_template_length_1\tDestruct_log_cdf\tDestruct_template_length_2\tDestruct_log_likelihood\tDestruct_template_length_min\tDestruct_num_reads\tDestruct_num_unique_reads\n")
	
	i = -1
	for file_ in input_files:
		i += 1
		dict_clean = {}
		with open(file_, 'r') as f:
			for line in f:
				line=line.strip()
				# print(line)
				p=line.split('\t')
				type_ = p[18]
				''' Selecting rows with Translocations'''
				if type_ == 'translocation':
					'''Extracting fields'''
					prediction_id = p[0]   
					chr1 = p[1]
					chr1 = chr1.replace('.','v')
					if is_number(chr1) or chr1 == 'X' or chr1 == 'Y':
						chr1 = 'chr'+chr1
					strand_1 = p[2]
					pos1 = int(p[3])      
					chr2 = p[4]
					chr2 = chr2.replace('.','v')
					if is_number(chr2) or chr2 == 'X' or chr2 == 'Y':
						chr2 = 'chr'+chr2
					strand_2 = p[5]        
					pos2 = int(p[6])      
					homology = p[7]        
					num_split = p[8]       
					inserted = p[9]        
					mate_score = p[10]      
					template_length_1 = p[11]       
					log_cdf = p[12] 
					template_length_2  = p[13]      
					log_likelihood  = p[14] 
					template_length_min = p[15]     
					num_reads = p[16]       
					num_unique_reads = p[17]
					
					'''Ignoring translocation events between same primary and alt chrs'''		
					chr1_list = chr1.split("_")
					chr2_list = chr2.split("_")
					if (chr1 != chr2) and (chr1 in chr_list or chr2 in chr_list) and (chr1_list[0] != chr2_list[0]) and ((chr1_list[0] in chr_list_all) and (chr2_list[0] in chr_list_all)) and (int(num_reads)>=2) and (int(num_split)>=2) and (float(log_likelihood) > -20.0) and (int(template_length_min)>70):
						'''Switch breakpoints if required, First BP MYC, BCL2, BCL6'''
						if chr_switch == 1:
							if chr1 in chr_list :
								ln = chr1+"\t"+str(pos1)+"\t"+chr2+"\t"+str(pos2)
							else:
								ln = chr2+"\t"+str(pos2)+"\t"+chr1+"\t"+str(pos1)
						else:
							ln = chr1+"\t"+str(pos1)+"\t"+chr2+"\t"+str(pos2)
						if not ln in dict_clean:
							dict_clean[ln]=1
							myfile.write(davelab_ids[i]+"\t")
							myfile.write(ln+"\t"+strand_1+"\t"+strand_2+"\t"+prediction_id+"\t"+homology+"\t"+num_split+"\t"+mate_score+"\t"+template_length_1+"\t"+log_cdf+"\t"+template_length_2)
							myfile.write("\t"+log_likelihood+"\t"+template_length_min+"\t"+num_reads+"\t"+num_unique_reads+"\n")
	myfile.close()

	''' Remove duplicates '''
	pos1 = -1
	pos2 = -1
	i = 0
	myfile = open(output_file, mode='wt')
	with open(temp_file, 'r') as f:
		for line in f:
			line_arr = line.strip().split()
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
			if line_arr[2] != pos2 or line_arr[4] != pos1:
				myfile.write(line)
			pos1 = line_arr[2]
			pos2 = line_arr[4]
	myfile.close()


if __name__ == "__main__":
	main()