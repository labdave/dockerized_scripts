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
	print(len(sys.argv))
	print(str(sys.argv))
	input_files = []
	output_file = sys.argv[1]
	davelab_ids = sys.argv[2].split('?')
	chr_switch = sys.argv[3]
	# Filter translocations by chr3, chr8, chr18
	chr_filter = sys.argv[4]
	for i in range(len(sys.argv)-5):
		input_files.append(sys.argv[i+5])

	print(output_file)
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
	myfile = open(output_file, mode='wt')
	
	'''Header'''
	myfile.write("dave_lab_id\t")
	myfile.write("Lumpy_CHROM1\tLumpy_POS1\tLumpy_CHROM2\tLumpy_POS2\tLumpy_PRECISE_status\tLumpy_QUAL\tLumpy_FILTER\t")
	myfile.write("Lumpy_GT\tLumpy_SU\tLumpy_PE\tLumpy_SR\tLumpy_GQ\tLumpy_SQ\tLumpy_GL\tLumpy_DP\tLumpy_RO\tLumpy_AO\tLumpy_QR\tLumpy_QA\tLumpy_RS\tLumpy_AS\tLumpy_ASC\tLumpy_RP\tLumpy_AP\tLumpy_AB\tLumpy_Distance1\tLumpy_Distance2\n")
	
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
					gt_list = sample.split(":")
					info_list = info.split(";")
					su_info = info_list[len(info_list)-3].replace("SU=","")
					pe_info = info_list[len(info_list)-2].replace("PE=","")
					sr_info = info_list[len(info_list)-1].replace("SR=","")
					gt_list[1]=su_info
					gt_list[2]=pe_info
					gt_list[3]=sr_info
					gt = str.join("\t",gt_list)
					precise='1'
					if "IMPRECISE" in info:
						precise='0'
					alt = alt.replace("[",'|')
					alt = alt.replace("]",'|')
					alt1 = alt.split('|')
					alt_list = alt1[1].split(':')
					chr2 = alt_list[0]
					pos2 = int(alt_list[1])
					'''Ignoring translocation events between same primary and alt chrs'''
					chr1_list = chr1.split("_")
					chr2_list = chr2.split("_")
					if (chr1 != chr2) and ( chr1 in chr_list or chr2 in chr_list) and ( chr1_list[0] != chr2_list[0]) and (chr1_list[0] in chr_list_all and chr2_list[0] in chr_list_all):
						'''Ignore the Distance, it is going to be recalculated in the later steps'''
						dist1=-1
						dist2=-1
						if chr1 in chr_list:
							indx1 = chr_list.index(chr1)
							indx_pos1 = pos1
							if indx_pos1 < gene_start_list[indx1]:
								dist1=indx_pos1-gene_start_list[indx1]
							elif indx_pos1 > gene_stop_list[indx1]:
								dist1=indx_pos1-gene_stop_list[indx1]
							else:
								dist1=0
						if chr2 in chr_list:
							indx2 = chr_list.index(chr2)
							indx_pos2 = pos2
							if indx_pos2 < gene_start_list[indx2]:
								dist2=indx_pos2-gene_start_list[indx2]
							elif indx_pos2 > gene_stop_list[indx2]:
								dist2=indx_pos2-gene_stop_list[indx2]
							else:
								dist2=0	
						'''Switch breakpoints if required, First BP MYC, BCL2, BCL6'''
						if chr_switch:
							if chr1 in chr_list :
								ln = chr1+"\t"+str(pos1)+"\t"+chr2+"\t"+str(pos2)
							else:
								ln = chr2+"\t"+str(pos2)+"\t"+chr1+"\t"+str(pos1)
						else:
							ln = chr1+"\t"+str(pos1)+"\t"+chr2+"\t"+str(pos2)
						if not ln in dict_clean:
							dict_clean[ln]=1	
							myfile.write(davelab_ids[i]+"\t"+ln+"\t")
							myfile.write(precise+"\t"+qual+"\t"+filter+"\t"+gt+"\t"+str(dist1)+"\t"+str(dist2)+"\n")
	myfile.close()
	
	
if __name__ == "__main__":
	main()