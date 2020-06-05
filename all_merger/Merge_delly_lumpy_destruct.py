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

	'''preparing input file for intersect bed'''
	distance_num = 100

	# devang's code for merging and fixing number of callers:
	# this code assumes the presence of only three callers,
	# in this case - delly, lumpy, destruct.
	
	'''
	========== heuristic for merging calls ==========
	-> Merge calls with both breakpoints within 100bp of each other
	-> Can only merge if called from same sample and different caller
	-> The BP positions are defined by the following hierarchy
		- The one with more split reads; if equal:
		- The one with more paired reads; if equal:
		- Lumpy > Delly > Destruct
	-> Collapse rows into one for merged reads
	'''

	# order chr1 < chr2
	def switch_chr_asc(line):
		arr = line.split()
		# print(arr[1], arr[3])
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
			# print(arr[1], arr[3])
			# print('-')
		return '\t'.join(arr), arr


	# check if two lines are valid calls from two different callers
	# works only when used from inside SV_calling due to difference in analysis IDs
	def check_proximity(line1, line2, line3=None):
		dist = distance_num
		line1 = line1.replace(':', ';').split(';')
		line2 = line2.replace(':', ';').split(';')
		if line3:
			line3 = line3.replace(':', ';').split(';')
			# check sample
			if line1[0] != line2[0] or line2[0] != line3[0] or line3[0] != line1[0]:
				return False
			# check chr1
			if line1[1] != line2[1] or line2[1] != line3[1] or line3[1] != line1[1]:
				return False
			# check chr2
			if line1[3] != line2[3] or line2[3] != line3[3] or line3[3] != line1[3]:
				return False
			# check pos1
			if abs(int(line1[2])-int(line2[2])) > dist or abs(int(line2[2])-int(line3[2])) > dist or abs(int(line3[2])-int(line1[2])) > dist:
				return False
			# check pos2
			if abs(int(line1[4])-int(line2[4])) > dist or abs(int(line2[4])-int(line3[4])) > dist or abs(int(line3[4])-int(line1[4])) > dist:
				return False
		else:
			# check sample
			if line1[0] != line2[0] or line2[0] != line1[0]:
				return False
			# check chr1
			if line1[1] != line2[1] or line2[1] != line1[1]:
				return False
			# check chr2
			if line1[3] != line2[3] or line2[3] != line1[3]:
				return False
			# check pos1
			if abs(int(line1[2])-int(line2[2])) > dist or abs(int(line2[2])-int(line1[2])) > dist:
				return False
			# check pos2
			if abs(int(line1[4])-int(line2[4])) > dist or abs(int(line2[4])-int(line1[4])) > dist:
				return False
		return True

	# parse a list of lines to get merged line
	def get_merged_line(joint_val, type_):
		if type_ == 0:
			# we know passed order is [delly, destruct, lumpy]
			delly_arr = joint_val[0].split('\t')
			destruct_arr = joint_val[1].split('\t')
			lumpy_arr = joint_val[2].split('\t')

			# check using hierarchy:
			delly_sr, delly_pe = delly_arr[6], delly_arr[5]
			destruct_sr, destruct_pe = destruct_arr[6], destruct_arr[5]
			lumpy_sr, lumpy_pe = lumpy_arr[6], lumpy_arr[5]

			# split reads
			if (delly_sr > destruct_sr) and (delly_sr > lumpy_sr):
				chosen = delly_arr
			elif (destruct_sr > delly_sr) and (destruct_sr > lumpy_sr):
				chosen = destruct_arr
			elif (lumpy_sr > delly_sr) and (lumpy_sr > destruct_sr):
				chosen = lumpy_arr
			# paired reads
			elif (delly_pe > destruct_pe) and (delly_pe > lumpy_pe):
				chosen = delly_arr
			elif (destruct_pe > delly_pe) and (destruct_pe > lumpy_pe):
				chosen = destruct_arr
			elif (lumpy_pe > delly_pe) and (lumpy_pe > destruct_pe):
				chosen = lumpy_arr
			else:
				chosen = lumpy_arr

			# create merged row:
			merged = chosen[:9]
			merged.extend(delly_arr[9:31])
			merged.extend(lumpy_arr[31:54])
			merged.extend(destruct_arr[54:])

			# return merged line
			merged = '\t'.join(merged).strip()+'\tDELLY, DESTRUCT, LUMPY\t3\n'
			return merged

		if type_ == 1:
			# we know passed order is [delly, destruct]
			delly_arr = joint_val[0].split('\t')
			destruct_arr = joint_val[1].split('\t')

			# check using hierarchy:
			delly_sr, delly_pe = delly_arr[6], delly_arr[5]
			destruct_sr, destruct_pe = destruct_arr[6], destruct_arr[5]

			# split reads
			if (delly_sr > destruct_sr):
				chosen = delly_arr
			elif (destruct_sr > delly_sr):
				chosen = destruct_arr
			# paired reads
			elif (delly_pe > destruct_pe):
				chosen = delly_arr
			elif (destruct_pe > delly_pe):
				chosen = destruct_arr
			else:
				chosen = delly_arr

			# create merged row:
			merged = chosen[:9]
			merged.extend(delly_arr[9:54])
			merged.extend(destruct_arr[54:])

			# return merged line
			merged = '\t'.join(merged).strip()+'\tDELLY, DESTRUCT\t2\n'
			return merged

		if type_ == 2:
			# we know passed order is [delly, lumpy]
			delly_arr = joint_val[0].split('\t')
			lumpy_arr = joint_val[1].split('\t')

			# check using hierarchy:
			delly_sr, delly_pe = delly_arr[6], delly_arr[5]
			lumpy_sr, lumpy_pe = lumpy_arr[6], lumpy_arr[5]

			# split reads
			if (delly_sr > lumpy_sr):
				chosen = delly_arr
			elif (lumpy_sr > delly_sr):
				chosen = lumpy_arr
			# paired reads
			elif (delly_pe > lumpy_pe):
				chosen = delly_arr
			elif (lumpy_pe > delly_pe):
				chosen = lumpy_arr
			else:
				chosen = lumpy_arr

			# create merged row:
			merged = chosen[:9]
			merged.extend(delly_arr[9:31])
			merged.extend(lumpy_arr[31:])

			# return merged line
			merged = '\t'.join(merged).strip()+'\tDELLY, LUMPY\t2\n'
			return merged

		if type_ == 3:
			# we know passed order is [destruct, lumpy]
			destruct_arr = joint_val[0].split('\t')
			lumpy_arr = joint_val[1].split('\t')

			# check using hierarchy:
			destruct_sr, destruct_pe = destruct_arr[6], destruct_arr[5]
			lumpy_sr, lumpy_pe = lumpy_arr[6], lumpy_arr[5]

			# split reads
			if (destruct_sr > lumpy_sr):
				chosen = destruct_arr
			elif (lumpy_sr > destruct_sr):
				chosen = lumpy_arr
			# paired reads
			elif (destruct_pe > lumpy_pe):
				chosen = destruct_arr
			elif (lumpy_pe > destruct_pe):
				chosen = lumpy_arr
			else:
				chosen = lumpy_arr

			# create merged row:
			merged = chosen[:9]
			merged.extend(lumpy_arr[9:54])
			merged.extend(destruct_arr[54:])

			# return merged line
			merged = '\t'.join(merged).strip()+'\tDESTRUCT, LUMPY\t2\n'
			return merged
			

	# read output_file and create dict
	with open(output_file, 'r') as f:
		i = True
		lumpy_dict = dict()
		delly_dict = dict()
		destruct_dict = dict()
		count = 0
		for line in f:
			count += 1
			if count % 1000 == 0:
				print(count)
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
			if arr[8] == 'DESTRUCT':
				destruct_dict[key] = line

	print('dict created')

	''' complex procedure to get merge-able rows '''
	delly_destruct_dict, delly_lumpy_dict, destruct_lumpy_dict = dict(), dict(), dict()
	delly_destruct_lumpy_dict = dict()

	# all three callers
	delly_remove, destruct_remove, lumpy_remove = [], [], []
	for delly_item in delly_dict:
		for destruct_item in destruct_dict:
			for lumpy_item in lumpy_dict:
				if check_proximity(delly_item, destruct_item, lumpy_item):
					print('yes')
					joint_key = delly_item+'|'+destruct_item+'|'+lumpy_item
					joint_val = [delly_dict[delly_item], destruct_dict[destruct_item], lumpy_dict[lumpy_item]]
					delly_destruct_lumpy_dict[joint_key] = joint_val
					lines += get_merged_line(joint_val, type_=0)
					delly_remove.append(delly_item)
					destruct_remove.append(destruct_item)
					lumpy_remove.append(lumpy_item)
	print(list(set(delly_remove)))
	print(list(set(destruct_remove)))
	print(list(set(lumpy_remove)))
	for item in delly_remove:
		delly_dict.pop(item)
	for item in destruct_remove:
		destruct_dict.pop(item)
	for item in lumpy_remove:
		lumpy_dict.pop(item)


	print('three callers done')
	# three pairs of two callers each
	delly_remove, destruct_remove, lumpy_remove = [], [], []

	# delly destruct
	for delly_item in delly_dict:
		for destruct_item in destruct_dict:
			if check_proximity(delly_item, destruct_item):
				joint_key = delly_item+'|'+destruct_item
				joint_val = [delly_dict[delly_item], destruct_dict[destruct_item]]
				delly_destruct_dict[joint_key] = joint_val
				lines += get_merged_line(joint_val, type_=1)
				delly_remove.append(delly_item)
				destruct_remove.append(destruct_item)
				continue
	print('two callers done')
	# delly lumpy
	for delly_item in delly_dict:
		for lumpy_item in lumpy_dict:
			if check_proximity(delly_item, lumpy_item):
				joint_key = delly_item+'|'+lumpy_item
				joint_val = [delly_dict[delly_item], lumpy_dict[lumpy_item]]
				delly_lumpy_dict[joint_key] = joint_val
				lines += get_merged_line(joint_val, type_=2)
				delly_remove.append(delly_item)
				lumpy_remove.append(lumpy_item)
				continue
	print('two callers done')
	# destruct lumpy
	for destruct_item in destruct_dict:
		for lumpy_item in lumpy_dict:
			if check_proximity(destruct_item, lumpy_item):
				joint_key = destruct_item+'|'+lumpy_item
				joint_val = [destruct_dict[destruct_item], lumpy_dict[lumpy_item]]
				destruct_lumpy_dict[joint_key] = joint_val
				lines += get_merged_line(joint_val, type_=3)
				destruct_remove.append(destruct_item)
				lumpy_remove.append(lumpy_item)
				continue
	print('two callers done')
	print(list(set(delly_remove)))
	print(list(set(destruct_remove)))
	print(list(set(lumpy_remove)))
	for item in delly_remove:
		delly_dict.pop(item)
	for item in destruct_remove:
		destruct_dict.pop(item)
	for item in lumpy_remove:
		lumpy_dict.pop(item)

	# print all other single call lines
	for item in delly_dict:
		# print(delly_dict[item])
		lines += delly_dict[item].strip()+'\tDELLY\t1\n'
	for item in lumpy_dict:
		print(item)
		lines += lumpy_dict[item].strip()+'\tLUMPY\t1\n'
	for item in destruct_dict:
		# print(item)
		lines += destruct_dict[item].strip()+'\tDESTRUCT\t1\n'

	# print output
	with open(output_cons_file, 'w') as f:
		f.write(lines)
	"""
	# read output file and check for matches
	with open(output_file, 'r') as f:
		i = True
		for line in f:
			if i:
				output_cons_text = line.strip()+'\tNum_Callers\tCallers\n'
				i = False
				continue
			arr = line.strip().split()
			callers = [arr[8]]
			key = arr[0]+';'+arr[1]+':'+arr[2]+';'+arr[3]+':'+arr[4]+';'+arr[8]
			for item in output_dict:
				if check_proximity(item, key):
					callers.append(item.split(';')[-1])
			callers = list(set(callers))
			n_callers = len(callers)
			callers_str = ','.join(callers)
			output_cons_text += (line.strip()+'\t{0}\t{1}\n'.format(n_callers, callers_str))




	# naresh's code for number of callers:
	fobj = open(output_file)
	tmp_bed = output_file+'.tmp.bed'
	tmp1_bed = output_file+'.tmp1.bed'
	tmp1_cp_bed = output_file+'.tmp1.cp.bed'
	tmp_all_bed = output_file+'.tmp.all.bed'
	myfile = open(tmp_bed, mode='wt')
	header = fobj.readline()
	linenum = 0
	for i in fobj:
		linenum = linenum+1
		i = i.strip()
		arr = i.split("\t")
		'''merging the sample id with chromosome to do intersect bed at sample level'''
		myfile.write(arr[0].split("_")[0]+"__"+arr[1]+"\t"+str(int(arr[2])-distance_num)+"\t"+str(int(arr[2])+distance_num)+"\t"+str(linenum)+"\t"+str(1)+"\t"+arr[8]+"\n")
		# print(arr[0]+"__"+arr[1]+"\t"+str(int(arr[2])-distance_num)+"\t"+str(int(arr[2])+distance_num)+"\t"+str(linenum)+"\t"+str(1)+"\t"+arr[8])
		myfile.write(arr[0].split("_")[0]+"__"+arr[3]+"\t"+str(int(arr[4])-distance_num)+"\t"+str(int(arr[4])+distance_num)+"\t"+str(linenum)+"\t"+str(2)+"\t"+arr[8]+"\n")
		# print(arr[0]+"__"+arr[3]+"\t"+str(int(arr[4])-distance_num)+"\t"+str(int(arr[4])+distance_num)+"\t"+str(linenum)+"\t"+str(2)+"\t"+arr[8])
	fobj.close()
	myfile.close()
	
	# debug
	# with open(tmp_bed, 'r') as f:
	# 	for line in f:
	# 		print(line)

	'''Intersect bed'''
	cmd = 'sort -k2,2n -k3,3n  ' + tmp_bed + ' > ' + tmp1_bed
	os.system(cmd)
	
	cmd = 'cp ' + tmp1_bed + ' ' + tmp1_cp_bed
	os.system(cmd)
	
	cmd = 'bedtools intersect -a ' + tmp1_bed + ' -b ' + tmp1_cp_bed + ' -wao > '+ tmp_all_bed
	os.system(cmd)
	
	'''Reading the intersect bed output file'''
	merged_bed={}
	fobj = open(tmp_all_bed)
	for i in fobj:
		i = i.strip()
		# print(i)
		arr = i.split("\t")
		'''separating sample & chr of first & second chr'''
		lst1 = arr[0].split("__")
		lst2 = arr[6].split("__")
		'''Overlap > 0, for same sample, and ignoring same caller'''
		if int(arr[12]) > 0 and lst1[0].split("_")[0]==lst2[0].split("_")[0] and arr[5]!=arr[11]:
			'''Creating key sample_caller_linenumber_BP(1 or 2)'''
			key = lst1[0]+' '+arr[5]+' '+arr[3]+' '+arr[4]
			'''Other caller information'''
			val = arr[11]
			# print(key)
			if key in merged_bed:
				merged_bed[key]=merged_bed[key]+','+val
			else:
				merged_bed[key]=val	
	fobj.close()
	
	'''Creating new output file for other caller information'''
	header = [line for line in open(output_file, 'r')][0].strip()

	myfile = open(output_cons_file, mode='wt')
	myfile.write(header + "\tCallers\tNum_Callers\n")
	linenum = 0
	with open(output_file, 'r') as f:
		for i in f:
			if linenum == 0:
				linenum += 1
				continue
			i = i.strip()
			# print(i)
			arr = i.split("\t")
			linenum = linenum + 1
			'''checking if both BPs in the dict'''
			key1 = arr[0].split('_')[0]+' '+arr[8]+' '+str(linenum)+' '+"1"
			key2 = arr[0].split('_')[0]+' '+arr[8]+' '+str(linenum)+' '+"2"
			val = "NA\tNA"
			# print(key1)
			# print(key2)
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
	"""
if __name__ == "__main__":
	main()
