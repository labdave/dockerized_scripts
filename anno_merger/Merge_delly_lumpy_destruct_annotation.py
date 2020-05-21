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
import numpy as np 

def input_file_validity(file):
	'''function to check if input files exists and valid''' 	
	if os.path.exists(file)==False:
		raise argparse.ArgumentTypeError( '\nERROR:Path:\n'+file+':Does not exist')
	if os.path.isfile(file)==False:
		raise argparse.ArgumentTypeError( '\nERROR:File expected:\n'+file+':is not a file')
	if os.access(file,os.R_OK)==False:
		raise argparse.ArgumentTypeError( '\nERROR:File:\n'+file+':no read access ')
	return file

	
def create_dir(f):
	'''Checks if a directory is present or not and creates one if not present'''
	if not os.path.exists(f):
		os.makedirs(f)

def argument_parse():
	'''Parses the command line arguments'''#required="True",
	parser=argparse.ArgumentParser(description='')
	parser.add_argument("-i","--input",help="Input merged file",type=input_file_validity)
	parser.add_argument("-c","--capture_kit",help="Capture kit file")
	parser.add_argument("-o","--output_file",help="Output file")
	parser.add_argument("-t","--temp_dir",help="Temp Directory")
	parser.add_argument("-g","--dac_gap",help="DAC, GAP, Delly blacklist")
	parser.add_argument("-r","--rep_mas",help="Repeat masker,seg dup blacklist")
	parser.add_argument("-O","--other_sample",help="Other Sample",type=bool)
	parser.add_argument("-n","--normal_samp",help="normal Samples")
	parser.add_argument("-l","--level1_bp",help="Level1 BP file")
	parser.add_argument("-G","--gtf_file",help="GTF File")
	parser.add_argument("-C","--chr_names",help="CHR Names")
	parser.add_argument("-v","--var_files",help="Variant files directory")
	parser.add_argument("-p","--paper_freq_pairs",help="Frequent gene pairs in the paper")
	parser.add_argument("-e","--gene_expr",help="Gene expr results")
	parser.add_argument("-f","--fusion_out_dir",help="Fusion catcher output per sample")
	parser.add_argument("-D","--detected_trans",help="Rachel detected_trans file")
	parser.add_argument("-F","--chr_filter",help="chr filter for myc,bcl2,bcl6")
	return parser

def shell_cmd(cmd):
	cmd = re.sub('\s+', ' ', cmd).strip()
	cmd = cmd.split(" ")
	m = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	stdout, stderr = m.communicate()
	exitCode = m.returncode
	return exitCode

def intersectBed_run(file1,file2,out,temp):
	'''Intersect bed'''
	cmd = 'sort -k2,2n -k3,3n  ' + file1 + ' > temp.bed'
	os.system(cmd)
	
	cmd = 'mv temp.bed ' + file1
	os.system(cmd)
	
	cmd = 'sort -k2,2n -k3,3n  ' + file2 + ' > temp.bed'
	os.system(cmd)
	
	cmd = 'mv temp.bed ' + file2
	os.system(cmd)
	
	cmd = 'bedtools intersect -a ' + file1 + ' -b ' + file2 + ' -wao >' + out
	os.system(cmd)
	
	return out
	
def main():	
	'''Create Looger'''
	logging.basicConfig()
	logger = logging.getLogger(__name__)
	logger.setLevel(logging.DEBUG)
	
	'''reading the config filename'''
	parser=argument_parse()
	arg=parser.parse_args()
	curr_path = os.path.abspath(__file__)
	
	'''Params'''
	input_file=arg.input
	out_file=arg.output_file
	temp_dir=arg.temp_dir
	current_file=input_file
	logger.info("Input File: "+input_file)
	logger.info("Input Output file: "+out_file)
	logger.info("Input temp directory: "+temp_dir)
	
	if arg.capture_kit != None :
		capture_kit=arg.capture_kit
		logger.info("Input capture_kit: "+capture_kit)
	
	if arg.dac_gap != None :
		dac_gap = arg.dac_gap
		logger.info("Input DAC,GAP,Delly blacklist file: "+dac_gap)

	if arg.rep_mas != None:
		rep_mas = arg.rep_mas
		logger.info("Input Rep_mask file: "+rep_mas)

	if arg.other_sample != None:
		other_sample = arg.other_sample
		logger.info("Other Sample: "+other_sample)
	
	if arg.normal_samp != None:
		normal_samp = arg.normal_samp
		logger.info("Normal Sample: "+normal_samp)

	if arg.level1_bp != None:
		level1_bp = arg.level1_bp
		logger.info("level1_bp file: "+level1_bp)
		
	if arg.gtf_file != None:
		gtf_file = arg.gtf_file
		logger.info("gtf_file file: "+gtf_file)
	
	if arg.chr_names != None:
		chr_names = arg.chr_names
		logger.info("Chr_names file: "+chr_names)		
		
	if arg.var_files != None:
		var_files = arg.var_files
		logger.info("var_files file: "+var_files)
		
	if arg.paper_freq_pairs != None:
		paper_freq_pairs = arg.paper_freq_pairs
		logger.info("paper_freq_pairs file: "+paper_freq_pairs)

	if arg.gene_expr != None:
		gene_expr = arg.gene_expr
		logger.info("gene_expr: "+gene_expr)

	if arg.fusion_out_dir != None:
		fusion_out_dir = arg.fusion_out_dir
		logger.info("fusion_out_dir: "+fusion_out_dir)

	if arg.detected_trans != None:
		detected_trans = arg.detected_trans
		logger.info("detected_trans file: "+detected_trans)			

	''''''
	gene_list=['bcl6','myc','bcl2']
	gene_start_list = [187721377,127735434,63123346]
	gene_stop_list = [187745725,127742951,63320128]
	if chr_filter:
		chr_list = ['chr3','chr8','chr18']
	else:
		chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	chr_list_all = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	
	'''printing the config param'''	
	if not os.path.exists(temp_dir):
		os.makedirs(temp_dir)
	if capture_kit != None or dac_gap != None or rep_mas != None:
		'''Reading the input files'''
		read_input = open(input_file)
		header = read_input.readline()
		write_tmp = open(temp_dir+'/tmp.bed', mode='wt')
		linenum = 0
		for line in read_input:
			line=line.strip()
			rw_lst = line.split("\t")
			linenum = linenum+1
			write_tmp.write(rw_lst[5]+"\t"+str(int(rw_lst[6])-1)+"\t"+rw_lst[6]+"\t"+str(linenum)+"\t"+str(1)+"\n")
			write_tmp.write(rw_lst[7]+"\t"+str(int(rw_lst[8])-1)+"\t"+rw_lst[8]+"\t"+str(linenum)+"\t"+str(2)+"\n")
		read_input.close()
		write_tmp.close()
	
		dict_bed={}
		if capture_kit != None :
			shutil.copyfile(capture_kit, temp_dir+'/tmp_cap.bed') 
			
			'''Run intersect bed'''
			intersectBed_run(temp_dir+'/tmp.bed',temp_dir+'/tmp_cap.bed',temp_dir+'/tmp_out.bed',temp_dir)
			
			'''reading the intersectBed output file'''
			
			read_tmp = open(temp_dir+'/tmp_out.bed')
			for line in read_tmp:
				line = line.strip()
				rw_lst = line.split("\t")
				'''ignoring rows with no overlap'''
				if rw_lst[5] != '.':
					'''selecting chr and pos with capture kit intersect''' 
					dict_bed[rw_lst[0]+' '+rw_lst[2]]=1
			read_tmp.close()

			'''removing temp files'''
			os.remove(temp_dir+'/tmp_cap.bed')
			os.remove(temp_dir+'/tmp_out.bed')
		
		dict_dac={}
		if dac_gap != None :
			'''Preparing DAC exclude file'''
			shutil.copyfile(dac_gap, temp_dir+'/tmp_dac.bed') 
				
			'''Run intersect bed'''
			intersectBed_run(temp_dir+'/tmp.bed',temp_dir+'/tmp_dac.bed',temp_dir+'/tmp_out.bed',temp_dir)
			
			'''reading the intersectBed output DAC file'''
			
			read_dac_out = open(temp_dir +'/tmp_out.bed')
			for line in read_dac_out:
				line = line.strip()
				rw_lst = line.split("\t")
				if int(rw_lst[9]) > 0:
					dict_dac[rw_lst[0]+' '+rw_lst[2]]=1
			read_dac_out.close()
			
			'''removing temp files'''
			os.remove(temp_dir+'/tmp_dac.bed')
			os.remove(temp_dir+'/tmp_out.bed')
		
		dict_rep={}
		if rep_mas != None:
			'''Preparing DAC exclude file'''
			shutil.copyfile(rep_mas, temp_dir+'/tmp_rep.bed') 
				
			'''Run intersect bed'''
			intersectBed_run(temp_dir+'/tmp.bed',temp_dir+'/tmp_rep.bed',temp_dir+'/tmp_out.bed',temp_dir)
			
			'''reading the intersectBed output RepeatMasker file'''	
			read_rep_out = open(temp_dir +'/tmp_out.bed')
			for line in read_rep_out:
				line = line.strip()
				rw_lst = line.split("\t")
				if int(rw_lst[9]) > 0:
					dict_rep[rw_lst[0]+' '+rw_lst[2]]=1
			read_rep_out.close()
			
			'''removing temp files'''
			os.remove(temp_dir+'/tmp_rep.bed')
			os.remove(temp_dir+'/tmp_out.bed')
		if capture_kit != None or dac_gap != None or rep_mas != None:
			os.remove(temp_dir+'/tmp.bed')
			'''Reading input'''
			read_input = open(input_file)
			write_out = open(out_file, mode='wt')
			header = read_input.readline()
			header = header.strip()
			write_out.write(header)
			if capture_kit != None:
				write_out.write("\tBP1_CapKit\tBP2_CapKit")
			
			if dac_gap != None:
				write_out.write("\tBP1_DAC_DellyEx_Gap\tBP2_DAC_DellyEx_Gap")
			
			if rep_mas != None:
				write_out.write("\tBP1_RepMaskSegDup\tBP2_RepMaskSegDup")
			
			write_out.write("\n")
			for line in read_input:
				line = line.strip()
				rw_lst = line.split("\t")
				
				write_out.write(line)
				
				if capture_kit != None :
					val1="0"
					'''if bp1 in capture kit dict'''
					if rw_lst[5]+' '+rw_lst[6] in dict_bed:
						val1="1"
					val2="0"
					'''if bp2 in capture kit dict'''
					if rw_lst[7]+' '+rw_lst[8] in dict_bed:
						val2="1"
					write_out.write("\t"+val1+"\t"+val2)
				
				if dac_gap != None :
					val3="0"
					'''if bp1 in dac kit dict'''
					if rw_lst[5]+' '+rw_lst[6] in dict_dac:
						val3="1"	
					val4="0"
					'''if bp2 in dac kit dict'''
					if rw_lst[7]+' '+rw_lst[8] in dict_dac:
						val4="1"
					write_out.write("\t"+val3+"\t"+val4)
				if rep_mas != None:
					val5="0"
					'''if bp1 in rep dict'''
					if rw_lst[5]+' '+rw_lst[6] in dict_rep:
						val5="1"
					val6="0"
					'''if bp2 in rep dict'''
					if rw_lst[7]+' '+rw_lst[8] in dict_rep:
						val6="1"
					write_out.write("\t"+val5+"\t"+val6)
				
				write_out.write("\n")
			read_input.close()
			write_out.close()
		input_file=out_file
		del dict_dac
		del dict_bed
		del dict_rep
	'''Same BP in Other samples'''
	if other_sample==True:
		'''Preparing file for other sample step'''
		shutil.copyfile(input_file, temp_dir+'/tmpfile')
		
		'''Reading input file'''
		fobj = open(temp_dir+'/tmpfile')
		header = fobj.readline()
		dict_othersamp_both = {}
		dict_othersamp_first = {}
		dict_othersamp_sec = {}
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			samp = rw_lst[1]
			chr1 = rw_lst[5]
			pos1 = rw_lst[6]
			chr2 = rw_lst[7]
			pos2 = rw_lst[8]
			caller = rw_lst[12]
			
			'''sort primary and mate chr to account for same bp with different chr order'''
			x=[chr1,chr2]
			x.sort()
			if x[0] != chr1:
				tmp1 = chr1
				tmp2 = pos1
				chr1 = chr2
				pos1 = pos2
				chr2 = tmp1
				pos2 = tmp2
			'''saving primate bp information in the dictionary'''
			if not chr1+'__'+pos1 in dict_othersamp_first:
				dict_othersamp_first[chr1+'__'+pos1]=samp
			else:
				if not samp in dict_othersamp_first[chr1+'__'+pos1]:
					dict_othersamp_first[chr1+'__'+pos1]=dict_othersamp_first[chr1+'__'+pos1]+'__'+samp
			'''saving mate bp information in the dictionary'''		
			if not chr2+'__'+pos2 in dict_othersamp_sec:
				dict_othersamp_sec[chr2+'__'+pos2]=samp
			else:
				if not samp in dict_othersamp_sec[chr2+'__'+pos2]:
					dict_othersamp_sec[chr2+'__'+pos2]=dict_othersamp_sec[chr2+'__'+pos2]+'__'+samp
			'''saving both bp end information in the dictionary'''
			if not chr1+'__'+pos1+'__'+chr2+'__'+pos2 in dict_othersamp_both:
				dict_othersamp_both[chr1+'__'+pos1+'__'+chr2+'__'+pos2]=samp
			else:
				if not samp in dict_othersamp_both[chr1+'__'+pos1+'__'+chr2+'__'+pos2]:
					dict_othersamp_both[chr1+'__'+pos1+'__'+chr2+'__'+pos2]=dict_othersamp_both[chr1+'__'+pos1+'__'+chr2+'__'+pos2]+'__'+samp
		fobj.close()

		'''write output file'''
		fobj = open(temp_dir+'/tmpfile')
		myfile = open(out_file, mode='wt')
		header = fobj.readline()
		header = header.strip()
		myfile.write(header+"\tOther_samp_num_pos1\tOther_samp_ID_pos1\tOther_samp_num_pos2\tOther_samp_ID_pos2\tOther_samp_num_pos_both\tOther_samp_ID_pos_both\n")
		linenum = 0
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			samp = rw_lst[1]
			chr1 = rw_lst[5]
			pos1 = rw_lst[6]
			chr2 = rw_lst[7]
			pos2 = rw_lst[8]
			caller = rw_lst[12]
			pos1_num=0
			pos2_num=0
			pos1_samp='NA'
			pos2_samp='NA'
			pos_num=0
			pos_samp='NA'
			'''checking if primary bp info in dict'''
			if chr1+'__'+pos1 in dict_othersamp_first:
				samp1 = dict_othersamp_first[chr1+'__'+pos1].split('__')
				if samp in samp1:
					samp1.remove(samp)
				samp11=list(set(samp1))
				pos1_num = len(samp11)
				if pos1_num > 0:
					pos1_samp = str.join("__",samp11)
			'''checking if mate bp info in dict'''
			if chr2+'__'+pos2 in dict_othersamp_sec:
				samp2 = dict_othersamp_sec[chr2+'__'+pos2].split('__')
				if samp in samp2:
					samp2.remove(samp)
				samp22=list(set(samp2))
				pos2_num = len(samp22)
				if pos2_num > 0:
					pos2_samp = str.join("__",samp2)
			'''checking if both ends bp info in dict'''
			if chr1+'__'+pos1+'__'+chr2+'__'+pos2 in dict_othersamp_both:
				samp3 = dict_othersamp_both[chr1+'__'+pos1+'__'+chr2+'__'+pos2].split('__')
				if samp in samp3:
					samp3.remove(samp)
				samp33=list(set(samp3))
				pos_num = len(samp3)
				if pos_num > 0:
					pos_samp = str.join("__",samp33)
			myfile.write(line+"\t"+str(pos1_num)+"\t"+pos1_samp+"\t"+str(pos2_num)+"\t"+pos2_samp+"\t"+str(pos_num)+"\t"+pos_samp+"\n")	
		fobj.close()
		myfile.close()
		os.remove(temp_dir+'/tmpfile')
		input_file = out_file
		del dict_othersamp_both
		del dict_othersamp_first
		del dict_othersamp_sec
		
	'''Normal samp freq annotation'''
	if normal_samp != None:
		shutil.copyfile(input_file, temp_dir+'/tmpfile')
		'''Preparing file for Normal sample step'''
		'''creating the bed coordinates for input file'''
		fobj = open(temp_dir+'/tmpfile')
		myfile = open(temp_dir+'/tmp1.bed', mode='wt')
		header = fobj.readline()
		distance=300
		linenum = 0
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			linenum = linenum+1
			tmp1=int(rw_lst[6])-distance
			if tmp1 < 1:
				tmp1=1
			tmp2=int(rw_lst[6])+distance                
			myfile.write(rw_lst[5]+"\t"+str(tmp1)+"\t"+str(tmp2)+"\t"+str(linenum)+"\t"+rw_lst[1]+'__'+str(1)+"\n")
			tmp1=int(rw_lst[8])-distance
			if tmp1 < 1:
				tmp1=1			
			tmp2=int(rw_lst[8])+distance
			myfile.write(rw_lst[7]+"\t"+str(tmp1)+"\t"+str(tmp2)+"\t"+str(linenum)+"\t"+rw_lst[1]+'__'+str(2)+"\n")
		fobj.close()
		myfile.close()
		'''creating the bed coordinates for normal sample breakpoint file'''
		fobj = open(normal_samp)
		myfile = open(temp_dir+'/tmp2.bed', mode='wt')
		header = fobj.readline()
		linenum = 0
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			linenum = linenum+1
			tmp1=int(rw_lst[6])-distance
			if tmp1 < 1:
				tmp1=1
			tmp2=int(rw_lst[6])+distance                
			myfile.write(rw_lst[5]+"\t"+str(tmp1)+"\t"+str(tmp2)+"\t"+str(linenum)+"\t"+rw_lst[1]+'__'+str(1)+"\n")
			tmp1=int(rw_lst[8])-distance
			if tmp1 < 1:
				tmp1=1			
			tmp2=int(rw_lst[8])+distance
			myfile.write(rw_lst[7]+"\t"+str(tmp1)+"\t"+str(tmp2)+"\t"+str(linenum)+"\t"+rw_lst[1]+'__'+str(2)+"\n")
		fobj.close()
		myfile.close()
		
		'''Run intersect bed'''
		intersectBed_run(temp_dir+'/tmp1.bed',temp_dir+'/tmp2.bed',temp_dir+'/tmp_out.bed',temp_dir)
		
		'''identifying the breakpoint intersected with normal sample bps and storing them in dict'''	
		dict_norm1={}
		dict_norm2={}
		final={}
		vl=""
		fobj = open(temp_dir +'/tmp_out.bed')
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			if rw_lst[5] != '.':
				line_num1 = rw_lst[3]
				rw_lst1 = rw_lst[4].split("__")
				samp1=rw_lst1[0]
				bp1=rw_lst1[1]
				line_num2 = rw_lst[8]
				rw_lst2 = rw_lst[9].split("__")
				samp2=rw_lst2[0]
				bp2=rw_lst2[1]
				if rw_lst[5] != '.' and samp1 != samp2:
					if int(samp1) > int(samp2):
						key = line_num1+'__'+samp1+'__'+line_num2+'__'+samp2
						val = bp1+bp2
					else:
						key = line_num2+'__'+samp2+'__'+line_num1+'__'+samp1
						val = bp2+bp1
					if not key in dict_norm1:	
						dict_norm1[key]=val
					else:
						val1=dict_norm1[key]
						if (val == "11" and val1 == "22") or (val == "22" and val1 == "11") or (val == "12" and val1 == "21") or (val == "21" and val1 == "12"):
							ky = line_num1+'__'+samp1
							vl = samp2
							if ky in dict_norm2:
								vl1 = dict_norm2[ky]
								if not samp2 in vl1:
									dict_norm2[ky] = dict_norm2[ky]+'__'+samp2
							else:
								dict_norm2[ky] = vl	
		fobj.close()
		os.remove(temp_dir+'/tmp1.bed')
		os.remove(temp_dir+'/tmp2.bed')
		os.remove(temp_dir+'/tmp_out.bed')
		'''Parsing input file and add the normal samp freq'''
		fobj = open(temp_dir+'/tmpfile')
		myfile = open(out_file, mode='wt')
		header = fobj.readline()
		header = header.strip()
		myfile.write(header+"\tNormal_samp_num\tNormal_samp_ID\n")
		linenum = 0
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			linenum = linenum+1
			if str(linenum)+'__'+rw_lst[1] in dict_norm2:
				oth_samp=dict_norm2[str(linenum)+'__'+rw_lst[1]]
				list_othsamp=dict_norm2[str(linenum)+'__'+rw_lst[1]].split('__')
				list_othsamp_len=len(list_othsamp)
			else:
				oth_samp="NA"
				list_othsamp_len=0
			myfile.write(line+"\t"+str(list_othsamp_len)+"\t"+oth_samp+"\n")
				
		fobj.close()
		myfile.close()
		os.remove(temp_dir+'/tmpfile')
		input_file = out_file
		del dict_norm1
		del dict_norm2
		
	'''L1 region annotation'''
	if level1_bp != None:
		shutil.copyfile(input_file, temp_dir+'/tmpfile')
		
		'''Reading the input files'''
		read_input = open(temp_dir+'/tmpfile')
		header = read_input.readline()
		write_tmp = open(temp_dir+'/tmp.bed', mode='wt')
		linenum = 0
		for line in read_input:
			line=line.strip()
			rw_lst = line.split("\t")
			linenum = linenum+1
			write_tmp.write(rw_lst[5]+"\t"+str(int(rw_lst[6])-1)+"\t"+rw_lst[6]+"\t"+str(linenum)+"\t"+str(1)+"\n")
			#write_tmp.write(rw_lst[7]+"\t"+str(int(rw_lst[8])-1)+"\t"+rw_lst[8]+"\t"+str(linenum)+"\t"+str(2)+"\n")
		read_input.close()
		write_tmp.close()
		
		shutil.copyfile(level1_bp, temp_dir+'/tmp_l1.bed') 
			
		'''Run intersect bed'''
		intersectBed_run(temp_dir+'/tmp.bed',temp_dir+'/tmp_l1.bed',temp_dir+'/tmp_out.bed',temp_dir)
		
		dict_rep={}
		'''reading the intersectBed output RepeatMasker file'''	
		read_rep_out = open(temp_dir +'/tmp_out.bed')
		for line in read_rep_out:
			line = line.strip()
			rw_lst = line.split("\t")
			if int(rw_lst[9]) > 0:
				dict_rep[rw_lst[0]+' '+rw_lst[2]]=1
		read_rep_out.close()
		
		'''removing temp files'''
		os.remove(temp_dir+'/tmp.bed')
		os.remove(temp_dir+'/tmp_l1.bed')
		os.remove(temp_dir+'/tmp_out.bed')

		fobj = open(temp_dir+'/tmpfile')
		write_out = open(out_file, mode='wt')
		header = fobj.readline()
		header = header.strip()
		write_out.write(header+"\tBP1_L1_Genic\n")
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			level1_col=len(rw_lst)
			val1="0"
			if rw_lst[5]+' '+rw_lst[6] in dict_rep:
				val1="1"
			write_out.write(line+"\t"+val1+"\n")
		fobj.close()
		write_out.close()
		del dict_rep
		os.remove(temp_dir+'/tmpfile')
		input_file = out_file

	'''Rachel's method detected_trans file'''
	if detected_trans != None:
		shutil.copyfile(input_file, temp_dir+'/tmpfile')
		'''Reading the input files'''
		read_input = open(temp_dir+'/tmpfile')
		header = read_input.readline()
		write_tmp = open(temp_dir+'/tmp.bed', mode='wt')
		linenum = 0
		for line in read_input:
			line=line.strip()
			rw_lst = line.split("\t")
			linenum = linenum+1
			write_tmp.write(rw_lst[1]+'__'+rw_lst[5]+"\t"+str(int(rw_lst[6])-1)+"\t"+rw_lst[6]+"\t"+str(linenum)+"\t"+str(1)+"\n")
		read_input.close()
		write_tmp.close()
		
		'''Reading the detected_trans file'''
		read_input = open(detected_trans)
		header = read_input.readline()
		write_tmp = open(temp_dir+'/tmp1.bed', mode='wt')
		linenum = 0
		for line in read_input:
			line=line.strip()
			rw_lst = line.split("\t")
			linenum = linenum+1
			write_tmp.write(rw_lst[0]+'__'+rw_lst[1]+"\t"+str(int(rw_lst[2])-1)+"\t"+rw_lst[3]+"\t"+line+"\n")
		read_input.close()
		write_tmp.close()
		
		'''Run intersect bed'''
		intersectBed_run(temp_dir+'/tmp.bed',temp_dir+'/tmp1.bed',temp_dir+'/tmp_out.bed',temp_dir)
		
		'''reading the intersectBed output DAC file'''
		
		read_dac_out = open(temp_dir +'/tmp_out.bed')
		dict_det_trans={}
		for line in read_dac_out:
			line = line.strip()
			rw_lst = line.split("\t")
			if int(rw_lst[15]) > 0:
				dict_det_trans[rw_lst[0]+' '+rw_lst[2]]=rw_lst[9]+"\t"+rw_lst[10]+"\t"+rw_lst[11]+"\t"+rw_lst[12]+"\t"+rw_lst[13]+"\t"+rw_lst[14]
		read_dac_out.close()
		
		'''Reading input'''
		read_input = open(temp_dir+'/tmpfile')
		write_out = open(out_file, mode='wt')
		header = read_input.readline()
		header = header.strip()
		write_out.write(header+"\tdetected_trans_chrom\tdetected_trans_start\tdetected_trans_stop\tdetected_trans_all\tdetected_trans_max_pct\tdetected_trans_max_pct_chrom\n")
		for line in read_input:
			line = line.strip()
			rw_lst = line.split("\t")
			write_out.write(line)
			val1="NA\tNA\tNA\tNA\tNA\tNA"
			'''if bp1 in capture kit dict'''
			if rw_lst[1]+'__'+rw_lst[5]+' '+rw_lst[6] in dict_det_trans:
				val1=dict_det_trans[rw_lst[1]+'__'+rw_lst[5]+' '+rw_lst[6]]
			write_out.write("\t"+val1+"\n")
		read_input.close()
		write_out.close()
			
		'''removing temp files'''
		os.remove(temp_dir+'/tmp.bed')
		os.remove(temp_dir+'/tmp1.bed')
		os.remove(temp_dir+'/tmp_out.bed')
		del dict_det_trans	
		os.remove(temp_dir+'/tmpfile')
		input_file = out_file
		
		
	'''GTF File annotation'''
	if gtf_file != None and chr_names != None:
		shutil.copyfile(input_file, temp_dir+'/tmpfile')	
		gene_list=['BCL6','MYC','BCL2']
		dist=1000000
		chrlist = {}
		'''loading chr names'''
		fobj = open(chr_names)
		for file0 in fobj:
			file0=file0.strip()
			rw_lst = file0.split("_")
			if len(rw_lst)>1:
				rw_lst[1]=rw_lst[1].replace("v1",".1")
				rw_lst[1]=rw_lst[1].replace("v2",".2")
				rw_lst[1]=rw_lst[1].replace("v3",".3")
				chrlist[rw_lst[1]]=file0
			else:
				chrlist[file0]=file0
		fobj.close()
			
		'''Gene names start and stop coord from gtf file'''	
		fobj = open(gtf_file)
		myfile = open(temp_dir+'/tmp1.bed', mode='wt')
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			if not line.startswith("##") and rw_lst[2] == "gene":
				gn=re.search('gene_name "(.+?)"', rw_lst[8]).group(1)
				chr=rw_lst[0]
				if chr in chrlist:
					chr = chrlist[chr]
				start=int(rw_lst[3])-dist
				stop=int(rw_lst[4])+dist
				if start<1:
					start=1
				if stop<1:
					stop=1
				myfile.write(chr+"\t"+str(start)+"\t"+str(stop)+"\t"+str(gn)+'___'+rw_lst[3]+'___'+rw_lst[4]+"\n")
		fobj.close()
		myfile.close()
		
		'''BP coords from input file'''	
		fobj = open(temp_dir+'/tmpfile')
		myfile = open(temp_dir+'/tmp2.bed', mode='wt')
		header = fobj.readline()
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			myfile.write(rw_lst[5]+"\t"+str(int(rw_lst[6])-1)+"\t"+rw_lst[6]+"\n")
			myfile.write(rw_lst[7]+"\t"+str(int(rw_lst[8])-1)+"\t"+rw_lst[8]+"\n")
		fobj.close()
		myfile.close()

		'''Run intersect bed'''
		intersectBed_run(temp_dir+'/tmp2.bed',temp_dir+'/tmp1.bed',temp_dir+'/tmp_out.bed',temp_dir)
		
		dict_bed={}
		dict_bed_dist={}
		dict_bed2={}
		dict_bed_dist2={}
		dict_dup={}
		'''Loading the intersect regions and adding the gene information to dict'''
		fobj = open(temp_dir +'/tmp_out.bed')
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			if int(rw_lst[7]) > 0:
				start1 = int(rw_lst[1])
				stop1 = int(rw_lst[2])
				rw_lst1 = rw_lst[6].split("___")
				start2 = int(rw_lst1[1])
				stop2 = int(rw_lst1[2])
				gn = rw_lst1[0]
				key = rw_lst[0]+' '+rw_lst[2]
				dist = 0
				if (start1 >= start2 and start1 <=stop2) or (stop1 >= start2 and stop1 <=stop2):
					dist=0
				else:
					dist=min([abs(start1-start2),abs(start1-stop2),abs(stop1-start2),abs(stop1-stop2)])
				if not key+' '+gn in dict_dup:
					dict_dup[key+' '+gn]=1
					if key in dict_bed_dist:
						dict_bed[key]=dict_bed[key]+','+gn
						dict_bed_dist[key]=dict_bed_dist[key]+','+str(dist)
					else:	
						dict_bed[key]=gn
						dict_bed_dist[key]=str(dist)
		fobj.close()

		'''removing temp files'''
		os.remove(temp_dir+'/tmp1.bed')
		os.remove(temp_dir+'/tmp2.bed')
		os.remove(temp_dir+'/tmp_out.bed')
		
		'''Output file with gene annotations'''
		fobj = open(temp_dir+'/tmpfile')
		myfile = open(out_file, mode='wt')
		header = fobj.readline()
		header = header.strip()
		myfile.write(header+"\tBP1_Gene\tBP2_Gene\tBP1_Dist\tBP2_Dist\n")
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			'''getting col number for gene1 and gene2 and will be used later in the next few steps'''
			gene_1_col_num=len(rw_lst)
			gene_2_col_num=len(rw_lst)+1
			gn1=["NA"]
			gn2=["NA"]
			if rw_lst[5]+' '+rw_lst[6] in dict_bed:
				gn1= dict_bed[rw_lst[5]+' '+rw_lst[6]].split(",")
				dst_list1= dict_bed_dist[rw_lst[5]+' '+rw_lst[6]].split(",")
				dst_list1 = list(map(int, dst_list1))
				sort_index = np.argsort(np.array(dst_list1))
				gna=[gn1[i] for i in sort_index if dst_list1[i] == 0]
				dista=[dst_list1[i] for i in sort_index if dst_list1[i] == 0]
				gnb=[gn1[i] for i in sort_index if dst_list1[i] != 0][:2]
				distb=[dst_list1[i] for i in sort_index if dst_list1[i] != 0][:2]
				if len(gna)<2:
					gn1=gna+gnb
					dist1=dista+distb
				else:
					gn1=gna
					dist1=dista
			if rw_lst[7]+' '+rw_lst[8] in dict_bed:
				gn2= dict_bed[rw_lst[7]+' '+rw_lst[8]].split(",")
				dst_list1= dict_bed_dist[rw_lst[7]+' '+rw_lst[8]].split(",")
				dst_list1 = list(map(int, dst_list1))
				sort_index = np.argsort(np.array(dst_list1))
				gna=[gn2[i] for i in sort_index if dst_list1[i] == 0]
				dista=[dst_list1[i] for i in sort_index if dst_list1[i] == 0]
				gnb=[gn2[i] for i in sort_index if dst_list1[i] != 0][:2]
				distb=[dst_list1[i] for i in sort_index if dst_list1[i] != 0][:2]
				if len(gna)<2:
					gn2=gna+gnb
					dist2=dista+distb
				else:
					gn2=gna
					dist2=dista
			for x in range(len(gn1)):
				for y in range(len(gn2)):
					myfile.write(line+"\t"+gn1[x]+"\t"+gn2[y]+"\t"+str(dist1[x])+"\t"+str(dist2[y])+"\n")
		fobj.close()
		myfile.close()
		del dict_bed
		del dict_bed_dist
		del dict_bed2
		del dict_bed_dist2
		del dict_dup
		os.remove(temp_dir+'/tmpfile')
		input_file = out_file
	'''Variant file annotation'''
	if gtf_file != None and chr_names != None and var_files != None:
		shutil.copyfile(input_file, temp_dir+'/tmpfile')
		'''loading chr names'''
		chrlist = {}
		fobj = open(chr_names)
		for file0 in fobj:
			file0=file0.strip()
			rw_lst = file0.split("_")
			if len(rw_lst)>1:
				rw_lst[1]=rw_lst[1].replace("v1",".1")
				rw_lst[1]=rw_lst[1].replace("v2",".2")
				rw_lst[1]=rw_lst[1].replace("v3",".3")
				chrlist[rw_lst[1]]=file0
			else:
				chrlist[file0]=file0
		fobj.close()
		
		'''Gene names start and stop coord from gtf file'''	
		dist=0	
		fobj = open(gtf_file)
		myfile = open(temp_dir+'/tmp.bed', mode='wt')
		for line in fobj:
			line = line.strip()
			arr = line.split("\t")
			if not line.startswith("##") and arr[2] == "exon":
				gn=re.search('gene_name "(.+?)"', arr[8]).group(1)
				chr=arr[0]
				if chr in chrlist:
					chr = chrlist[chr]
				start=int(arr[3])-dist
				stop=int(arr[4])+dist
				if start<1:
					start=1
				if stop<1:
					stop=1
				myfile.write(chr+"\t"+str(start)+"\t"+str(stop)+"\t"+str(gn)+"\n")
		fobj.close()
		myfile.close()
		
		'''Reading summarized vcf files'''
		vcffiles=os.listdir(var_files)
		myfile = open(temp_dir+'/tmp_l1.bed', mode='wt')
		for line in vcffiles:
			lst = line.split('.')
			samp = lst[0]
			fobj = open(var_files+'/'+line)
			header = fobj.readline()
			for line in fobj:
				line = line.strip()
				if line!="":
					rw_lst = line.split("\t")
					rw_lst = rw_lst[:3]+rw_lst
					'''selecting snvs passed GATK best practices cutoff'''
					if len(rw_lst[6]) < 2 and len(rw_lst[7]) < 2 and (rw_lst[8] == '.' or float(rw_lst[8])>30) \
					and (rw_lst[25] == '.' or float(rw_lst[25])<3)  \
					and (rw_lst[17] == '.' or float(rw_lst[17])<60) \
					and (rw_lst[21] == '.' or float(rw_lst[21])>40) \
					and (rw_lst[22] == '.' or float(rw_lst[22])> -12.5) \
					and (rw_lst[24] == '.' or float(rw_lst[24]) > -8)  \
					and ( rw_lst[8] == '.' or rw_lst[14] == '.' or (float(rw_lst[8])/float(rw_lst[14]))>2 )  \
					and ( rw_lst[33] == '.' or float(rw_lst[33]) < 0.01) \
					and (rw_lst[41] == '.' or float(rw_lst[41])<0.01)  \
					and ( rw_lst[50] == '.' or float(rw_lst[50]) <0.01) \
					and (rw_lst[60] == '.' or float(rw_lst[60])<0.01)  \
					and (rw_lst[104] == '.' or float(rw_lst[104])>10):
						tmp = int(rw_lst[1])-1
						myfile.write(rw_lst[0]+"\t"+str(tmp)+"\t"+rw_lst[1]+"\t"+samp+"\n")
					'''selecting indels passed GATK best practices cutoff'''
					if (len(rw_lst[6]) > 1 or len(rw_lst[7]) > 1) \
					 and (rw_lst[8] == '.' or float(rw_lst[8]) > 30) \
					 and (rw_lst[17] == '.' or float(rw_lst[17]) < 200) \
					 and (rw_lst[24] == '.' or float(rw_lst[24]) > -20) \
					 and (rw_lst[8] == '.' or rw_lst[14] == '.' or (float(rw_lst[8]) / float(rw_lst[14])) > 2) \
					 and (rw_lst[33] == '.' or float(rw_lst[33]) < 0.01) \
					 and (rw_lst[41] == '.' or float(rw_lst[41]) < 0.01) \
					 and (rw_lst[50] == '.' or float(rw_lst[50]) < 0.01) \
					 and (rw_lst[60] == '.' or float(rw_lst[60]) < 0.01) \
					 and (rw_lst[104] == '.' or float(rw_lst[104]) > 10):
						tmp = int(rw_lst[1])-1
						myfile.write(rw_lst[0]+"\t"+str(tmp)+"\t"+rw_lst[1]+"\t"+samp+"\n")
			fobj.close()	
		myfile.close()	
		
		'''Run intersect bed'''
		intersectBed_run(temp_dir+'/tmp.bed',temp_dir+'/tmp_l1.bed',temp_dir+'/tmp_out.bed',temp_dir)
		
		'''reading the intersectBed output variant file'''	
		dict_bed={}
		dict_gene={}
		fobj = open(temp_dir+'/tmp_out.bed')
		for line in fobj:
			line = line.strip()
			arr = line.split("\t")
			if int(arr[8]) > 0:
				gene = arr[3]
				samp = arr[7]
				chrom = arr[4]
				pos = arr[6]
				key1 = gene+'__'+samp+'__'+chrom+'__'+pos
				key2 = samp+'__'+gene
				if not key1 in dict_bed:
					dict_bed[key1]=1
					if key2 in dict_gene:
						dict_gene[key2]=dict_gene[key2]+1
					else:	
						dict_gene[key2]=1
		fobj.close()
		
		'''removing temp files'''
		os.remove(temp_dir+'/tmp.bed')
		os.remove(temp_dir+'/tmp_l1.bed')
		os.remove(temp_dir+'/tmp_out.bed')
		
		'''Parsing input file and add the num variants per gene'''
		fobj = open(temp_dir+'/tmpfile')
		myfile = open(out_file, mode='wt')
		header = fobj.readline()
		header = header.strip()
		myfile.write(header+"\tBP1_NumVar\tBP2_NumVar\n")
		for line in fobj:
			line = line.strip()
			arr = line.split("\t")
			samp=arr[1]
			gn1=samp+'__'+arr[gene_1_col_num]
			gn2=samp+'__'+arr[gene_2_col_num]
			val1 = 0
			if gn1 in dict_gene:
				val1 = dict_gene[gn1]
			val2 = 0
			if gn2 in dict_gene:
				val2 = dict_gene[gn2]
			
			myfile.write(line+"\t"+str(val1)+"\t"+str(val2)+"\n")
		fobj.close()
		myfile.close()
		del dict_bed
		del dict_gene
		os.remove(temp_dir+'/tmpfile')
		input_file = out_file

	'''paper_freq_pairs file'''
	if gtf_file != None and chr_names != None and paper_freq_pairs != None:
		shutil.copyfile(input_file, temp_dir+'/tmpfile')
		'''Reading the paper freq file'''
		dict_mate_freq = {}
		fobj = open(paper_freq_pairs)
		header = fobj.readline()
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			dict_mate_freq[rw_lst[0]+'__'+rw_lst[1]]=1
		fobj.close()				
		
		fobj = open(temp_dir+'/tmpfile')
		myfile = open(out_file, mode='wt')
		header = fobj.readline()
		header = header.strip()
		myfile.write(header+"\tChong_et_al_freq_mate\n")
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			val3=0
			if rw_lst[gene_1_col_num]+'__'+rw_lst[gene_2_col_num] in dict_mate_freq:
				val3 = dict_mate_freq[rw_lst[gene_1_col_num]+'__'+rw_lst[gene_2_col_num]]
			
			myfile.write(line+"\t"+str(val3)+"\n")
		fobj.close()
		myfile.close()

		del dict_mate_freq
		os.remove(temp_dir+'/tmpfile')
		input_file = out_file

	'''gene expr file'''
	if gtf_file != None and chr_names != None and gene_expr != None:
		shutil.copyfile(input_file, temp_dir+'/tmpfile')
		
		gnexp = {}
		'''Reading gene expression file'''
		fobj = open(gene_expr)
		header = fobj.readline()
		for file0 in fobj:
			file0=file0.strip()
			rw_lst = file0.split("\t")
			samp = rw_lst[7]
			genename = rw_lst[0]
			fc = rw_lst[1]
			pval = rw_lst[4]
			gnexp[samp+'__'+genename] = fc+"\t"+pval
		fobj.close()
		'''Reading input file'''
		fobj = open(temp_dir+'/tmpfile')
		myfile = open(out_file, mode='wt')
		header = fobj.readline()
		header = header.strip()
		myfile.write(header+"\tBP1_FC\tBP1_Pval\tBP2_FC\tBP2_Pval\n")
		for line in fobj:
			line = line.strip()
			rw_lst = line.split("\t")
			samp=rw_lst[1]
			#flag=rw_lst[87]
			gn1=rw_lst[gene_1_col_num]
			gn2=rw_lst[gene_2_col_num]
			val1 = "NA\tNA"
			if samp+'__'+gn1 in gnexp:
				val1 = gnexp[samp+'__'+gn1]
			val2 = "NA\tNA"
			if samp+'__'+gn2 in gnexp:
				val2 = gnexp[samp+'__'+gn2]				
			myfile.write(line+"\t"+str(val1)+"\t"+str(val2)+"\n")
		fobj.close()
		myfile.close()
		del gnexp
		os.remove(temp_dir+'/tmpfile')
		input_file = out_file
	
	'''Fusions'''
	if gtf_file != None and chr_names != None and fusion_out_dir!=None:
		shutil.copyfile(input_file, temp_dir+'/tmpfile')
		'''list fusion catcher output files'''
		fusion_files=os.listdir(fusion_out_dir)
		dict_fusion_genes={}
		dist=100
		hed = 0
		'''parsing each file and creating dictionary with fusions per sample'''
		for file in fusion_files:
			lst = file.split('.')
			samp = lst[0]
			'''reading each fusion out file'''
			fobj = open(fusion_out_dir+'/'+file)
			header = fobj.readline()
			linenum = 0
			
			for line in fobj:
				linenum = linenum+1
				line = line.strip()
				rw_lst = line.split("\t")
				gn1 = rw_lst[0]
				gn2 = rw_lst[1]
				'''converting all iGH family genes as IGH'''
				if gn1.startswith("IGH"):
					gn1="IGH"
				if gn1.startswith("IGK"):
					gn1="IGK"
				if gn1.startswith("IGL"):
					gn1="IGL"
				if gn2.startswith("IGH"):
					gn2="IGH"
				if gn2.startswith("IGK"):
					gn2="IGK"
				if gn2.startswith("IGL"):
					gn2="IGL"
				if gn2 == "BCL2" or gn2 == "MYC"  or gn2 == "BCL6" :
					tmpgn = gn1
					gn1 = gn2
					gn2 = tmpgn
				dict_fusion_genes[samp+'__'+gn1+'__'+gn2]=1
				dict_fusion_genes[samp+'__'+gn2+'__'+gn1]=1
			fobj.close()	
		'''Reading input file and writing fusion information'''
		fobj = open(temp_dir+'/tmpfile')
		myfile = open(out_file, mode='wt')
		header = fobj.readline()
		header = header.strip()
		myfile.write(header+"\tFusion_Present_Gene\n")
		for line in fobj:
			line = line.strip()
			arr = line.split("\t")
			linenum = linenum+1
			gn1=arr[gene_1_col_num]
			gn2=arr[gene_2_col_num]
			if gn1.startswith("IGH"):
				gn1="IGH"
			if gn1.startswith("IGK"):
				gn1="IGK"
			if gn1.startswith("IGL"):
				gn1="IGL"
			if gn2.startswith("IGH"):
				gn2="IGH"
			if gn2.startswith("IGK"):
				gn2="IGK"
			if gn2.startswith("IGL"):
				gn2="IGL"
			if arr[1]+'__'+gn1+'__'+gn2 in dict_fusion_genes:
				fusion_gene="1"
			else:
				fusion_gene="0"
			myfile.write(line+"\t"+fusion_gene+"\n")
		fobj.close()
		myfile.close()
		del dict_fusion_genes
		os.remove(temp_dir+'/tmpfile')
		input_file = out_file
	

		
	'''trimmed version with only BPs with BCL2,BCL6,MYC'''
	if gtf_file != None and chr_names != None :
		fobj = open(out_file)
		myfile = open(out_file+'.bcl2_bcl6_myc.xls', mode='wt')
		header = fobj.readline()
		header = header.strip()
		myfile.write(header+"\n")
		for i in fobj:
			i = i.strip()
			arr = i.split("\t")
			samp=arr[1]
			
			gn1=arr[gene_1_col_num]
			gn2=arr[gene_2_col_num]
			if level1_bp != None :
				flag=arr[level1_col]
				if flag == "1" or gn1 in gene_list:
					myfile.write(i+"\n")
			else:
				if gn1 in gene_list:
					myfile.write(i+"\n")				
		fobj.close()
		myfile.close()

if __name__ == "__main__":
	main()