"""
Author:
Devang Thakkar

Date:
5 May 2021

Input:
This script takes in results from gatk_callcopyratiosegments__dna module in the 
CNV pipeline.

Output:
Removes parts of segments matching the blacklist

Rationale: 
The segments described the blaclist region are appearing as gains in multiple 
normal samples (9/20) and are obvious artefacts of PON denoising.

Usage:

python normalize.py <gatk_seg_file> <blacklist_seg_file>
"""

import math
import sys

in_file = sys.argv[1]
out_file = sys.argv[2]
blacklist_file = sys.argv[3]

blacklist = []
with open(blacklist_file, "r") as f:
	for line in f:
		line_arr = line.strip().split()
		blacklist.append(line_arr)

file_str = ""
with open(in_file, "r") as f:
	for line in f:

		# skip headers
		if "@" in line or "MEAN_LOG2_COPY_RATIO" in line or "Segment_Mean" in line:
			file_str += line
			continue
			
		line_arr = line.strip().split()
		line_arr1 = [i for i in line_arr]  # copy

		# check for overlap with blacklist
		for item in blacklist:

			# set flags
			edit_line = False
			skip_line = False
			copy_line = False

			# chr match
			if line_arr[1] == item[0]:

				# seg is around the start of blacklist region
				if (int(line_arr[2]) <= int(item[1])) and (int(line_arr[3]) >= int(item[1])) and (int(line_arr[3]) <= int(item[2])):
					# scale num segments assuming uniform bait distribution - works well enough
					line_arr[4] = str(int(float(line_arr[4])*(float(item[1])-float(line_arr[2]))/(float(line_arr[3])-float(line_arr[2])))+1)
					line_arr[3] = item[1]
					edit_line = True

				# seg is around stop of blacklist region
				if (int(line_arr[2]) >= int(item[1])) and (int(line_arr[2]) <= int(item[2])) and (int(line_arr[3]) >= int(item[2])):
					# scale num segments assuming uniform bait distribution - works well enough
					line_arr[4] = str(int(float(line_arr[4])*(float(line_arr[3])-float(item[2]))/(float(line_arr[3])-float(line_arr[2])))+1)
					line_arr[2] = item[2]
					edit_line = True

				# seg is inside blacklist region
				if (int(line_arr[2]) >= int(item[1])) and (int(line_arr[3]) <= int(item[2])):
					skip_line = True

				# seg is around blacklist region
				if (int(line_arr[2]) <= int(item[1])) and (int(line_arr[3]) >= int(item[2])):
					# scale num segments assuming uniform bait distribution - works well enough
					line_arr[4] = str(int(float(line_arr[4])*(float(item[1])-float(line_arr[2]))/(float(line_arr[3])-float(line_arr[2])))+1)
					line_arr[3] = item[1]
					line_arr1[4] = str(int(float(line_arr1[4])*(float(line_arr1[3])-float(item[2]))/(float(line_arr1[3])-float(line_arr1[2])))+1)
					line_arr1[2] = item[2]
					copy_line = True

			if edit_line:
				print(line_arr)
				print("edit")
				line = "\t".join(line_arr)
				break
			if skip_line:
				print(line_arr)
				print("skip")
				line = ""
				break
			if copy_line:
				print(line_arr)
				print("copy")
				line = "\t".join(line_arr) + "\n" + "\t".join(line_arr1)
				break

		if line:
			file_str += line.strip()+"\n"

# write file
with open(out_file, "w") as f:
	f.write(file_str)