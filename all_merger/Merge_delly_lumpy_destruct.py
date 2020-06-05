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
	
if __name__ == "__main__":
	main()
