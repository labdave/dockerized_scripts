import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Filter rows/columns as needed")
parser.add_argument("-i", "--input", help="input vcf")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t')

new_str = "SAMPLE\tCHROM\tPOS\n"

for index, row in df.iterrows():
	chrom = row["CHROM"]
	pos = row["POS"]
	new_str += f"Sample\t{chrom}\t{pos}\n"

with open(args.output, "w") as f:
	f.write(new_str)