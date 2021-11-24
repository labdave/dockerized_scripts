# Goal: merge DELLY and LUMPY non-translocation calls together
# Rachel Kositsky
# Created: 2021-11-19

import argparse
import pandas as pd
import sys
import vcf


def create_vcf_readers(delly_file, lumpy_file):
	"""Create PyVCF readers and return sample ID"""

	# Set up DELLY VCF reader
	delly_vcf_reader = vcf.Reader(open(delly_file, "r"))
	if len(delly_vcf_reader.samples) != 1: 
		raise Exception("DELLY VCF is multi-sample")
	sample_id = delly_vcf_reader.samples[0]

	# Set up LUMPY VCF reader
	lumpy_vcf_reader = vcf.Reader(open(lumpy_file, "r"))
	if len(lumpy_vcf_reader.samples) != 1:
		raise Exception("LUMPY VCF is multi-sample")
	if lumpy_vcf_reader.samples[0] != sample_id:
		raise Exception("LUMPY sample ID does not match DELLY sample ID")

	return delly_vcf_reader, lumpy_vcf_reader, sample_id


def parse_args(args=None):
	"""Parse command line arguments and return constants"""
	parser = argparse.ArgumentParser(
		description="Merge DELLY and LUMPY non-translocation calls together" \
		"into a matrix with one row per (shared or unique) call.")

	parser.add_argument("delly_input", help="DELLY VCF input")

	parser.add_argument("lumpy_input", help="LUMPY VCF input")

	parser.add_argument("output_file",
		help="Tab-delimited output table of merged structural variants.")

	args = parser.parse_args(args)

	return args


def process_record(r, source_vcf):
	"""Return a dictionary with values filled in for the record columns given
	Arguments: r, PyVCF Record from DELLY or LUMPY VCF
	Output: dict of {column_name: value} pairs for final merged file"""
	d = {}
	d["CHROM"] = r.CHROM
	d["POS"] = r.POS
	d["REF"] = r.REF

	# Fill in site information from INFO
	for k in r.INFO.keys():
		d[f"{k}_info_{source_vcf}"] = r.INFO[k]
	
	# Fill in genotype information following keys in FORMAT
	call = r.samples[0]
	for k in r.FORMAT.split(":"):
		d[f"{k}_fmt_{source_vcf}"] = call[k]

	return d


def main(args):
	"""Merge DELLY and LUMPY non-translocation calls together into a matrix with
	one row per (shared or unique) call."""

	delly_vcf_reader, lumpy_vcf_reader, sample_id = create_vcf_readers(
		args.delly_input, args.lumpy_input)

	# Set up output pandas dataframe
	fixed_cols = ["sample_id", "CHROM", "POS", "REF"]
	delly_info_cols = [f"{i}_info_delly" for i in delly_vcf_reader.infos.keys()]
	delly_format_cols = [f"{i}_fmt_delly" for i in delly_vcf_reader.formats.keys()]
	lumpy_info_cols = [f"{i}_info_lumpy" for i in lumpy_vcf_reader.infos.keys()]
	lumpy_format_cols = [f"{i}_fmt_lumpy" for i in lumpy_vcf_reader.formats.keys()]
	all_cols = fixed_cols + delly_info_cols + delly_format_cols + lumpy_info_cols + lumpy_format_cols

	df = pd.DataFrame(columns = all_cols)

	# Add all DELLY variants as-is
	delly_idx = 0
	for record in delly_vcf_reader:
		# Fill in sample ID
		df.at[delly_idx, "sample_id"] = sample_id

		record_dict = process_record(record, "delly")
		for k,v in record_dict.items():
			df.at[delly_idx, k] = v
		
		delly_idx += 1

	print(f"Parsed {delly_idx} records from DELLY")

	# Construct list for easy lookup: CHROM-POS-SVTYPE
	# then just call FIND on the list of DELLY SVs when reading through LUMPY SVs
	chrom_pos_svtype_delly = ['-'.join(i) for i in zip(df["CHROM"], df["POS"].map(str),
		df["SVTYPE_info_delly"])]

	# Merge in LUMPY variants and add new rows if needed
	next_idx = delly_idx
	n_lumpy = 0
	for record in lumpy_vcf_reader:
		# Search df to see if you've already seen this variant already
		# If so, then write it to that index
		chrom_pos_svtype = f"{record.CHROM}-{record.POS}-{record.INFO['SVTYPE']}"
		if chrom_pos_svtype in chrom_pos_svtype_delly:
			write_idx = chrom_pos_svtype_delly.index(chrom_pos_svtype)
		else:
			write_idx = next_idx
			next_idx += 1
		
		# Fill in sample ID
		df.at[write_idx, "sample_id"] = sample_id

		record_dict = process_record(record, "lumpy")
		for k,v in record_dict.items():
			df.at[write_idx, k] = v
		
		n_lumpy += 1

	n_new = next_idx - delly_idx
	print(f"Parsed {n_lumpy} records from LUMPY: {n_lumpy-n_new} in common with DELLY and {n_new} new")

	# Write output
	df.to_csv(path_or_buf = args.output_file, sep = "\t", index = False)


if __name__ == '__main__':
	main(parse_args(sys.argv[1:]))

