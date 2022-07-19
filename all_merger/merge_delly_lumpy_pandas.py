# Rachel Kositsky
# 2022-06-24

# Goal: merge outputs of dellymerger and lumpymerger into one output table

import argparse
import pandas as pd
import pybedtools
import os

CHROM_ORDERS = ["chr" + str(i) for i in list(range(1,23)) + ["X", "Y", "M"]]

def get_sorted_translocation_key(chr1, pos1, chr2, pos2):
	"""
	Order chr1 < chr2
	Given two breakpoints, return a string with the chromosomes sorted
	e.g. for input chr18,123000,chr14,900400, return 'chr14:900400-chr18:123000'
	"""
	# Check inputs for right format
	assert(chr1[:3] == "chr")
	assert(chr2[:3] == "chr")
	assert(isinstance(pos1, int))
	assert(isinstance(pos2, int))

	# First sort by position in case of the same chromosome
	sorted_keys = sorted([(chr1, pos1), (chr2, pos2)], key = lambda x: x[1])
	# Then do the real sort by chromosome
	sorted_keys = sorted(sorted_keys, key = lambda x: CHROM_ORDERS.index(x[0]))

	# Construct translocation key
	return f"{sorted_keys[0][0]}:{sorted_keys[0][1]}-{sorted_keys[1][0]}:{sorted_keys[1][1]}"


def initialize_files(delly_file, lumpy_file):
	"""Initialize data frames
	1. Reads in DELLY merger output and outputs pandas dataframe ready for merge
	2. Reads in LUMPY merger output and outputs pandas dataframe ready for merge
	3. Initialize output dataframe with appropriate columns
	"""
	# Read inputs into pandas dataframes
	delly_df = pd.read_csv(delly_file, sep = "\t")
	lumpy_df = pd.read_csv(lumpy_file, sep = "\t")
	print("-------------- initialize_files 2 ----------------------------")

	# Initialize output dataframe
	output_columns = ["dave_lab_id", "chr1", "pos1", "chr2", "pos2",
		"pe", "sr", "pe_sr", "caller"] + list(delly_df.columns[5:]) + \
		list(lumpy_df.columns[5:]) + ["Callers", "Num_callers"]
	out_df = pd.DataFrame(columns = output_columns)

	print("-------------- initialize_files 2 ----------------------------")

	# Add translocation key column for merging
	delly_df.loc[:, "trl_key"] = delly_df.apply(lambda row: 
		get_sorted_translocation_key(row["Delly_CHR1"], row["Delly_POS1"],
			row["Delly_CHR2"], row["Delly_POS2"]), axis = 1)
	print("-------------- initialize_files 3 ----------------------------")

	lumpy_df.loc[:, "trl_key"] = lumpy_df.apply(lambda row: 
		get_sorted_translocation_key(row["Lumpy_CHROM1"], row["Lumpy_POS1"],
			row["Lumpy_CHROM2"], row["Lumpy_POS2"]), axis = 1)

	return delly_df, lumpy_df, out_df
	

def make_bed(sv_df, caller, distance=100):
	"""Make a BED file with the sorted key for later merging

	Arguments:
		sv_df: DataFrame output from initialize_files with all fields in 
			dellymerger or lumpymerger output files
		caller: string, "DELLY" or "LUMPY"
		distance: how close in basepairs both breakpoints need to be in order to
			be considered the same translocation.
	Returns:
		pybedtools BedTool object with each line being a breakpoint and the
		4th column being the full sorted trl_key
	"""

	print("-------------- make_bed 0 ----------------------------")
	return_list = [None, None]

	if caller == "DELLY":
		colnames_dict = {"chr1": "Delly_CHR1", "chr2": "Delly_CHR2", 
			"pos1": "Delly_POS1", "pos2": "Delly_POS2"}
	elif caller == "LUMPY":
		colnames_dict = {"chr1": "Lumpy_CHROM1", "chr2": "Lumpy_CHROM2",
			"pos1": "Lumpy_POS1", "pos2": "Lumpy_POS2"}
	else:
		raise Exception(f"Unknown input for 'caller' in make_bed: {caller}")

	out_df = pd.DataFrame(columns = ["chrom", "start", "end", "trl_key"])

	print("-------------- make_bed 1 ----------------------------")
	# Copy both breakpoints into output bedtools DataFrame
	for idx in ["1", "2"]:
		bp_df = sv_df.loc[:, [colnames_dict["chr" + idx], 
			colnames_dict["pos" + idx], "trl_key"]]
		bp_df.loc[:, "chrom"] = bp_df.loc[:, colnames_dict["chr" + idx]]
		bp_df.loc[:, "start"] = bp_df.loc[:, colnames_dict["pos" + idx]].map(lambda x: int(x - distance/2))
		bp_df.loc[:, "end"] = bp_df.loc[:, colnames_dict["pos" + idx]].map(lambda x: int(x + distance/2))

		# Copy to out_df
		out_df = pd.concat([out_df, bp_df.loc[:, ["chrom", "start", "end", "trl_key"]]])

	print("-------------- make_bed 2 ----------------------------")
	# Save to temporary BED file for merging
	file_name = f"tmp.{caller}.{distance}bp.bed"
	out_df.to_csv(file_name, sep="\t", header=False, index=False)

	# TODO: improve runtime by sorting this input file

	# Return the merged BED file as a pybedtools BedTool object
	return pybedtools.BedTool(file_name)


def check_both_breakpoints(row, distance = 100):
	"""Checks the DELLY and LUMPY translocation keys to make sure they align on
	both breakpoints. Return True if the align within the distance, False otherwise.
	Assumes that the chromosomes are in the same order in DELLY and LUMPY
	because they've been sorted before."""

	# Extract chromosomes and positions
	delly_trl_key = row["delly_trl_key"]
	lumpy_trl_key = row["lumpy_trl_key"]

	delly_chr1, delly_pos1 = delly_trl_key.split("-")[0].split(":")
	delly_chr2, delly_pos2 = delly_trl_key.split("-")[1].split(":")

	lumpy_chr1, lumpy_pos1 = lumpy_trl_key.split("-")[0].split(":")
	lumpy_chr2, lumpy_pos2 = lumpy_trl_key.split("-")[1].split(":")

	if (delly_chr1 != lumpy_chr1) or (delly_chr2 != lumpy_chr2):
		return False

	if abs(int(delly_pos1) - int(lumpy_pos1)) > distance:
		return False

	if abs(int(delly_pos2) - int(lumpy_pos2)) > distance:
		return False

	return True


def intersect_trl_beds(delly_bed, lumpy_bed, distance = 100):
	"""Given the BED files with positions centered around each translocation 
	breakpoint, return translocation positions called by both DELLY and LUMPY.
	Arguments:
		delly_bed: pybedtools object from make_bed, dataframe in BED4 style with trl_key as 4th column
	 	lumpy_bed: pybedtools object from make_bed, dataframe in BED4 style with trl_key as 4th column

	Outputs:
		shared transloctions as a DataFrame with delly and lumpy paired translocation keys
	"""
	
	print("-------------- intersect_trl_beds 0 ----------------------------")
	# Intersect DELLY and LUMPY BED files and save to a temporary output file
	intersect_file_name = "tmp.intersected.bed"
	delly_bed.intersect(lumpy_bed, wa=True, wb=True).moveto(intersect_file_name)

	# Read in output BED file into dataframe
	intersect_df = pd.read_csv(intersect_file_name, sep = "\t", 
		names = ["delly_chrom", "delly_start", "delly_end", "delly_trl_key",
		"lumpy_chrom", "lumpy_start", "lumpy_end", "lumpy_trl_key"])

	print("-------------- intersect_trl_beds 1 ----------------------------")
	# Go through one-way intersected translocations and check for two-way compatibility
	intersect_df["keep"] = intersect_df.apply(check_both_breakpoints, axis = 1)


	print("-------------- intersect_trl_beds 2 ----------------------------")
	# Filter for rows that are kept
	intersect_df = intersect_df.loc[intersect_df["keep"], :]

	print("-------------- intersect_trl_beds 3 ----------------------------")
	# Get DELLY and LUMPY translocation keys that passed two-way inspection
	shared_keys = intersect_df.loc[:, ["delly_trl_key", "lumpy_trl_key"]].copy()
	
	return shared_keys


def prepare_single_caller_dataframe(single_caller_df, caller):
	"""Does the nitty-gritty of cleaning up columns in delly_df or lumpy_df
	before merging in the single-caller translocations

	Arguments:
		single_caller_df: pandas DataFrame from merge_delly_lumpy_translocations
		caller: string, "DELLY" or "LUMPY"
	"""
	single_caller_df.loc[:, "chr1"] = single_caller_df.loc[:, "trl_key"].map(lambda x: x.split("-")[0].split(":")[0])
	single_caller_df.loc[:, "pos1"] = single_caller_df.loc[:, "trl_key"].map(lambda x: x.split("-")[0].split(":")[1])
	single_caller_df.loc[:, "chr2"] = single_caller_df.loc[:, "trl_key"].map(lambda x: x.split("-")[1].split(":")[0])
	single_caller_df.loc[:, "pos2"] = single_caller_df.loc[:, "trl_key"].map(lambda x: x.split("-")[1].split(":")[1])

	if caller == "DELLY":
		single_caller_df.loc[:, ["pe", "sr"]] = single_caller_df.loc[:, ["Delly_PE_NReads", "Delly_SR_NReads"]]
		single_caller_df.loc[:, "caller"] = "DELLY"
		single_caller_df.loc[:, "Callers"] = "DELLY"
	elif caller == "LUMPY":
		single_caller_df.loc[:, ["pe", "sr"]] = single_caller_df.loc[:, ["Lumpy_PE", "Lumpy_SR"]]
		single_caller_df.loc[:, "caller"] = "LUMPY"
		single_caller_df.loc[:, "Callers"] = "LUMPY"
	else:
		raise Exception(f"Unknown input for 'caller', expect 'DELLY' or 'LUMPY': {caller}")

	single_caller_df.loc[:, "pe_sr"] = single_caller_df.apply(lambda row: row["pe"] + row["sr"], axis = 1)
	single_caller_df.loc[:, "Num_callers"] = 1

	# Drop temporary columns before copying to output
	single_caller_df = single_caller_df.drop(columns = ["trl_key", "shared"])

	return single_caller_df


def prepare_shared_caller_dataframe(shared_caller_df, caller):
	"""Does the nitty-gritty of cleaning up columns in delly_df or lumpy_df
	before merging is done in get_merged_row().

	Arguments:
		shared_caller_df: pandas DataFrame from merge_delly_lumpy_translocations
		caller: string, "DELLY" or "LUMPY"
	"""

	print("-------------- prepare_shared_caller_dataframe 0 ----------------------------")
	trl_key_df = shared_caller_df["trl_key"].copy()
	print("-------------- prepare_shared_caller_dataframe 0.5 ----------------------------")
	shared_caller_df["chr1"] = trl_key_df.map(lambda x: x.split("-")[0].split(":")[0])
	print("-------------- prepare_shared_caller_dataframe 0.6 ----------------------------")
	shared_caller_df.loc[:, "pos1"] = trl_key_df.map(lambda x: x.split("-")[0].split(":")[1])
	print("-------------- prepare_shared_caller_dataframe 0.7 ----------------------------")
	shared_caller_df.loc[:, "chr2"] = trl_key_df.map(lambda x: x.split("-")[1].split(":")[0])
	print("-------------- prepare_shared_caller_dataframe 0.8 ----------------------------")
	shared_caller_df.loc[:, "pos2"] = trl_key_df.map(lambda x: x.split("-")[1].split(":")[1])

	print("-------------- prepare_shared_caller_dataframe 1 ----------------------------")
	if caller == "DELLY":
		shared_caller_df[["pe", "sr"]] = shared_caller_df.loc[:, ["Delly_PE_NReads", "Delly_SR_NReads"]]
		shared_caller_df.loc[:, "caller"] = "DELLY"
		print("-------------- prepare_shared_caller_dataframe 2D ----------------------------")
	elif caller == "LUMPY":
		shared_caller_df[["pe", "sr"]] = shared_caller_df.loc[:, ["Lumpy_PE", "Lumpy_SR"]]
		shared_caller_df["caller"] = "LUMPY"
		print("-------------- prepare_shared_caller_dataframe 2L ----------------------------")
	else:
		raise Exception(f"Unknown input for 'caller', expect 'DELLY' or 'LUMPY': {caller}")

	print("-------------- prepare_shared_caller_dataframe 3 ----------------------------")
	shared_caller_df["pe_sr"] = shared_caller_df.apply(lambda row: row["pe"] + row["sr"], axis = 1)
	shared_caller_df["Callers"] = "DELLY,LUMPY"
	shared_caller_df["Num_callers"] = 2
	print("-------------- prepare_shared_caller_dataframe 4 ----------------------------")

	# Drop temporary columns before copying to output
	shared_caller_df = shared_caller_df.drop(columns = ["shared"])

	print("-------------- prepare_shared_caller_dataframe 5 ----------------------------")

	# Set trl_key as the index so that you can use it for lookup later. This 
	# also conveniently drops it.
	shared_caller_df = shared_caller_df.set_index("trl_key")

	return shared_caller_df


def get_merged_row(delly_row, lumpy_row, output_columns):
	""" Merge shared translocation from DELLY and LUMPY into one line in output.
	Arguments: 
		delly_row: row from delly_df of shared translocation
		lumpy_row: row from lumpy_df of shared translocation
	Returns:
		row with columns of out_df
	"""

	# Get each caller's read support
	delly_reads = int(delly_row["pe"]) + int(delly_row["sr"])
	lumpy_reads = int(lumpy_row["pe"]) + int(lumpy_row["sr"])

	# Choose caller with top split + paired end reads and create merged row
	if (lumpy_reads > delly_reads):
		output_row = lumpy_row[["dave_lab_id"] + lumpy_row.index[5:].tolist()]

		# Fill in all other columns from DELLY
		for idx in delly_row.index:
			if idx not in output_row.index:
				output_row[idx] = delly_row[idx]
	else:
		output_row = delly_row[["dave_lab_id"] + delly_row.index[5:].tolist()]

		# Fill in all other columns from LUMPY
		for idx in lumpy_row.index:
			if idx not in output_row.index:
				output_row[idx] = lumpy_row[idx]

	# Make output row match the column order of the output dataframe
	# removes caller-specific columns, e.g. LUMPY_CHROM1, DELLY_POS2
	output_row = output_row[output_columns]

	return output_row


def merge_delly_lumpy_translocations(delly_df, lumpy_df, out_df, shared_trls):
	"""
	Given DELLY and LUMPY DataFrames, add translocations to output DataFrame.
	Copy over single-caller translocations and merge shared translocations.
	Arguments:
		delly_df: pandas DataFrame with DELLY calls, with 'shared' column
		lumpy_df: pandas DataFrame with DELLY calls, with 'shared' column
		out_df: output pandas DataFrame; is empty at first, has column names

	========== heuristic for merging calls ==========
	-> Merge calls with both breakpoints within 100bp of each other
	-> Can only merge if called from same sample and different caller
	-> The BP positions are defined by the following hierarchy
		- The one with more reads; if equal:
		- Lumpy > Delly 
	-> Collapse rows into one for merged reads
	"""

	print("-------------- merge_delly_lumpy_translocations 0 ----------------------------")
	# Get shared translocations from DELLY
	delly_shared = prepare_shared_caller_dataframe(delly_df.loc[delly_df["shared"], :], "DELLY")
	print(f"DELLY shared: {len(delly_shared.index)} out of {len(delly_df.index)}")


	print("-------------- merge_delly_lumpy_translocations 1 ----------------------------")
	# Add missing columns to delly_only before concatenation
	delly_only = delly_df.loc[delly_df["shared"] == False, :].copy()
	delly_only = prepare_single_caller_dataframe(delly_only, "DELLY")

	print("-------------- merge_delly_lumpy_translocations 2 ----------------------------")
	# Copy DELLY-only translocations to output DataFrame
	out_df = pd.concat([out_df, delly_only.loc[:, ["dave_lab_id"] + delly_only.columns[5:].tolist()]])


	print("-------------- merge_delly_lumpy_translocations 3 ----------------------------")
	# Get shared translocations from LUMPY
	lumpy_shared = prepare_shared_caller_dataframe(lumpy_df.loc[lumpy_df["shared"], :], "LUMPY")
	print(f"LUMPY shared: {len(lumpy_shared.index)} out of {len(lumpy_df.index)}")

	print("-------------- merge_delly_lumpy_translocations 4 ----------------------------")
	# Add missing columns to lumpy_only before concatenation
	lumpy_only = lumpy_df.loc[lumpy_df["shared"] == False, :].copy()
	lumpy_only = prepare_single_caller_dataframe(lumpy_only, "LUMPY")

	# Copy LUMPY-only translocations to output DataFrame
	out_df = pd.concat([out_df, lumpy_only.loc[:, ["dave_lab_id"] + lumpy_only.columns[5:].tolist()]])

	print("-------------- merge_delly_lumpy_translocations 5 ----------------------------")
	# Add shared DELLY and LUMPY calls with merging procedure
	for index, row in shared_trls.iterrows():
		x = get_merged_row(
			delly_shared.loc[shared_trls.loc[index, "delly_trl_key"], :],
			lumpy_shared.loc[shared_trls.loc[index, "lumpy_trl_key"], :],
			out_df.columns)
		out_df = pd.concat([out_df, pd.DataFrame(x).transpose()])
	
	return out_df


def main(args):
	"""
	Use pandas to read in the DELLY and LUMPY tables from dellymerger and
	lumpymerger. Then construct the output merged calls.
	"""
	
	print("-------------- main 0 ----------------------------")
	delly_df, lumpy_df, out_df = initialize_files(args.delly_file, args.lumpy_file)
	
	print("-------------- main 1 ----------------------------")
	# Get shared translocations
	shared_trls = intersect_trl_beds(make_bed(delly_df, caller="DELLY"), 
		make_bed(lumpy_df, caller = "LUMPY"))

	print("-------------- main 2 ----------------------------")
	# Drop duplicate rows
	print(f"Before dropping duplicates, shared_trls has {len(shared_trls.index)} rows")
	shared_trls = shared_trls.drop_duplicates()
	print(f"After dropping duplicates, shared_trls has {len(shared_trls.index)} rows")

	print("-------------- main 3 ----------------------------")
	# Mark which translocations are shared
	delly_df["shared"] = delly_df.loc[:, "trl_key"].map(lambda x: x in shared_trls["delly_trl_key"].tolist())
	lumpy_df["shared"] = lumpy_df.loc[:, "trl_key"].map(lambda x: x in shared_trls["lumpy_trl_key"].tolist())

	print("-------------- main 4 ----------------------------")
	# Write shared and single-caller translocations to output
	out_df = merge_delly_lumpy_translocations(delly_df, lumpy_df, out_df, shared_trls)

	# Write output file
	out_df.to_csv(args.output_file, sep = "\t", index = False)
	print("Wrote {0} calls to output file {1}".format(len(out_df.index), args.output_file))
	return


if __name__ == "__main__":
	# Read in arguments
	parser = argparse.ArgumentParser(description="DELLY and LUMPY merger")
	parser.add_argument('--output_file', required=True, 
		help = "Output for DELLY + LUMPY merged tab-separated table")
	parser.add_argument("--delly_file", required = True,
		help = "dellymerger input file")
	parser.add_argument("--lumpy_file", required = True,
		help = "lumpymerger input file")

	args = parser.parse_args()

	main(args)
