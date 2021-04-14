# R script to merge translocations that came from SV_CALLING
# Rachel Kositsky
# 2021-04-14

library("optparse")
library("tidyverse")

# Make command line options and parse arguments ####
options_list <- list(
  make_option(c("-i", "--input_files"), dest = "input_files",
              type = "character", help = "comma-separated list of input files from SV_CALLING pipeline"),
  make_option(c("-o", "--output"), dest = "output", 
              type = "character", help = "path to merged output file (TSV)"),
  make_option(c("-d", "--merge_distance"), dest = "merge_distance", 
              type = "integer", help = "distance in base pairs for merging adjacent translocation calls"))

opt_parser <- OptionParser(option_list = options_list)
opt <- parse_args(opt_parser)

print(opt)

input_files_list = unlist(strsplit(opt$input_file, ","))
merge_distance = opt$merge_distance
output = opt$output

#### Download all translocation positions ####

parse_translocation_tsv <- function(path) {
  return(read_tsv(path,
                  col_types = 
                    cols(
                      dave_lab_id = col_character(),
                      chr1 = col_character(),
                      pos1 = col_double(),
                      chr2 = col_character(),
                      pos2 = col_double(),
                      pe = col_double(),
                      sr = col_double(),
                      pe_sr = col_double(),
                      Num_callers = col_double(),
                      Callers = col_character(),
                      discowave_chrom = col_character(),
                      discowave_start = col_double(),
                      discowave_stop = col_double(),
                      discowave_pe = col_double(),
                      discowave_depth = col_double(),
                      discowave_pct_to_partner = col_double(),
                      discowave_partner_chrom = col_character(),
                      discowave_evenness = col_double(),
                      BP1_gene = col_character(),
                      BP2_gene = col_character(),
                      BP1_repeats_200bp = col_character(),
                      BP2_repeats_200bp = col_character(),
                      matching_repeats = col_character(),
                      BP1_polynt_200bp = col_character(),
                      BP2_polynt_200bp = col_character(),
                      BP1_segdup_200bp = col_character(),
                      BP2_segdup_200bp = col_character(),
                      segdup_100k = col_character(),
                      segdup_1M = col_character(),
                      segdup_10M = col_character(),
                      segdup_100M = col_character(),
                      Ig = col_character(),
                      FISH_capture = col_character(),
                      PON = col_character())
  ))
}

trl_df <- parse_translocation_tsv(input_files_list[1])
stop_for_problems(trl_df)
n_files <- length(input_files_list)
print(sprintf("Parsed %s [1/%d]", input_files_list[1], n_files))

if (n_files > 1) {
  for (i in 2:n_files) {
    trl_df <- rbind(trl_df, parse_translocation_tsv(input_files_list[i]))
    stop_for_problems(trl_df)
    print(sprintf("Parsed %s [%d/%d]", input_files_list[i], i, n_files))
  }
}

#### Make the dataframe for the merged list ####
# Pare down the columns we use: just location, sample ID, and pe_sr
# TODO: try keeping all columns
df <- trl_df %>% rename("sample_id" = "dave_lab_id")
#df <- trl_df %>% rename("sample_id" = "dave_lab_id") %>% 
#  select(c("sample_id", "chr1", "pos1", "chr2", "pos2", "pe_sr"))


# To account for both ways of ordering a translocation's breakpoints, flip the
# first and second positions
flipped_df <- df
flipped_df[,c("chr1", "pos1", "chr2", "pos2")] <- df[,c("chr2", "pos2", "chr1", "pos1")]
df <- bind_rows(df, flipped_df)
rm(flipped_df)

# Sort dataframe by chr1, pos1, chr2, pos2, sample_id, in that order
chrom_order <- map_chr(c(1:22, "X", "Y", "M"), ~ paste0("chr", .x))
df <- df %>% arrange(sample_id) %>% 
  arrange(pos2) %>% arrange(match(chr2, chrom_order)) %>% 
  arrange(pos1) %>% arrange(match(chr1, chrom_order))

print("Created normal dataframe; now collapsing windows...")

# Collapse windows ####
# Copy over first window
output_df <- df[1,]

# Logging to retain sanity in long runs: every 10%
log_points = ceiling(1:10 * nrow(df) / 10)

# We assume that the input windows are sorted by sample and chromosomal position.
for(i in 2:nrow(df)) {
  j <- nrow(output_df)
  
  # Print log message every 10%
  if (i %in% log_points) {
    print(sprintf("%d%% done.. [%d/%d rows]", 
                  which(log_points == i)*10, i, nrow(df)))
  } 
  
  if ((df[i, "chr1"] == output_df[j, "chr1"]) &&
      ((df[i, "pos1"] - output_df[j, "pos1"]) <= merge_distance) &&
      (df[i, "chr2"] == output_df[j, "chr2"]) &&
      ((df[i, "pos2"] - output_df[j, "pos2"]) <= merge_distance)) {
    # Condition 1: merge current row into previous one
    
    # 1) Update the samples that have this so that "sample_id" is a
    # comma-separated list of samples with this breakpoint
    df_samp_id = df[i, "sample_id"]
    output_samp_ids = output_df[j, "sample_id"]
    if (!grepl(df_samp_id, output_samp_ids)) {
      output_df[j, "sample_id"] = paste(output_samp_ids, df_samp_id, sep = ",")
    }
    
    # 2) Update fields if the new window has more supporting reads than current
    # breakpoint. This means the output has stats from the best-supported breakpoint.
    if (df[i, "pe_sr"] > output_df[j, "pe_sr"]) {
      # Replace everything after sample_id column
      output_df[j, 2:ncol(output_df)] <- df[i, 2:ncol(output_df)]
    }
  } else {
    # Condition 2: add a new row without merging
    output_df[j+1,] <- df[i,]
  }
}

# Add n_samples column for each windows
output_df$n_samples <- str_count(output_df$sample_id, pattern = ",") + 1

# TODO: reorder columns

# Save output
write_delim(output_df, file = output, delim = "\t")
print(sprintf("Wrote output at %s", output))


print("Done!")
