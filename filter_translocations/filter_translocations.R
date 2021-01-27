# R script to do filtering of translocations in SV_CALLING pipeline
# Rachel Kositsky
# 2021-01-27

library("optparse")
library("tidyverse")

# Make command line options and parse arguments ####
options_list <- list(
  make_option(c("-i", "--input_file"),
              type = "character", 
              help = "SV_CALLING tab-separated 'VCF' that's been annotated with a PON"),
  make_option(c("-t", "--filtered_translocation_output"),
              type = "character", 
              help = "filtered SV_CALLING tab-separated 'VCF'"),
  make_option(c("-g", "--gene_calls_output"),
              type = "character", 
              help = "For MYC, BCL2, and BCL6, report best-supported filtered translocations for each; TSV"),
  make_option(c("-r", "--min_reads"), default = 5,
              type = "integer", 
              help = "minimum number of reads for a translocation to be kept")
  )

opt_parser <- OptionParser(option_list = options_list)
opt <- parse_args(opt_parser)

print(opt)

# Check mandatory arguments
if(is.null(opt$input_file)) stop("input_file not set")
if(is.null(opt$filtered_translocation_output)) stop("filtered_translocation_output not set")
if(is.null(opt$gene_calls_output)) stop("gene_calls_output not set")

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

trl_df <- parse_translocation_tsv(opt$input_file)
problems(trl_df)
stop_for_problems(trl_df)
print(sprintf("Parsed %s", opt$input_file))

#### Filter translocations ####

print(sprintf("Unfiltered calls: %d", nrow(trl_df)))

trl_df <- trl_df %>% filter(!is.na(FISH_capture))
print(sprintf("After filtering for FISH capture regions (MYC, BCL2, BCL6): %d", nrow(trl_df)))

trl_df <- trl_df %>% filter(is.na(PON))
print(sprintf("After filtering out PON: %d", nrow(trl_df)))

trl_df <- trl_df %>% filter(is.na(matching_repeats))
print(sprintf("After filtering out matching repeats: %d", nrow(trl_df)))

trl_df <- trl_df %>% filter(is.na(BP1_polynt_200bp)) %>% filter(is.na(BP2_polynt_200bp)) 
print(sprintf("After filtering out polynucleotide repeats: %d", nrow(trl_df)))

trl_df <- trl_df %>% filter(is.na(segdup_100k))
print(sprintf("After filtering out segmental duplications within 100kb: %d", nrow(trl_df)))

trl_df <- trl_df %>% filter(pe_sr >= opt$min_reads)
print(sprintf("After filtering for %d+ reads: %d", opt$min_reads, nrow(trl_df)))

# not included: gene-specific filters for CCDC26, CASC8
# filter(!(BP1_Gene %in% c("CCDC26", "CASC8"))) %>% filter(!(BP2_Gene %in% c("CCDC26", "CASC8")))

# Write output
write_tsv(trl_df, file = opt$filtered_translocation_output)

#### Make FISH calls ####

# Output: the best-supported translocation for each gene in FISH_capture

# Initialize tibble with the same columns
gene_calls <- trl_df[0,]

for (gene in sort(unique(trl_df$FISH_capture))) {
  gene_trls <- trl_df %>% filter(FISH_capture == gene) %>% arrange(desc(pe_sr))
  gene_calls <- rbind(gene_calls, gene_trls[1,])
}

print(sprintf("Found translocations in %d genes: %s", nrow(gene_calls), 
              paste(sort(unique(trl_df$FISH_capture)), collapse = ", ")))

write_tsv(gene_calls, file = opt$gene_calls_output)

