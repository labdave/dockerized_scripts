#!/usr/bin/env Rscript

library("stringr")
library(ggplot2)
library("GenVisR")

args = commandArgs(trailingOnly=TRUE)
output_plot = args[1]
output_data = args[2]

# read cnv segment files
files <- Sys.glob("/data/*cns")
print(files)
# files <- Sys.glob("*cns")

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

log_fn <- function(x){
  # read data and set column names
  df <- read.csv(x, header=TRUE, sep="\t")
  # rename
  names(df)[names(df) == "log2"] <- "segmean"
  # add sample name
  df$sample <- x
  # drop unwanted columns
  drops <- c("gene", "depth probes", "weight", "ci_lo", "ci_hi")
  df <- df[, !(names(df) %in% drops)]
  
  # return the data
  return(df)
}

# run the anonymous function defined above
cnData <- lapply(files, log_fn)

# turn the list of data frames into a single data frame
cnData <- do.call("rbind", cnData)

# remove chr prefix from chromosomes
cnData <- data.frame(lapply(cnData, function(x) {gsub("chr", "", x)}))

# convert from factor to character
cnData <- data.frame(lapply(cnData, as.character), stringsAsFactors=FALSE)

# set variables, plot graph, return to variable
plot_title = paste("Landscape of Copy Number Variation, N=", as.character(length(files)), sep="")
CN_low_cutoff =  -0.3
CN_high_cutoff = 0.3
p1 <- cnFreq(cnData, genome="hg38", CN_low_cutoff = CN_low_cutoff, CN_high_cutoff = CN_high_cutoff, plot_title = plot_title)
# remove X and Y from the figure
p1[["data"]] <- p1[["data"]][(p1[["data"]]$chromosome != "chrX" & p1[["data"]]$chromosome != "chrY") & p1[["data"]]$chromosome != "chrM", ]
# plot modified graph
ggsave(output_plot, plot=p1, device="png", width=90, height=20, units="cm")
save(p1, cnData, file=output_data)