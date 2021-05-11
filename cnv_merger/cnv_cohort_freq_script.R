#!/usr/bin/env Rscript

install.packages("stringr")
library("stringr")
install.packages("stringr")
library(ggplot2)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenVisR")
library("GenVisR")

args = commandArgs(trailingOnly=TRUE)
output_plot = args[1]
output_data = args[2]

# read normalized files
files_global <- Sys.glob("/data/*global*.norm.seg")
files_focal <- Sys.glob("/data/*focal*.norm.seg")

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

focal_blacklist <- function(df){
  # remove centromere-proximal no-coverage segments
  df <- df[!(df$chromosome == "1" & as.numeric(as.character(df$start)) < 123400000 & as.numeric(as.character(df$end)) > 123400000), ]
  df <- df[!(df$chromosome == "1" & as.numeric(as.character(df$start)) < 130000000 & as.numeric(as.character(df$end)) > 130000000), ]
  df <- df[!(df$chromosome == "10" & as.numeric(as.character(df$start)) < 39800000 & as.numeric(as.character(df$end)) > 39800000), ]
  df <- df[!(df$chromosome == "11" & as.numeric(as.character(df$start)) < 53400000 & as.numeric(as.character(df$end)) > 53400000), ]
  df <- df[!(df$chromosome == "12" & as.numeric(as.character(df$start)) < 35500000 & as.numeric(as.character(df$end)) > 35500000), ]
  df <- df[!(df$chromosome == "13" & as.numeric(as.character(df$start)) < 17700000 & as.numeric(as.character(df$end)) > 17700000), ]
  df <- df[!(df$chromosome == "14" & as.numeric(as.character(df$start)) < 17200000 & as.numeric(as.character(df$end)) > 17200000), ]
  df <- df[!(df$chromosome == "15" & as.numeric(as.character(df$start)) < 19000000 & as.numeric(as.character(df$end)) > 19000000), ]
  df <- df[!(df$chromosome == "16" & as.numeric(as.character(df$start)) < 36800000 & as.numeric(as.character(df$end)) > 36800000), ]
  df <- df[!(df$chromosome == "17" & as.numeric(as.character(df$start)) < 25100000 & as.numeric(as.character(df$end)) > 25100000), ]
  df <- df[!(df$chromosome == "18" & as.numeric(as.character(df$start)) < 18500000 & as.numeric(as.character(df$end)) > 18500000), ]
  df <- df[!(df$chromosome == "19" & as.numeric(as.character(df$start)) < 26200000 & as.numeric(as.character(df$end)) > 26200000), ]
  df <- df[!(df$chromosome == "2" & as.numeric(as.character(df$start)) < 91800000 & as.numeric(as.character(df$end)) > 91800000), ]
  df <- df[!(df$chromosome == "20" & as.numeric(as.character(df$start)) < 28100000 & as.numeric(as.character(df$end)) > 28100000), ]
  df <- df[!(df$chromosome == "21" & as.numeric(as.character(df$start)) < 12000000 & as.numeric(as.character(df$end)) > 12000000), ]
  df <- df[!(df$chromosome == "22" & as.numeric(as.character(df$start)) < 15000000 & as.numeric(as.character(df$end)) > 15000000), ]
  df <- df[!(df$chromosome == "3" & as.numeric(as.character(df$start)) < 90900000 & as.numeric(as.character(df$end)) > 90900000), ]
  df <- df[!(df$chromosome == "4" & as.numeric(as.character(df$start)) < 50000000 & as.numeric(as.character(df$end)) > 50000000), ]
  df <- df[!(df$chromosome == "5" & as.numeric(as.character(df$start)) < 48800000 & as.numeric(as.character(df$end)) > 48800000), ]
  df <- df[!(df$chromosome == "6" & as.numeric(as.character(df$start)) < 59800000 & as.numeric(as.character(df$end)) > 59800000), ]
  df <- df[!(df$chromosome == "7" & as.numeric(as.character(df$start)) < 60100000 & as.numeric(as.character(df$end)) > 60100000), ]
  df <- df[!(df$chromosome == "8" & as.numeric(as.character(df$start)) < 45200000 & as.numeric(as.character(df$end)) > 45200000), ]
  df <- df[!(df$chromosome == "9" & as.numeric(as.character(df$start)) < 45500000 & as.numeric(as.character(df$end)) > 45500000), ]
  # remove igh
  df <- df[!(df$chromosome == "14" & as.numeric(as.character(df$start)) > 105500000), ]
  
  # return modified df
  return(df) 
}

log_fn <- function(x){
  # read data and set column names
  df <- read.csv(x, header=TRUE, sep='\t')
  names(df)[names(df) == 'Chromosome'] <- 'chromosome'
  names(df)[names(df) == 'Start'] <- 'start'
  names(df)[names(df) == 'End'] <- 'end'
  names(df)[names(df) == 'MEAN_LOG2_COPY_RATIO'] <- 'segmean'
  # remove unwanted chromosomes
  df <- df[!grepl("KI", df$chromosome),]
  df <- df[!grepl("GL", df$chromosome),]
  df <- df[!grepl("X", df$chromosome),]
  df <- df[!grepl("Y", df$chromosome),]
  # add sample name
  df$sample <- x
  # drop unwanted columns
  drops <- c("Call")
  df <- df[, !(names(df) %in% drops)]
  
  # return the data
  return(df)
}

# run the anonymous function defined above
cnData_global <- lapply(files_global, log_fn)
cnData_focal <- lapply(files_focal, log_fn)

# turn the list of data frames into a single data frame
cnData_global <- do.call("rbind", cnData_global)
cnData_focal <- do.call("rbind", cnData_focal)

# remove chr prefix from chromosomes
cnData_global <- data.frame(lapply(cnData_global, function(x) {gsub("chr", "", x)}))
cnData_focal <- data.frame(lapply(cnData_focal, function(x) {gsub("chr", "", x)}))

# blacklist centromere and igh regions
cnData_focal <- focal_blacklist(cnData_focal)

# convert from factor to character
cnData_global <- data.frame(lapply(cnData_global, as.character), stringsAsFactors=FALSE)
cnData_focal <- data.frame(lapply(cnData_focal, as.character), stringsAsFactors=FALSE)

# bring both to the same scale
cnData_global$segmean[as.numeric(as.character(cnData_global$segmean)) < log2(0.75)] <- -2
cnData_global$segmean[as.numeric(as.character(cnData_global$segmean)) > log2(1.25)] <- 2

# merge into one
cnData <- rbind(cnData_focal, cnData_global)

# set variables, plot graph, return to variable
plot_title = 'DLBCL Mini Project (high coverage >10x) N=38, Exome Panel'
plot_title = 'Normals, N=25, Noncoding Panel'
CN_low_cutoff =  -1
CN_high_cutoff = 1
p1 <- cnFreq(cnData, genome="hg38", CN_low_cutoff = CN_low_cutoff, CN_high_cutoff = CN_high_cutoff, plot_title = plot_title)
# remove X and Y from the figure
p1[["data"]] <- p1[["data"]][(p1[["data"]]$chromosome != "chrX" & p1[["data"]]$chromosome != "chrY") & p1[["data"]]$chromosome != "chrM", ]
# plot modified graph
ggsave(output_plot, plot=p1, device="png", width=40, height=10, units="cm")
