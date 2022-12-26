#!/usr/bin/env Rscript
#Script to merge all MSI pipeline outputs
#Veronica Russell
#12/16/22
#Dave Lab

#get command line arguments
#input order must be 
#1) tabbed file with sample ID's and their MSI file directory outputs
#2) defined output .csv file for unstable sites
#3) defined output .csv file for all sites
#4) defined output .csv file for summary
args = commandArgs(trailingOnly=TRUE)

#get table of files
file_tbl <- read.table(args[1], header=TRUE)
head(file_tbl)
#define output file
output_file_all <- args[2]
output_file_unstable <- args[3]
output_file_summary <- args[4]

final_table_all <- data.frame(matrix(ncol=length(colnames(read.table(file_tbl[grep("all",file_tbl$V1)[1],1], header=TRUE)))+1, nrow=0))
names(final_table_all) <- c(colnames(read.table(file_tbl[grep("all",file_tbl$V1)[1],1], header=TRUE)), "sample_id")

final_table_unstable <- data.frame(matrix(ncol=length(colnames(read.table(file_tbl[grep("unstable",file_tbl$V1)[1],1], header=TRUE)))+1,nrow=0))
names(final_table_unstable) <- c(colnames(read.table(file_tbl[grep("unstable",file_tbl$V1)[1],1], header=TRUE)), "sample_id")

final_table_summary <- data.frame(matrix(ncol=length(colnames(read.table(file_tbl[!grepl("all",file_tbl$V1) & !grepl("unstable",file_tbl$V1),]$V1[1], header=TRUE, check.names = FALSE)))+1,nrow=0))
names(final_table_summary) <- c(colnames(read.table(file_tbl[!grepl("all",file_tbl$V1) & !grepl("unstable",file_tbl$V1),]$V1[1], header=TRUE, check.names = FALSE)), "sample_id")

#loop through all MSI sites files
for(i in 1:nrow(file_tbl)){
  print(i)
  #check type
  if(grepl("all",file_tbl[i,1])){
    tmp_file_all <- read.table(file_tbl[i,1], header=TRUE)
    tmp_file_all$sample_id <- file_tbl[i,2]
    final_table_all <- rbind(final_table_all, tmp_file_all)
  }
  if(grepl("unstable",file_tbl[i,1])){
    tmp_file_unstable <- read.table(file_tbl[i,1], header=TRUE)
    tmp_file_unstable$sample_id <- file_tbl[i,2]
    final_table_unstable <- rbind(final_table_unstable, tmp_file_unstable)
  }
  if(!grepl("unstable",file_tbl[i,1]) & !grepl("all",file_tbl[i,1])){
    print(i)
    tmp_file_summary <- read.table(file_tbl[i,1], header=TRUE, check.names = FALSE)
    tmp_file_summary$sample_id <- file_tbl[i,2]
    final_table_summary <- rbind(final_table_summary, tmp_file_summary)
  }
}

#write to outputs
write.table(final_table_all, output_file_all, row.names = FALSE)
write.table(final_table_unstable, output_file_unstable, row.names = FALSE)
write.table(final_table_summary, output_file_summary, row.names = FALSE)
