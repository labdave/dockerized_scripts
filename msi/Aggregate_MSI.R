#!/usr/bin/env Rscript
#Script to merge all MSI pipeline outputs
#Veronica Russell
#12/16/22
#Dave Lab

#get command line arguments
#input order must be 
#1) tabbed file with header "files"\t"samples" containing their MSI file directory outputs and sample IDs
#2) defined output .txt file for all sitesdefined output .csv file for unstable sites
#3) defined output .txt file for unstable sites
#4) defined output .txt file for summary
args = commandArgs(trailingOnly=TRUE)

#get table of files
file_tbl <- read.table(args[1], header=TRUE)
head(file_tbl)
#define output file
output_file_all <- args[2]
output_file_unstable <- args[3]
output_file_summary <- args[4]

final_table_all <- data.frame(matrix(ncol=length(colnames(read.table(file_tbl[grep("all",file_tbl$files)[1],1], header=TRUE)))+1, nrow=0))
names(final_table_all) <- c(colnames(read.table(file_tbl[grep("all",file_tbl$files)[1],1], header=TRUE)), "sample_id")
head(final_table_all)

final_table_unstable <- data.frame(matrix(ncol=length(colnames(read.table(file_tbl[grep("unstable",file_tbl$files)[1],1], header=TRUE)))+1,nrow=0))
names(final_table_unstable) <- c(colnames(read.table(file_tbl[grep("unstable",file_tbl$files)[1],1], header=TRUE)), "sample_id")
head(final_table_unstable)

final_table_summary <- data.frame(matrix(ncol=length(colnames(read.table(file_tbl[!grepl("all",file_tbl$files) & !grepl("unstable",file_tbl$files),]$files[1], header=TRUE, check.names = FALSE)))+1,nrow=0))
names(final_table_summary) <- c(colnames(read.table(file_tbl[!grepl("all",file_tbl$files) & !grepl("unstable",file_tbl$files),]$files[1], header=TRUE, check.names = FALSE)), "sample_id")
head(final_table_summary)

#loop through all MSI sites files
for(i in 1:nrow(file_tbl)){
  print(i)
  #check type
  if(grepl("all",file_tbl[i,1])){
    tmp_file_all <- read.table(file_tbl[i,1], header=TRUE)
    if(nrow(tmp_file_all)>0){
      tmp_file_all$sample_id <- file_tbl[i,2]
      final_table_all <- rbind(final_table_all, tmp_file_all)
    }
  }
  if(grepl("unstable",file_tbl[i,1])){
    tmp_file_unstable <- read.table(file_tbl[i,1], header=TRUE)
    if(nrow(tmp_file_unstable)>0){
      tmp_file_unstable$sample_id <- file_tbl[i,2]
      final_table_unstable <- rbind(final_table_unstable, tmp_file_unstable)   
    }
  }
  if(!grepl("unstable",file_tbl[i,1]) & !grepl("all",file_tbl[i,1])){
    print(i)
    tmp_file_summary <- read.table(file_tbl[i,1], header=TRUE, check.names = FALSE)
    if(nrow(tmp_file_summary)>0){
      tmp_file_summary$sample_id <- file_tbl[i,2]
      final_table_summary <- rbind(final_table_summary, tmp_file_summary)
    }
  }
}

#write to outputs
write.table(final_table_all, output_file_all, row.names = FALSE)
write.table(final_table_unstable, output_file_unstable, row.names = FALSE)
write.table(final_table_summary, output_file_summary, row.names = FALSE)
