
#multisample genotyping merger


args = commandArgs(trailingOnly=TRUE)

input.dir<-args[[1]]

merged_file<-args[[2]]

geno.file.list<-list.files(path=input.dir, pattern="*geno*", full.names=TRUE)
print(geno.file.list)

# remove files if exist
system(paste("rm", args[2]))

all_vars<-NULL
for( i in 1:length(geno.file.list)){
  print(i)
  curr.samp<-read.delim(geno.file.list[i], header=TRUE, sep="\t", check.names = FALSE)
  
  
  all_vars<-rbind(all_vars, curr.samp)
}

write.table(all_vars, file=merged_file, sep="\t", row.names=FALSE, quote=FALSE)
