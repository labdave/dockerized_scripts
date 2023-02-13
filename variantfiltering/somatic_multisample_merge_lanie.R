
#multisample somatic variant calling merger


args = commandArgs(trailingOnly=TRUE)

input.dir<-args[[1]]

merged_file<-args[[3]]

sv.file.list<-list.files(path=input.dir, pattern="*sv_var*", full.names = TRUE)
print(sv.file.list)

mpileup.file.list<-list.files(path=input.dir, pattern="*cleaned_var*", full.names=TRUE)
print(mpileup.file.list)

# remove files if exist
system(paste("rm", args[2]))
system(paste("rm", args[3]))


for(i in 1:length(sv.file.list)) {
  print(i)
  curr.samp<-read.delim(sv.file.list[i], header=TRUE, sep="\t", check.names = FALSE)

  
  ### DEVANG'S CODE ###
  if(i == 1){
    write.table(curr.samp, "temp.txt", sep="\t", row.names=FALSE, quote=FALSE)
  }
  else {
    write.table(curr.samp, "temp.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
  system(paste("cat temp.txt >> ", args[2]))
  ######
} 

all_vars<-NULL
for( i in 1:length(mpileup.file.list)){
  print(i)
  curr.samp<-read.delim(mpileup.file.list[i], header=TRUE, sep="\t", check.names = FALSE)

  
all_vars<-rbind(all_vars, curr.samp)
}

write.table(all_vars, file=merged_file, sep="\t", row.names=FALSE, quote=FALSE)
