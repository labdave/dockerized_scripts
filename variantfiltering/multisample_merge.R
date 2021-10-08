
#multisample variant calling merger


args = commandArgs(trailingOnly=TRUE)

input.dir<-args[[1]]

filt.file.list<-list.files(path=input.dir, pattern="*filt_var*", full.names = TRUE)
print(filt.file.list)

wl.file.list<-list.files(path=input.dir, pattern="*wl_var*", full.names=TRUE)
print(wl.file.list)

cleaned.file.list<-list.files(path=input.dir, pattern="*cleaned_filtered_var*", full.names=TRUE)
print(cleaned.file.list)

mpileup.file.list<-list.files(path=input.dir, pattern="*cleaned_var*", full.names=TRUE)
print(mpileup.file.list)

#define misbehaving columns - columns that only exist for indels --> put at end of dataframe

bad_cols<-c("SNVHPOL", "HC_BaseQRankSum", "HC_ClippingRankSum", "HC_MQRankSum", "HC_ReadPosRankSum", "S2_SNVHPOL", "S2_CIGAR", "S2_RU", "S2_REFREP", "S2_IDREP", "S2_DPI", "S2_AD", "S2_ADF", "S2_ADR", "S2_FT", "S2_PL", "S2_PS")
really_bad_col<-c("DV_VAF.PL")

# remove files if exist
system(paste("rm", args[2]))
system(paste("rm", args[3]))
system(paste("rm", args[4]))
system(paste("rm", args[5]))

for(i in 1:length(filt.file.list)) {
  print(i)
  curr.samp<-read.delim(filt.file.list[i], header=TRUE, sep="\t", check.names = FALSE)
  curr.samp<-curr.samp[,!(colnames(curr.samp) %in% really_bad_col)]
 
  curr.samp.good<-curr.samp[,!(colnames(curr.samp) %in% bad_cols)]
  curr.samp.bad<-curr.samp[,colnames(curr.samp) %in% bad_cols]
  
  if("SNVHPOL" %in% colnames(curr.samp.bad)){
    #colnames(curr.samp.bad)[colnames(curr.samp.bad)=="SNVHPOL"]<-"S2_SNVHPOL"
    
    correct_order<-cbind(curr.samp.bad[,2:5], curr.samp.bad[,1], curr.samp.bad[,11:15], curr.samp.bad[,10], curr.samp.bad[,16])
    curr.samp.bad<-correct_order
    colnames(curr.samp.bad)[5]<-"S2_SNVHPOL"
    colnames(curr.samp.bad)[11]<-"S2_DPI"
    colnames(curr.samp.bad)[12]<-"S2_PS"
  }
  
  curr.samp<-cbind(curr.samp.good, curr.samp.bad)

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


for( i in 1:length(wl.file.list)){
  print(i)
  curr.samp<-read.delim(wl.file.list[i], header=TRUE, sep="\t", check.names = FALSE)
  curr.samp<-curr.samp[,!(colnames(curr.samp) %in% really_bad_col)]
  
  curr.samp.good<-curr.samp[,!(colnames(curr.samp) %in% bad_cols)]
  curr.samp.bad<-curr.samp[,colnames(curr.samp) %in% bad_cols]
  
  if("SNVHPOL" %in% colnames(curr.samp.bad)){
    #colnames(curr.samp.bad)[colnames(curr.samp.bad)=="SNVHPOL"]<-"S2_SNVHPOL"
    
    correct_order<-cbind(curr.samp.bad[,2:5], curr.samp.bad[,1], curr.samp.bad[,11:15], curr.samp.bad[,10], curr.samp.bad[,16])
    curr.samp.bad<-correct_order
    colnames(curr.samp.bad)[5]<-"S2_SNVHPOL"
    colnames(curr.samp.bad)[11]<-"S2_DPI"
    colnames(curr.samp.bad)[12]<-"S2_PS"
  }
  
  curr.samp<-cbind(curr.samp.good, curr.samp.bad)
  
  ### DEVANG'S CODE ###
  if(i == 1){
    write.table(curr.samp, "temp.txt", sep="\t", row.names=FALSE, quote=FALSE)
  }
  else {
    write.table(curr.samp, "temp.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
  system(paste("cat temp.txt >> ", args[3]))
  ######
}


for( i in 1:length(cleaned.file.list)){
  print(i)
  curr.samp<-read.delim(cleaned.file.list[i], header=TRUE, sep="\t", check.names = FALSE)
  curr.samp<-curr.samp[,!(colnames(curr.samp) %in% really_bad_col)]
 
  curr.samp.good<-curr.samp[,!(colnames(curr.samp) %in% bad_cols)]
  curr.samp.bad<-curr.samp[,colnames(curr.samp) %in% bad_cols]
  
  if("SNVHPOL" %in% colnames(curr.samp.bad)){
    #colnames(curr.samp.bad)[colnames(curr.samp.bad)=="SNVHPOL"]<-"S2_SNVHPOL"
    
    correct_order<-cbind(curr.samp.bad[,2:5], curr.samp.bad[,1], curr.samp.bad[,11:15], curr.samp.bad[,10], curr.samp.bad[,16])
    curr.samp.bad<-correct_order
    colnames(curr.samp.bad)[5]<-"S2_SNVHPOL"
    colnames(curr.samp.bad)[11]<-"S2_DPI"
    colnames(curr.samp.bad)[12]<-"S2_PS"
  }
  
  curr.samp<-cbind(curr.samp.good, curr.samp.bad)
  
  ### DEVANG'S CODE ###
  if(i == 1){
    write.table(curr.samp, "temp.txt", sep="\t", row.names=FALSE, quote=FALSE)
  }
  else {
    write.table(curr.samp, "temp.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
  system(paste("cat temp.txt >> ", args[4]))
  ######
}


for( i in 1:length(mpileup.file.list)){
  print(i)
  curr.samp<-read.delim(mpileup.file.list[i], header=TRUE, sep="\t", check.names = FALSE)
  curr.samp<-curr.samp[,!(colnames(curr.samp) %in% really_bad_col)]
 
  curr.samp.good<-curr.samp[,!(colnames(curr.samp) %in% bad_cols)]
  curr.samp.bad<-curr.samp[,colnames(curr.samp) %in% bad_cols]
  
  if("SNVHPOL" %in% colnames(curr.samp.bad)){
    #colnames(curr.samp.bad)[colnames(curr.samp.bad)=="SNVHPOL"]<-"S2_SNVHPOL"
    
    correct_order<-cbind(curr.samp.bad[,2:5], curr.samp.bad[,1], curr.samp.bad[,11:15], curr.samp.bad[,10], curr.samp.bad[,16])
    curr.samp.bad<-correct_order
    colnames(curr.samp.bad)[5]<-"S2_SNVHPOL"
    colnames(curr.samp.bad)[11]<-"S2_DPI"
    colnames(curr.samp.bad)[12]<-"S2_PS"
  }
  
  curr.samp<-cbind(curr.samp.good, curr.samp.bad)
  
  ### DEVANG'S CODE ###
  if(i == 1){
    write.table(curr.samp, "temp.txt", sep="\t", row.names=FALSE, quote=FALSE)
  }
  else {
    write.table(curr.samp, "temp.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
  system(paste("cat temp.txt >> ", args[5]))
  ######    
}

