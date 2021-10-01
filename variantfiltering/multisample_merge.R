
#multisample variant calling merger


args = commandArgs(trailingOnly=TRUE)

input.dir<-args[[1]]

filt.file.list<-list.files(path=input.dir, pattern="*filt_var*", full.names = TRUE)
print(filt.file.list)

wl.file.list<-list.files(path=input.dir, pattern="*wl_var*", full.names=TRUE)
print(wl.file.list)

cleaned.file.list<-list.files(path=input.dir, pattern="*cleaned_var*", full.names=TRUE)
print(cleaned.file.list)

#define misbehaving columns - columns that only exist for indels --> put at end of dataframe

bad_cols<-c("SNVHPOL", "HC_BaseQRankSum", "HC_ClippingRankSum", "HC_MQRankSum", "HC_ReadPosRankSum", "S2_SNVHPOL", "S2_CIGAR", "S2_RU", "S2_REFREP", "S2_IDREP", "S2_DPI", "S2_AD", "S2_ADF", "S2_ADR", "S2_FT", "S2_PL", "S2_PS")
really_bad_col<-c("DV_VAF.PL")

for( i in 1:length(filt.file.list)){
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
  
    
  #long form - just concatenate all rows, each row represents one sample-variant pair
  if(exists("long.form")){
    
    ##check to make sure correct number of columns, if not, report but continue merging other files 
    if(ncol(curr.samp)==ncol(long.form)){
      long.form<-rbind(long.form, curr.samp)
    }
    if(ncol(curr.samp)!=ncol(long.form)){
      print(paste0("Wrong number of columns in output file: ", filt.file.list[i]))
    }
  }
  if(!(exists("long.form"))){
    long.form<-curr.samp
  }  

  
  #wide form - merge variants - each row represents one variants and each sample has their own columns (nCallers, AF, depth)
  if(exists("wide.form.fixed")){
    
    fixed<-curr.samp[,1:139]
    samp<-cbind(curr.samp$CHROM_POS_REF_ALT, curr.samp[,grepl("nCallers|afMax|dpMax", colnames(curr.samp))])
    id<-curr.samp$Sample_ID[1]
    
    colnames(samp)[1]<-"CHROM_POS_REF_ALT"
    colnames(samp)[2:4]<-paste(id, colnames(samp[,2:4]))
    
    wide.form.fixed<-unique(rbind(wide.form.fixed, fixed))
    
    wide.form.samp<-merge(x=wide.form.samp, y=samp, by="CHROM_POS_REF_ALT", all.x=TRUE, all.y=TRUE)
    
  }
   if(!(exists("wide.form.fixed"))){
    
    wide.form.fixed<-curr.samp[,1:139]
    samp<-cbind(curr.samp$CHROM_POS_REF_ALT, curr.samp[,grepl("nCallers|afMax|dpMax", colnames(curr.samp))])
    id<-curr.samp$Sample_ID[1]
    
    colnames(samp)[1]<-"CHROM_POS_REF_ALT"
    colnames(samp)[2:4]<-paste(id, colnames(samp[,2:4]))
    
    wide.form.samp<-samp
  }  
  
  
}

wide.form.fixed<-wide.form.fixed[order(wide.form.fixed$CHROM_POS_REF_ALT),]
wide.form.samp<-wide.form.samp[order(wide.form.samp$CHROM_POS_REF_ALT),]

wide.form<-merge(x=wide.form.fixed, y=wide.form.samp, by="CHROM_POS_REF_ALT")

write.table(long.form, file=args[2], row.names = FALSE, quote = FALSE, sep="\t")
write.table(wide.form, file=args[3], row.names = FALSE, quote = FALSE, sep="\t")



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
  
  if(nrow(curr.samp)>0){
    #long form - just concatenate all rows, each row represents one sample-variant pair
    if(exists("wl.long.form")){
      wl.long.form<-rbind(wl.long.form, curr.samp)
    }
    if(!(exists("wl.long.form"))){
      wl.long.form<-curr.samp
    }  
    
    
    #wide form - merge variants - each row represents one variants and each sample has their own columns (nCallers, AF, depth)
    if(exists("wl.wide.form.fixed")){
      
      fixed<-curr.samp[,1:139]
      samp<-cbind(curr.samp$CHROM_POS_REF_ALT, curr.samp[,grepl("nCallers|afMax|dpMax|DNA_", colnames(curr.samp))])
      id<-curr.samp$Sample_ID[1]
      
      colnames(samp)[1]<-"CHROM_POS_REF_ALT"
      colnames(samp)[2:7]<-paste(id, colnames(samp[,2:7]))
      
      wl.wide.form.fixed<-unique(rbind(wl.wide.form.fixed, fixed))
      
      wl.wide.form.samp<-merge(x=wl.wide.form.samp, y=samp, by="CHROM_POS_REF_ALT", all.x=TRUE, all.y=TRUE)
      
    }
    if(!(exists("wl.wide.form.fixed"))){
      
      wl.wide.form.fixed<-curr.samp[,1:139]
      samp<-cbind(curr.samp$CHROM_POS_REF_ALT, curr.samp[,grepl("nCallers|afMax|dpMax|DNA_", colnames(curr.samp))])
      id<-curr.samp$Sample_ID[1]
      
      colnames(samp)[1]<-"CHROM_POS_REF_ALT"
      colnames(samp)[2:7]<-paste(id, colnames(samp[,2:7]))
      
      wl.wide.form.samp<-samp
    }  
  }
  
  
  
}

wl.wide.form.fixed<-wl.wide.form.fixed[order(wl.wide.form.fixed$CHROM_POS_REF_ALT),]
wl.wide.form.samp<-wl.wide.form.samp[order(wl.wide.form.samp$CHROM_POS_REF_ALT),]

wl.wide.form<-merge(x=wl.wide.form.fixed, y=wl.wide.form.samp, by="CHROM_POS_REF_ALT")

write.table(wl.long.form, file=args[4], row.names = FALSE, quote = FALSE, sep="\t")
write.table(wl.wide.form, file=args[5], row.names = FALSE, quote = FALSE, sep="\t")



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
  
    
  #long form - just concatenate all rows, each row represents one sample-variant pair
  if(exists("long.form")){
    
    ##check to make sure correct number of columns, if not, report but continue merging other files 
    if(ncol(curr.samp)==ncol(long.form)){
      long.form<-rbind(long.form, curr.samp)
    }
    if(ncol(curr.samp)!=ncol(long.form)){
      print(paste0("Wrong number of columns in output file: ", cleaned.file.list[i]))
    }
  }
  if(!(exists("long.form"))){
    long.form<-curr.samp
  }  

  
  #wide form - merge variants - each row represents one variants and each sample has their own columns (nCallers, AF, depth)
  if(exists("wide.form.fixed")){
    
    fixed<-curr.samp[,1:144]
    samp<-cbind(curr.samp$CHROM_POS_REF_ALT, curr.samp[,grepl("nCallers|afMax|dpMax", colnames(curr.samp))])
    id<-curr.samp$Sample_ID[1]
    
    colnames(samp)[1]<-"CHROM_POS_REF_ALT"
    colnames(samp)[2:4]<-paste(id, colnames(samp[,2:4]))
    
    wide.form.fixed<-unique(rbind(wide.form.fixed, fixed))
    
    wide.form.samp<-merge(x=wide.form.samp, y=samp, by="CHROM_POS_REF_ALT", all.x=TRUE, all.y=TRUE)
    
  }
   if(!(exists("wide.form.fixed"))){
    
    wide.form.fixed<-curr.samp[,1:144]
    samp<-cbind(curr.samp$CHROM_POS_REF_ALT, curr.samp[,grepl("nCallers|afMax|dpMax", colnames(curr.samp))])
    id<-curr.samp$Sample_ID[1]
    
    colnames(samp)[1]<-"CHROM_POS_REF_ALT"
    colnames(samp)[2:4]<-paste(id, colnames(samp[,2:4]))
    
    wide.form.samp<-samp
  } 
}

wide.form.fixed<-wide.form.fixed[order(wide.form.fixed$CHROM_POS_REF_ALT),]
wide.form.samp<-wide.form.samp[order(wide.form.samp$CHROM_POS_REF_ALT),]

wide.form<-merge(x=wide.form.fixed, y=wide.form.samp, by="CHROM_POS_REF_ALT")

write.table(long.form, file=args[6], row.names = FALSE, quote = FALSE, sep="\t")
write.table(wide.form, file=args[7], row.names = FALSE, quote = FALSE, sep="\t")
