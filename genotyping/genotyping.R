#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

hotspots_file<-args[1]
dna_mpileup_file<-args[2]
rna_mpileup_file<-args[3]
samp<-args[4]


hotspots<-read.table(hotspots_file, sep="\t")
colnames(hotspots)<-c("Chrom", "Start", "End", "Ref", "Alt", "Gene", "AA_change", "OncoKB_annot")
dna_mpileup<-read.csv(dna_mpileup_file, header=F, sep="\t")
colnames(dna_mpileup)<-c("Chrom", "Pos", "Ref", "Depth", "mpileup", "Qual")
rna_mpileup<-read.csv(rna_mpileup_file, header=F, sep="\t")
colnames(rna_mpileup)<-c("Chrom", "Pos", "Ref", "Depth", "mpileup", "Qual")



all_chars<-unique(c(unlist(strsplit(dna_mpileup$mpileup, "")), unlist(strsplit(rna_mpileup$mpileup, "")), "-", "+"))
dna_char_counts<-matrix(nrow=nrow(hotspots), ncol=length(all_chars))
rownames(dna_char_counts)<-paste0(hotspots$Chrom, "_", hotspots$Start, "_", hotspots$Ref, "_", hotspots$Alt)
colnames(dna_char_counts)<-all_chars

rna_char_counts<-matrix(nrow=nrow(hotspots), ncol=length(all_chars))
rownames(rna_char_counts)<-paste0(hotspots$Chrom, "_", hotspots$Start, "_", hotspots$Ref, "_", hotspots$Alt)
colnames(rna_char_counts)<-all_chars

dna_cov<-NULL
rna_cov<-NULL

for (i in 1:nrow(hotspots)){
  sel<-hotspots[i,]
  
  dna_cov<-c(dna_cov, dna_mpileup$Depth[dna_mpileup$Chrom==sel$Chrom & dna_mpileup$Pos==sel$Start])
  dna<-table(unlist(strsplit(dna_mpileup$mpileup[dna_mpileup$Chrom==sel$Chrom & dna_mpileup$Pos==sel$Start], "")))
  for(j in 1:length(dna)){
    dna_char_counts[i,names(dna[j])]<-dna[j]
  }
  
  rna_cov<-c(rna_cov, rna_mpileup$Depth[rna_mpileup$Chrom==sel$Chrom & rna_mpileup$Pos==sel$Start])
  rna<-table(unlist(strsplit(rna_mpileup$mpileup[rna_mpileup$Chrom==sel$Chrom & rna_mpileup$Pos==sel$Start], "")))
  for(j in 1:length(rna)){
    rna_char_counts[i,names(rna[j])]<-rna[j]
  }
  
}

dna_char_counts[is.na(dna_char_counts)]<-0
rna_char_counts[is.na(rna_char_counts)]<-0


hotspots$sample<-rep(samp, nrow(hotspots))
hotspots$dna_cov<-dna_cov
hotspots$dna_ref_count<-0
hotspots$dna_alt_count<-0
hotspots$dna_AF<-0
hotspots$dna_warning<-"."


hotspots$rna_cov<-rna_cov
hotspots$rna_ref_count<-0
hotspots$rna_alt_count<-0
hotspots$rna_AF<-0
hotspots$rna_warning<-"."

hotspots$dna_full<-0
hotspots$rna_full<-0


for (i in 1:nrow(hotspots)){
  ref<-hotspots$Ref[i]
  alt<-hotspots$Alt[i]
  
  hotspots$dna_ref_count[i]<-dna_char_counts[i,ref]+dna_char_counts[i,tolower(ref)]
  hotspots$rna_ref_count[i]<-rna_char_counts[i,ref]+rna_char_counts[i,tolower(ref)]
  
  hotspots$dna_alt_count[i]<-dna_char_counts[i,alt]+dna_char_counts[i,tolower(alt)]
  hotspots$rna_alt_count[i]<-rna_char_counts[i,alt]+rna_char_counts[i,tolower(alt)]
  
  hotspots$dna_AF[i]<-hotspots$dna_alt_count[i]/hotspots$dna_cov[i]
  hotspots$rna_AF[i]<-hotspots$rna_alt_count[i]/hotspots$rna_cov[i]
  
  if(dna_char_counts[i, "-"]>0 | dna_char_counts[i,"+"]>0){
    hotspots$dna_warning[i]<-"Indel possible"
  }
  if(rna_char_counts[i, "-"]>0 | rna_char_counts[i,"+"]>0){
    hotspots$rna_warning[i]<-"Indel possible"
  }
  
  sel_dna<-dna_mpileup[dna_mpileup$Chrom==hotspots$Chrom[i] & dna_mpileup$Pos==hotspots$Start[i],]
  hotspots$dna_full[i]<-sel_dna$mpileup
  sel_rna<-rna_mpileup[rna_mpileup$Chrom==hotspots$Chrom[i] & rna_mpileup$Pos==hotspots$Start[i],]
  hotspots$rna_full[i]<-sel_rna$mpileup
  
}

hotspots[hotspots=="NaN"]<-0

