

geno_file<-args[[1]]
filt_snv_file<-args[[2]]
filt_indel_file<-args[[3]]
var_coverage_out<-[[4]]

all_genotype<-read.delim(geno_file)


all_genotype$CHROM_POS_REF_ALT<-paste0(all_genotype$Chrom, "-", all_genotype$Start, "-", all_genotype$Ref, "-", all_genotype$Alt)

myeloid_genes<-unlist(read.csv("myeloid_gene_panel_v2.1_121622.txt", header=F, sep="\t"))
lymphoid_genes<-unlist(read.csv("lymphoid_gene_panel_v2.1_121622.txt", header=F, sep="\t"))
new_myeloid<-unlist(read.csv("myeloid_gene_panel_v3.1_020723.txt", header=F, sep="\t"))

all_genotype$myeloid<-NA
all_genotype$myeloid[all_genotype$Gene %in% myeloid_genes]<-"yes"
all_genotype$lymphoid<-NA
all_genotype$lymphoid[all_genotype$Gene %in% lymphoid_genes]<-"yes"
all_genotype$new_myeloid<-NA
all_genotype$new_myeloid[all_genotype$Gene %in% new_myeloid_genes]<-"yes"




###reprocess indel calcs

all_genotype$ins<-0
all_genotype$ins_af<-0
all_genotype$del<-0
all_genotype$del_af<-0
snv_sel<-all_genotype[!(grepl("indel|fullgene", all_genotype$AA_change)),]
indel_sel<-all_genotype[grepl("indel|fullgene", all_genotype$AA_change),]

all_chars<-unique(c(unlist(strsplit(indel_sel$dna_full, "")), "-", "+"))
char_counts<-matrix(nrow=nrow(indel_sel), ncol=length(all_chars))
rownames(char_counts)<-paste0(indel_sel$Chrom, "_", indel_sel$Start, "_", indel_sel$Ref, "_", indel_sel$Alt)
colnames(char_counts)<-all_chars


for (i in 1:nrow(indel_sel)){
  sel<-indel_sel[i,]
  
  counts<-table(unlist(strsplit(sel$dna_full, "")))
  for(j in 1:length(counts)){
    char_counts[i,names(counts[j])]<-counts[j]
  }
  
  
}

char_counts[is.na(char_counts)]<-0




for (i in 1:nrow(indel_sel)){
  indel_sel$ins[i]<-char_counts[i, "+"]
  indel_sel$ins_af[i]<-indel_sel$ins[i]/indel_sel$dna_cov[i]
  indel_sel$del[i]<-char_counts[i, "-"]
  indel_sel$del_af[i]<-indel_sel$del[i]/indel_sel$dna_cov[i]
  
}



#write.table(snv_sel, file="genotype_merged_snv.txt",row.names = F, quote = F, sep="\t")
#write.table(indel_sel, file="genotype_merged_indel.txt", row.names = F, quote = F, sep="\t")

snv_filt<-snv_sel
snv_filt<-snv_filt[snv_filt$dna_alt_count>=10 | (snv_filt$dna_alt_count>3 & snv_filt$dna_AF>0.02),]
write.table(snv_filt, file=filt_snv_file, row.names = F, quote = F, sep="\t")

indel_filt<-indel_sel
ins_filt<-indel_filt[(indel_filt$ins>3&indel_filt$ins_af>0.1) | indel_filt$ins>=20,]
ins_filt$type<-"insertion"
del_filt<-indel_filt[(indel_filt$del>3&indel_filt$del_af>0.1) | indel_filt$del>=20,]
del_filt$type<-"deletion"

indel_filt<-rbind(ins_filt, del_filt)

write.table(indel_filt, file=filt_indel_file, row.names = F, quote = F, sep="\t")


moldx_variants<-read.table("moldx_variant_mapping.txt", header=T, sep="\t")

all_genotype$Chrom_Pos<-paste0(all_genotype$Chrom, " ", all_genotype$Start)
moldx_variants$Chrom_Pos<-paste0(moldx_variants$chr, " ", moldx_variants$pos)

samps<-unique(all_genotype$sample)
unique_moldx<-unique(moldx_variants$moldx)

moldx_per_sample<-data.frame(matrix(ncol = length(samps)+2, nrow=length(unique_moldx)))

moldx_per_sample[,1]<-unlist(strsplit(unique_moldx, "-"))[c(T,F)]
moldx_per_sample[,2]<-unlist(strsplit(unique_moldx, "-"))[c(F,T)]

for(i in 1:length(samps)){
  print(i)
  sel_g<-all_genotype[all_genotype$sample==samps[i] & all_genotype$Chrom_Pos %in% moldx_variants$Chrom_Pos,]
  for(j in 1:length(unique_moldx)){
    print(j)
    sel<-moldx_variants$Chrom_Pos[moldx_variants$moldx==unique_moldx[j]]
    sel_val<-unique(sel_g[sel_g$Chrom_Pos %in% sel,c("Chrom_Pos", "dna_cov")])
    m<-mean(sel_val$dna_cov)
    moldx_per_sample[j,i+2]<-m
    
    
  }
}

colnames(moldx_per_sample)<-c("Gene", "Variant",samps)
#rownames(moldx_per_sample)<-unique_moldx

write.table(moldx_per_sample, file=var_coverage_out, sep="\t", quote = F, row.names = F)

