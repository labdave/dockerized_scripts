library(vcfR)
library(parallel)
library(stringr)
library(plyr)
library(bedr)

args = commandArgs(trailingOnly=TRUE)

####
samp<-args[2]
path<-args[1]

####


all.data<-read.vcfR(file=path)


all.fix<-data.frame(all.data@fix, stringsAsFactors = FALSE)
all.gt<-data.frame(all.data@gt, stringsAsFactors = FALSE)
all.info<-all.fix$INFO

##only keep "PASS" variants
all.PASS<-all.fix$FILTER == "PASS"
all.info<-all.info[all.PASS]
all.gt<-all.gt[all.PASS,]
all.fix<-all.fix[all.PASS,]



CHROM_POS_REF_ALT<-paste(all.fix$CHROM, all.fix$POS, all.fix$REF, all.fix$ALT, sep="-")




###parse INFO column
tmp<-mclapply(all.info, function(x){
  
  num.values<-str_count(x, ";")
  parsed<-str_split_fixed(x, ";", n=num.values+1)
  
  colnames(parsed)<-unlist(lapply(parsed, function(x){strsplit(x, "=")[[1]][1]}))
  parsed<-gsub(".*=", "", parsed)
  parsed<-data.frame(parsed, stringsAsFactors = FALSE)
  return(parsed)
})

all.info.parsed<-rbind.fill(tmp)


if(!("ANNOVAR_DATE" %in% colnames(all.info.parsed))){
  all.info.hc<-cbind(all.info.parsed, all.fix[,6:7])
  colnames(all.info.hc)<-paste0(colnames(all.info.hc))
  all.fix.hc<-all.fix[,1:5]
}

if("ANNOVAR_DATE" %in% colnames(all.info.parsed)){
  annovar.start<-which(colnames(all.info.parsed)=="ANNOVAR_DATE")
  annovar.end<-which(colnames(all.info.parsed)=="ALLELE_END")
  if(annovar.end == ncol(all.info.parsed)){
    all.info.annovar<-all.info.parsed[,annovar.start:ncol(all.info.parsed)]
    all.info.hc<-cbind(all.info.parsed[,1:(annovar.start-1)], all.fix[,6:7])
  }
  if(annovar.end != ncol(all.info.parsed)){
    all.info.annovar<-all.info.parsed[,annovar.start:annovar.end]
    all.info.hc<-cbind(all.info.parsed[,1:(annovar.start-1)], all.info.parsed[,(annovar.end+1):ncol(all.info.parsed)])
    all.info.hc<-cbind(all.info.hc, all.fix[,6:7])
  }
  
  colnames(all.info.hc)<-paste0(colnames(all.info.hc))
  all.fix.hc<-cbind(all.fix[,1:5], all.info.annovar)
  
  
}


tmp<-mclapply(1:nrow(all.gt), function(x){
  sel<-all.gt[x,]
  num.values<-str_count(sel[2], ":")
  parsed<-as.vector(t(str_split_fixed(sel[1,2:ncol(sel)], ":", n=num.values+1)))
  col.type<-str_split_fixed(sel[1], ":", n=num.values+1)  
  parsed<-as.data.frame(t(parsed), stringsAsFactors = FALSE)
  new.names<-paste("HC", col.type, sep="_")
  colnames(parsed)<-new.names
  
  parsed<-data.frame(parsed, stringsAsFactors = FALSE)
  
  return(parsed)
})

all.gt.parsed.hc<-rbind.fill(tmp)
print(paste("Parsed ", ncol(all.gt.parsed.hc), " GT columns."))

all.fix.merged<-cbind(CHROM_POS_REF_ALT, all.fix.hc)
all.info.merged<-cbind(CHROM_POS_REF_ALT, all.info.hc)
all.gt.merged<-cbind(CHROM_POS_REF_ALT, all.gt.parsed.hc)

print("After initial parse:")
print("fix.filt:")
print(head(all.fix.merged$CHROM_POS_REF_ALT))
print("info.filt:")
print(head(all.info.merged$CHROM_POS_REF_ALT))
print("gt.filt:")
print(head(all.gt.merged$CHROM_POS_REF_ALT))

all.info.merged<-all.info.merged[!grepl("ERCC|GL|KI|KN19|KQ0|KN5|KV8|KV7|KZ5|KZ2|JH", all.info.merged$CHROM_POS_REF_ALT),]
all.gt.merged<-all.gt.merged[!grepl("ERCC|GL|KI|KN19|KQ0|KN5|KV8|KV7|KZ5|KZ2|JH", all.gt.merged$CHROM_POS_REF_ALT),]
all.fix.merged<-all.fix.merged[!grepl("ERCC|GL|KI|KN19|KQ0|KN5|KV8|KV7|KZ5|KZ2|JH", all.fix.merged$CHROM_POS_REF_ALT),]

all.fix.merged$REF<-as.character(all.fix.merged$REF)
all.fix.merged$ALT<-as.character(all.fix.merged$ALT)

####SEPARATE OUT SV CALLS#####
all.gt.sv<-all.gt.merged[!(is.na(all.fix.merged$ID)),]
all.info.sv<-all.info.merged[!(is.na(all.fix.merged$ID)),]
all.fix.sv<-all.fix.merged[!(is.na(all.fix.merged$ID)),]


all.gt.merged<-all.gt.merged[(is.na(all.fix.merged$ID)),]
all.info.merged<-all.info.merged[(is.na(all.fix.merged$ID)),]
all.fix.merged<-all.fix.merged[(is.na(all.fix.merged$ID)),]

print("After removing SVs:")
print("fix.filt:")
print(head(all.fix.merged$CHROM_POS_REF_ALT))
print("info.filt:")
print(head(all.info.merged$CHROM_POS_REF_ALT))
print("gt.filt:")
print(head(all.gt.merged$CHROM_POS_REF_ALT))



pop.freq.data<-suppressWarnings(sapply(all.fix.merged[,grep("gnomad|ExAC|esp6500|X1000g",colnames(all.fix.merged),value=T)], as.numeric))
pop.freq.data[is.na(pop.freq.data)] <- 0
all.fix.merged$pop.freq.max.all <- apply(pop.freq.data,1,max)
rm(pop.freq.data)

all.fix.merged<-all.fix.merged[order(all.fix.merged$CHROM_POS_REF_ALT),]
all.gt.merged<-all.gt.merged[order(all.gt.merged$CHROM_POS_REF_ALT),]
all.info.merged<-all.info.merged[order(all.info.merged$CHROM_POS_REF_ALT),]


print("After ordering before repeatmasker:")
print("fix.filt:")
print(head(all.fix.merged$CHROM_POS_REF_ALT))
print("info.filt:")
print(head(all.info.merged$CHROM_POS_REF_ALT))
print("gt.filt:")
print(head(all.gt.merged$CHROM_POS_REF_ALT))


all.filters<-data.frame(all.fix.merged$CHROM_POS_REF_ALT, stringsAsFactors = FALSE)
all.filters$exonic <- all.fix.merged$Func.refGene=="exonic"
load("filtered_whitelist_09302020.RData")
#load("variantFilter/filtered_whitelist_09302020.RData")
all.filters$whitelist<-all.fix.merged$CHROM_POS_REF_ALT %in% wl$CHROM_POS_REF_ALT
all.filters$rare.variants <-all.fix.merged$pop.freq.max.all<=0.01
all.filters$not.syn <- all.fix.merged[,'ExonicFunc.refGene']!="synonymous_SNV" 
all.filters$not.superdups <- all.fix.merged$genomicSuperDups=="."
all.filters$basic.filters<- (all.filters$whitelist) | (all.filters$rare & all.filters$not.superdups)  
all.filters$pop.freq.max.all<-all.fix.merged$pop.freq.max.all

load("repeat_masker.RData")
options(scipen = 999)

print(head(all.fix.merged$REF))

end.pos<-unlist(lapply(1:nrow(all.fix.merged), function(x){as.numeric(paste(all.fix.merged$POS[x])) + max(nchar(as.character(all.fix.merged$REF[x])), nchar(as.character(all.fix.merged$ALT[x])))}))

var.locs<-c(paste0(all.fix.merged$CHROM, ":", all.fix.merged$POS,"-", end.pos))
#var.locs<-cbind(all.fix.merged[,1:2], all.fix.merged[,2]+1)
table(is.valid.region(var.locs))
var.locs<-bedr.sort.region(var.locs)


not.repeatmasker<-!(in.region(var.locs, bed.file.filt))

names(not.repeatmasker)<-var.locs
not.repeatmasker<-not.repeatmasker[order(names(not.repeatmasker))]
head(all.filters$all.fix.merged.CHROM_POS_REF_ALT)
head(names(not.repeatmasker))
all.filters$not.repeatmasker<-not.repeatmasker

all.info.merged<-cbind(all.info.merged, all.filters[,2:9])


all.fix.filt<-all.fix.merged[ all.info.merged$basic.filters & all.info.merged$not.repeatmasker,]
all.gt.filt<-all.gt.merged[ all.info.merged$basic.filters & all.info.merged$not.repeatmasker,]
all.info.filt<-all.info.merged[ all.info.merged$basic.filters & all.info.merged$not.repeatmasker,]

single.sample.merged<-list(all.fix.merged, all.info.merged, all.gt.merged, all.fix.filt, all.info.filt, all.gt.filt)

if(nrow(single.sample.merged[[4]])==0){
  
  print("No variants pass filtering")
  tmp<-read.csv(text="RNA_depth_total,RNA_depth_alt,RNA_AF,RNA_evidence,Sample_ID")
  all_filt_variants<-cbind(single.sample.merged[[4]], single.sample.merged[[5]], single.sample.merged[[6]], tmp)
}

if(nrow(single.sample.merged[[4]])>0){
  
  all_filt_variants<-cbind(single.sample.merged[[4]], single.sample.merged[[5]], single.sample.merged[[6]])
  
  all_filt_variants <- data.frame(lapply(all_filt_variants, function(x) gsub("\\\\x3d", ":", x)), stringsAsFactors = FALSE)
  all_filt_variants <- data.frame(lapply(all_filt_variants, function(x) gsub("\\\\x3b", "=", x)), stringsAsFactors = FALSE)
  all_filt_variants$Sample_ID<-samp
}  

single.sample.all<-cbind(single.sample.merged[[1]], single.sample.merged[[2]], single.sample.merged[[3]])

single.sample.sv<-cbind(all.fix.sv, all.info.sv, all.gt.sv)

single.sample.all <- data.frame(lapply(single.sample.all, function(x) gsub("\\\\x3d", ":", x)), stringsAsFactors = FALSE)
single.sample.all <- data.frame(lapply(single.sample.all, function(x) gsub("\\\\x3b", "=", x)), stringsAsFactors = FALSE)

single.sample.all$Sample_ID<-samp


clean_filtered<-all_filt_variants
print(paste("clean_filtered size in the beginning:", dim(clean_filtered)))
clean_filtered <- clean_filtered[grepl("exonic", clean_filtered$Func.refGene) | 
                                   grepl("splicing", clean_filtered$Func.refGene) |
                                   grepl("exonic", clean_filtered$Func.ensGene) | 
                                   grepl("splicing", clean_filtered$Func.ensGene), ]
print(paste("clean_filtered size after exonic & splicing:", dim(clean_filtered)))
clean_filtered <- clean_filtered[!(clean_filtered$ExonicFunc.refGene=="synonymous_SNV"),]
print(paste("clean_filtered size after synonymous_SNV:", dim(clean_filtered)))

alt_depth<-as.numeric(unlist(lapply(clean_filtered$HC_AD, function(x){strsplit(x, ",")[[1]][2]})))
print(paste("clean_filtered size before max depth:", dim(clean_filtered)))
clean_filtered<-clean_filtered[alt_depth>=5,]
print(paste("clean_filtered size after max depth:", dim(clean_filtered)))

########put problem columns at end of file
bad_cols<-c("SVTYPE", "MATEID", "SOMATIC", "HC_BaseQRankSumPS", "HC_ClippingRankSumPS", "HC_MQRankSumPS")

single.sample.all<-cbind(single.sample.all[,!(colnames(single.sample.all) %in% bad_cols)], single.sample.all[,bad_cols])
single.sample.sv<-cbind(single.sample.sv[,!(colnames(single.sample.sv) %in% bad_cols)], single.sample.sv[,bad_cols])
clean_filtered<-cbind(clean_filtered[,!(colnames(clean_filtered) %in% bad_cols)], clean_filtered[, bad_cols])


write.table(single.sample.all, file=args[3], row.names  =FALSE, quote=FALSE, sep = "\t")
write.table(single.sample.sv, file=args[4], row.names  =FALSE, quote=FALSE, sep = "\t")
write.table(clean_filtered, file=args[5], row.names  =FALSE, quote=FALSE, sep = "\t")
