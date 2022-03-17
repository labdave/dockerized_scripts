library("FARDEEP")
library("ggplot2")
library("RColorBrewer")
library("colorspace")
library("ComplexHeatmap")
library("circlize")
library("tidyverse")

args = commandArgs(trailingOnly=TRUE)

######cleanup######
print("read file")
cts <- read.table(args[1], sep="\t", header=1, check.names=FALSE)
sample.names <- colnames(cts)[3:ncol(cts)]
cts <- cts[rowSums(cts[,sample.names]) > 0,]
cts[,3:ncol(cts)] <- cts[,3:ncol(cts)] + 1
cts <- mutate_if(cts, is.numeric, log2)

dup_genes <- names(table(cts$gene_symbol)[table(cts$gene_symbol)>1])
dup_cts <- cts %>% filter(gene_symbol %in% dup_genes)

# Keep row with maximum expression?
dup_cts <- dup_cts[order(rowSums(dup_cts[,sample.names]), decreasing = TRUE),] %>% 
  distinct(gene_symbol, .keep_all = TRUE)

cts <- cts %>% filter(!(gene_symbol %in% dup_genes))
cts <- rbind(cts, dup_cts)

rownames(cts) <- cts$gene_symbol
cts <- cts[,sample.names]

colnames(cts)<-unlist(strsplit(colnames(cts), "_S"))[c(TRUE,FALSE)]

diagnoses <- read.table("who_diagnoses.csv", sep=",", header=1, check.names=FALSE)


### hierarchical cell type classification
print("read markers")
signatures <- list()

who<-read.table("tables/dummy.tsv", sep="\t", header = T, row.names = 1, check.names = F)
signatures[["who"]]<-who

level1<-read.table("tables/level1.tsv", sep="\t", header = T, row.names = 1, check.names = F)
signatures[["level1"]]<-level1

level2_b<-read.table("tables/level2_b.tsv", sep="\t", header = T, row.names = 1, check.names = F)
signatures[["level2_b"]]<-level2_b

level2_mye<-read.table("tables/level2_mye.tsv", sep="\t", header = T, row.names = 1, check.names = F)
signatures[["level2_mye"]]<-level2_mye

level2_t<-read.table("tables/level2_t.tsv", sep="\t", header = T, row.names = 1, check.names = F)
signatures[["level2_t"]]<-level2_t

level3_b_immature<-read.table("tables/level3_b_immature.tsv", sep="\t", header = T, row.names = 1, check.names = F)
signatures[["level3_b_immature"]]<-level3_b_immature

level3_b_mature<-read.table("tables/level3_b_mature.tsv", sep="\t", header = T, row.names = 1, check.names = F)
signatures[["level3_b_mature"]]<-level3_b_mature

level3_mye_other<-read.table("tables/level3_mye_other.tsv", sep="\t", header = T, row.names = 1, check.names = F)
signatures[["level3_mye_other"]]<-level3_mye_other

level3_t_cd4<-read.table("tables/level3_t_cd4.tsv", sep="\t", header = T, row.names = 1, check.names = F)
signatures[["level3_t_cd4"]]<-level3_t_cd4
level3_t_other<-read.table("tables/level3_t_other.tsv", sep="\t", header = T, row.names = 1, check.names = F)
signatures[["level3_t_other"]]<-level3_t_other

level4_b_mature_germinal<-read.table("tables/level4_b_mature_germinal.tsv", sep="\t", header = T, row.names = 1, check.names = F)
signatures[["level4_b_mature_germinal"]]<-level4_b_mature_germinal

# run FARDEEP
print("run FARDEEP")
sh_results <- lapply(signatures, function(sig) {fardeep(sig, cts, permn = 0)})

print("postprocess")
# create a dataframe of zeros
sh_results$who$relative.beta[is.na(sh_results$who$relative.beta)] <- 0
sh_results$who$relative.beta[sh_results$who$relative.beta > -1] <- 0

# set all NAs to zeros
sh_results$level1$relative.beta[is.na(sh_results$level1$relative.beta)] <- 0
sh_results$level2_b$relative.beta[is.na(sh_results$level2_b$relative.beta)] <- 0
sh_results$level2_t$relative.beta[is.na(sh_results$level2_t$relative.beta)] <- 0
sh_results$level2_mye$relative.beta[is.na(sh_results$level2_mye$relative.beta)] <- 0
sh_results$level3_b_immature$relative.beta[is.na(sh_results$level3_b_immature$relative.beta)] <- 0
sh_results$level3_b_mature$relative.beta[is.na(sh_results$level3_b_mature$relative.beta)] <- 0
sh_results$level3_t_nk$relative.beta[is.na(sh_results$level3_t_nk$relative.beta)] <- 0
sh_results$level3_t_cd4$relative.beta[is.na(sh_results$level3_t_cd4$relative.beta)] <- 0
sh_results$level3_t_other$relative.beta[is.na(sh_results$level3_t_other$relative.beta)] <- 0
sh_results$level3_mye_other$relative.beta[is.na(sh_results$level3_mye_other$relative.beta)] <- 0
sh_results$level4_b_mature_germinal$relative.beta[is.na(sh_results$level4_b_mature_germinal$relative.beta)] <- 0
sh_results$level4_mye_other_mast$relative.beta[is.na(sh_results$level4_mye_other_mast$relative.beta)] <- 0
sh_results$level4_t_cd4_memory$relative.beta[is.na(sh_results$level4_t_cd4_memory$relative.beta)] <- 0

# multiply hierarchical values and set in full dataframe
for (row in 1:nrow(sh_results$who$relative.beta)) {
  
  # b cell subgroups
  sh_results$who$relative.beta[row, "B.immature..pro-b"]                = sh_results$level1$relative.beta[row, "B"]*sh_results$level2_b$relative.beta[row, "B.immature"]*sh_results$level3_b_immature$relative.beta[row, "B.immature..pro-b"]
  sh_results$who$relative.beta[row, "B.immature..pre-b"]                = sh_results$level1$relative.beta[row, "B"]*sh_results$level2_b$relative.beta[row, "B.immature"]*sh_results$level3_b_immature$relative.beta[row, "B.immature..pre-b"]
  sh_results$who$relative.beta[row, "B.mature..naive"]                  = sh_results$level1$relative.beta[row, "B"]*sh_results$level2_b$relative.beta[row, "B.mature"]*sh_results$level3_b_mature$relative.beta[row, "B.mature..naive"]
  sh_results$who$relative.beta[row, "B.mature..germinal...dz"]          = sh_results$level1$relative.beta[row, "B"]*sh_results$level2_b$relative.beta[row, "B.mature"]*sh_results$level3_b_mature$relative.beta[row, "B.mature..germinal"]*sh_results$level4_b_mature_germinal$relative.beta[row, "B.mature..germinal...dz"]
  sh_results$who$relative.beta[row, "B.mature..germinal...lz"]          = sh_results$level1$relative.beta[row, "B"]*sh_results$level2_b$relative.beta[row, "B.mature"]*sh_results$level3_b_mature$relative.beta[row, "B.mature..germinal"]*sh_results$level4_b_mature_germinal$relative.beta[row, "B.mature..germinal...lz"]
  sh_results$who$relative.beta[row, "B.mature..memory"]                 = sh_results$level1$relative.beta[row, "B"]*sh_results$level2_b$relative.beta[row, "B.mature"]*sh_results$level3_b_mature$relative.beta[row, "B.mature..memory"]
  sh_results$who$relative.beta[row, "B.mature..plasma"]                 = sh_results$level1$relative.beta[row, "B"]*sh_results$level2_b$relative.beta[row, "B.mature"]*sh_results$level3_b_mature$relative.beta[row, "B.mature..plasma"]
  
  # t cell subgroups
  sh_results$who$relative.beta[row, "T.CD8"]                            = sh_results$level1$relative.beta[row, "T"]*sh_results$level2_t$relative.beta[row, "T.CD8"]
  sh_results$who$relative.beta[row, "T.follicular_helper"]              = sh_results$level1$relative.beta[row, "T"]*sh_results$level2_t$relative.beta[row, "T.follicular_helper"]
  sh_results$who$relative.beta[row, "T.CD4..naive"]                     = sh_results$level1$relative.beta[row, "T"]*sh_results$level2_t$relative.beta[row, "T.CD4"]*sh_results$level3_t_cd4$relative.beta[row, "T.CD4..naive"]
  sh_results$who$relative.beta[row, "T.CD4..memory"]                    = sh_results$level1$relative.beta[row, "T"]*sh_results$level2_t$relative.beta[row, "T.CD4"]*sh_results$level3_t_cd4$relative.beta[row, "T.CD4..memory"]
  sh_results$who$relative.beta[row, "T.NK"]                             = sh_results$level1$relative.beta[row, "T"]*sh_results$level2_t$relative.beta[row, "T.NK"]
  sh_results$who$relative.beta[row, "T.other..regulatory"]              = sh_results$level1$relative.beta[row, "T"]*sh_results$level2_t$relative.beta[row, "T.other"]*sh_results$level3_t_other$relative.beta[row, "T.other..regulatory"]
  sh_results$who$relative.beta[row, "T.other..gamma_delta"]             = sh_results$level1$relative.beta[row, "T"]*sh_results$level2_t$relative.beta[row, "T.other"]*sh_results$level3_t_other$relative.beta[row, "T.other..gamma_delta"]
  sh_results$who$relative.beta[row, "T.plasmacytoid"]                   = sh_results$level1$relative.beta[row, "T"]*sh_results$level2_t$relative.beta[row, "T.plasmacytoid"]
  
  # myeloid cell subgroups
  sh_results$who$relative.beta[row, "Myeloid.langerhans"]               = sh_results$level1$relative.beta[row, "Myeloid"]*sh_results$level2_mye$relative.beta[row, "Myeloid.langerhans"]
  sh_results$who$relative.beta[row, "Myeloid.dendritic"]                = sh_results$level1$relative.beta[row, "Myeloid"]*sh_results$level2_mye$relative.beta[row, "Myeloid.dendritic"]
  sh_results$who$relative.beta[row, "Myeloid.other..eosinophils"]       = sh_results$level1$relative.beta[row, "Myeloid"]*sh_results$level2_mye$relative.beta[row, "Myeloid.other"]*sh_results$level3_mye_other$relative.beta[row, "Myeloid.other..eosinophils"]
  sh_results$who$relative.beta[row, "Myeloid.other..neutrophils"]       = sh_results$level1$relative.beta[row, "Myeloid"]*sh_results$level2_mye$relative.beta[row, "Myeloid.other"]*sh_results$level3_mye_other$relative.beta[row, "Myeloid.other..neutrophils"]
  sh_results$who$relative.beta[row, "Myeloid.other..monocytes"]         = sh_results$level1$relative.beta[row, "Myeloid"]*sh_results$level2_mye$relative.beta[row, "Myeloid.other"]*sh_results$level3_mye_other$relative.beta[row, "Myeloid.other..monocytes"]
  sh_results$who$relative.beta[row, "Myeloid.other..mast"]              = sh_results$level1$relative.beta[row, "Myeloid"]*sh_results$level2_mye$relative.beta[row, "Myeloid.other"]*sh_results$level3_mye_other$relative.beta[row, "Myeloid.other..mast"]

  # scale up in case of NA induced zeros
  if(sum(sh_results$who$relative.beta[row,]) < 0.99 & sum(sh_results$who$relative.beta[row,]) > 0.01) {
    print(row)
    print(sh_results$who$relative.beta[row,])
    scale = 1/sum(sh_results$who$relative.beta[row,])
    for(i in colnames(sh_results$who$relative.beta)){
      sh_results$who$relative.beta[row, i] = sh_results$who$relative.beta[row, i]*scale
    }
  }
}

# add diagnoses
new_df <- data.frame(sh_results$who$relative.beta, check.names=FALSE)
new_df["Diagnosis"] <- ""
new_df <- new_df[,c(ncol(new_df),1:(ncol(new_df)-1))]

rownames(new_df) <- sub("^X", "", rownames(new_df))
rownames(diagnoses) <- diagnoses$Patient

for(i in 1:nrow(new_df)){
  new_df[i, "Diagnosis"] = data.frame(diagnoses)[rownames(new_df[i,]), "Diagnosis"]
}

write.table(t(new_df), file=args[2], sep="\t", quote = F, row.names = T, col.names=NA)

heatmap_df = t(sh_results$who$relative.beta)
mappings = read.table("name_mapping.txt", sep="\t", header=1, check.names=FALSE)
rownames(mappings) = mappings[,1]
mappings[,1] = NULL

for (i in 1:nrow(heatmap_df)) {
  row.names(heatmap_df)[i] = mappings[row.names(heatmap_df)[i], "name"]
}

print("create plot")
diagnosis_str = ""
if(length(unique(new_df[, "Diagnosis"])) == 1) {
  diagnosis_str = paste("for", unique(new_df[, "Diagnosis"]), "samples")
}
print(nchar(diagnosis_str))
print(64*ncol(heatmap_df)/50)
if(nchar(diagnosis_str)>64*ncol(heatmap_df)/50) {
  diagnosis_str = ""
}
pdf(file=args[3], width=max(8, 4+16*ncol(heatmap_df)/50), height=15)
col_fun = colorRamp2(c(0, 0.5, 1), c("beige", "darkorchid3", "darkorchid4"))
h1 = Heatmap(heatmap_df,
  name="Cell type fraction",
  column_title = paste("Cell type deconvolution", diagnosis_str),
  column_title_gp = gpar(fontsize = 24, title_position="lefttop"),
  heatmap_legend_param=list(
    title_gp=gpar(fontsize=18),
    labels_gp=gpar(fontsize=14),
    legend_direction="horizontal",
    legend_width = unit(1.75, "inch")),
  cluster_rows=FALSE,
  show_column_dend=FALSE,
  col=col_fun,
  row_split=c("B cells", "B cells", "B cells", "B cells", "B cells", "B cells", "B cells", "T cells", "T cells", "T cells", "T cells", "T cells", "T cells", "T cells", "T cells", "Myeloid cells", "Myeloid cells", "Myeloid cells", "Myeloid cells", "Myeloid cells", "Myeloid cells"),
  row_names_gp = gpar(fontsize = 16),
  row_title_gp = gpar(fontsize = 24),
  cluster_row_slices=FALSE,
  row_gap = unit(0.3, "inch"))
draw(h1)
dev.off()
