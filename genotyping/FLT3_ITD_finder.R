#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

samp<-args[1]
out_file<-args[2]
out_filt_file<-args[3]

plot_file<-args[4]
table_file<-args[5]



out<-read.csv(out_file, header=F, sep="\t", quote = "")
out_filt<-read.csv(out_filt_file, header=F, sep="\t", quote = "")
  
out$percent_notclip<-out_filt$V4/out$V4
min_notclip<-round(min(out$percent_notclip), digits = 3)
  
pdf(file=plot_file)
g<-ggplot(out, aes(x=V2, y=percent_notclip))+geom_point()+ggtitle(paste0(samp, ": ", min_notclip))+geom_hline(yintercept=min_notclip, linetype=2)+ylim(c(0,1))
print(g)
dev.off()
  
write.table(out, file=table_file, quote = F, row.names = F, sep="\t")
    
  
  
