#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

samp<-args[1]
out_file<-args[2]
out_filt_file<-args[3]



out<-read.csv(out_file, header=F, sep="\t", quote = "")
out_filt<-read.csv(out_filt_file, header=F, sep="\t", quote = "")
  
out$percent_notclip<-out_filt$V4/out$V4
min_notclip<-round(min(out$percent_notclip), digits = 3)
  
pdf(file=paste0(samp, "_plot.pdf"))
g<-ggplot(out, aes(x=V2, y=percent_notclip))+geom_point()+ggtitle(paste0(samp, ": ", min_notclip))+geom_hline(yintercept=min_notclip, linetype=2)+ylim(c(0,1))
print(g)
dev.off()
  
write.table(out, paste0(samp, "_flt3_itd_parsed_vals.txt"), quote = F, row.names = F, sep="\t")
  
  
  
