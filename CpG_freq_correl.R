library(dplyr)
library(data.table)
library(ggplot2)
library(plyr)
#library(stats)
library(methods)
library(tidyr)
combined = fread("overlapping_with_chr_pos.txt",fill=T)
#combined = select(combined,-c("chrCHR"))
SNP_Gene = fread("SNP_gene_pairs.txt")
combined = filter(combined,SNP%in%SNP_Gene$SNP)
#CpG_number_per_SNP = as.data.frame(table(combined$SNP))
#head(CpG_number_per_SNP)
#for_plot = as.data.frame(table(CpG_number_per_SNP$Freq))
#for_plot$Var1 = as.numeric(for_plot$Var1)
#for_plot$Freq = as.numeric(for_plot$Freq)
#head(for_plot)
#for_plot = fread("for_plot2.txt")
#weight = (for_plot$Freq)/sum(for_plot$Freq)
#for_plot$w = rep(w,nrow(for_plot))
#w = weighted.mean(for_plot$Var1,weight)
#pdf("CpG_per_SNP_plot.pdf")
#plot = ggplot(for_plot,aes(x=Var1,y = Freq)) + geom_bar(stat="identity",color = "black",fill="white") + labs(title = "Distribution of number of CpG Sites", x = "Number of CpG Sites", y = "Number of SNPs")+geom_vline(aes(xintercept = w))+ geom_text(aes(x=4, label="Mean number of CpG per SNP", y=9e5), colour="black", angle=90, vjust = 1.2, size =5)+geom_label(aes(x=w, label=paste("X = ", round(w,2)), y=5.7e5))+geom_label(aes(label=Freq))
#print(plot)
#dev.off()
#write.table(for_plot,"for_plot2.txt",quote=F,row.names=F)
#combined = fread("overlapping_with_chr_pos.txt",fill=T)
head(combined)
combined = select(combined,c("SNP","dist","gene"))
combined = group_by(combined,"SNP")
func <- function(xx){
  xx = data.frame(gene=xx$gene,d=xx$dist)
  row.names(xx) = xx$gene
  xx = data.frame(d = xx$d,row.names=row.names(xx))
if(nrow(xx) ==1){
corr = NA
}
 else{
   for(i in 1:nrow(xx)){
     xx = cbind(xx,sample(xx$d,nrow(xx)))
     corr = cor(xx)
   }
 }
 corr[corr==1]<-NA
 return(mean(corr,na.rm=T))
 }
 x = ddply(combined, .(SNP), func)
 
  head(x)
# x[is.na(x)] <- 1
 mean(x$V1,na.rm=T)
 write.table(x,"correlations_gene.txt",quote=F,row.names=F)

# dist = data.frame(gene=combined$gene,d=combined$dist)
# dist = dist[4:6,]
# row.names(dist) = dist$gene
# dist = data.frame(d = dist$d,row.names=row.names(dist))
# for(i in 2:nrow(dist)){
#   dist = cbind(dist,sample(dist$d,3))
# }
# mean(cor(dist))

