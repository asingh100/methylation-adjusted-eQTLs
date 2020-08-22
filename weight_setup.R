#setwd("/home/asi8235/Documents/2nd_Iteration")
library(dplyr)
library(data.table)
library(methods)
library(plyr)
#pos = fread("yz_hepatocyte_snp.position",nrow=30000)
meth = fread("csp_methnorm_Combat_final_53.txt",fill=T)
combined = fread("overlapping_with_chr_pos_1.txt",fill=T)
colnames(meth)[1]="gene"
combined = merge(combined,meth,by="gene")
combined = select(combined,-c("chr","pos","pos_2500","neg_2500","MAPINFO"))
combined = group_by(combined,SNP)
combined$dist = abs(combined$dist)
func <- function(xx){
  if(nrow(xx) ==1){
    xy = select(xx,-c("gene","SNP"))
    dist = 1-(xy$dist/2500)
    xy = select(xy,-c("dist"))
    methyl = dist*xy[,1:53]
  }
  else{
      xy = select(xx,-c("gene","SNP"))
      dist = 1-(xy$dist/2500)
      dist = dist/(sum(dist))
      xy = select(xy,-c("dist"))
      i = 1
      methyl=xy[1:length(dist),1:53]
      while(i<=length(dist)){
      methylation = data.frame(dist[i]*xy[,1:53])
      methyl[i,1:53] = methylation
      i = i+1
      }
      methyl_sum = colSums(methyl)
      methyl= rbind(methyl,methyl_sum)
      methyl = last(methyl)
      }
  print("d")
  return(methyl)
}
x = ddply(combined, .(SNP), func)
count(x>=1)
count(x<=0)
write.table(x,"SNP_Methylation_2_1.txt",row.names=F,quote=F)
#x = x[,-c(1:54)]
# xy = combined[4:6,4:56]
# dist = (sum(combined$dist[4:6])-abs(combined$dist[4:6]))/sum(combined$dist[4:6])
# xy = select(xy,-c("SNP","gene","dist"))
# i=4
# methyl= xy[4:6,1:53]
# while(i<=6){
# methylation = data.frame(dist[i]*xy[i,1:53])
# methyl[i-3,1:53] = methylation
# i = i+1
# }
# methyl_sum= colSums(methyl)
# methyl = rbind(methyl,methyl_sum)
# methyl=data.frame(methyl[-c(1:length()),])
#meth = fread("SNP_Methylation.txt")
