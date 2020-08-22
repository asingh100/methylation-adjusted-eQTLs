library(dplyr)
library(data.table)
library(methods)
library(plyr)
pc = fread("pc.out")
pc = group_by(pc,gene)
#meth = fread("output_10_peer.txt")
##meth = group_by(meth,gene)
func <- function(xx){
  min_pc=xx$SNP[xx$`p-value`==min(xx$`p-value`)]
  #print("d")
  return(min_pc)
}
min_pc = ddply(pc, .(gene), func)
min_pc$snpg = paste0(min_pc$SNP,"_",min_pc$gene)

meth_pc_snpg = meth[meth$snpg%in%min_pc$snpg,]
meth_pc_snpg$group = 'Common'
meth_pc_snpg$group[meth_pc_snpg$snpg%in%PC_sig$snpg&meth_pc_snpg$snpg%in%Sig_Meth$snpg]="Significant with Both"
meth_pc_snpg$group[(!meth_pc_snpg$snpg%in%PC_sig$snpg)&meth_pc_snpg$snpg%in%Sig_Meth$snpg]="Significant with Only Methylation Adjustment"
meth_pc_snpg$group[meth_pc_snpg$snpg%in%PC_sig$snpg&!meth_pc_snpg$snpg%in%Sig_Meth$snpg]="Significant with Only PC Adjustment"
meth_pc_snpg$group[!meth_pc_snpg$snpg%in%PC_sig$snpg&!meth_pc_snpg$snpg%in%Sig_Meth$snpg]="Not Significant with Both"
pc_meth_snpg = min_pc[min_pc$snpg%in%meth_pc_snpg$snpg,]
pc_meth_snpg = pc_meth_snpg[order(meth_pc_snpg$snpg),]
meth_pc_snpg= merge(meth_pc_snpg,pc_meth_snpg,by="snpg")
#meth_pc_snpg$beta.x=pc_meth_snpg$beta
min_pc_2 = min_pc[!min_pc$snpg%in%meth_pc_snpg$snpg,]
min_pc_2$group = "Only found with PC Adjustment"
min_pc_2$`p-value.x` = min_pc_2$`p-value`
min_pc_2$beta.x = min_pc_2$beta
colnames(min_pc_2)[c(3,5)]=c('beta.y','p-value.y')
min_pc_2 = select(min_pc_2,c('beta.y','beta.x','p-value.y','p-value.x','snpg','group'))
meth_pc_snpg=select(meth_pc_snpg,c('beta.y','beta.x','p-value.y','p-value.x','snpg','group'))
meth_pc_snpg = rbind(meth_pc_snpg,min_pc_2)

meth_pc_snp = meth[meth$SNP%in%min_pc$SNP,]
meth_pc_snp$group = 'Common'
meth_pc_snp$group[meth_pc_snp$snpg%in%PC_sig$snpg&meth_pc_snp$snpg%in%Sig_Meth$snpg]="Significant with Both"
meth_pc_snp$group[(!meth_pc_snp$snpg%in%PC_sig$snpg)&meth_pc_snp$snpg%in%Sig_Meth$snpg]="Significant with Only Methylation Adjustment"
meth_pc_snp$group[meth_pc_snp$snpg%in%PC_sig$snpg&!meth_pc_snp$snpg%in%Sig_Meth$snpg]="Significant with Only PC Adjustment"
meth_pc_snp$group[!meth_pc_snp$snpg%in%PC_sig$snpg&!meth_pc_snp$snpg%in%Sig_Meth$snpg]="Not Significant with Both"
meth_pc_snp= merge(meth_pc_snp,pc_meth_snp,by="SNP")
#meth_pc_snpg$beta.x=pc_meth_snpg$beta
min_pc_2 = min_pc[!min_pc$snpg%in%meth_pc_snp$snpg,]
min_pc_2$group = "Only found with PC Adjustment"
min_pc_2$`p-value.x` = min_pc_2$`p-value`
min_pc_2$beta.x = min_pc_2$beta
colnames(min_pc_2)[c(3,5)]=c('beta.y','p-value.y')
min_pc_2 = select(min_pc_2,c('beta.y','beta.x','p-value.y','p-value.x','SNP','group'))
meth_pc_snp=select(meth_pc_snp,c('beta.y','beta.x','p-value.y','p-value.x','SNP','group'))
meth_pc_snp = rbind(meth_pc_snp,min_pc_2)

plot=ggplot()+geom_point(data=meth_pc_snpg,aes(x=-log10(`p-value.y`),y=-log10(`p-value.x`),color=group),size=4)+scale_color_manual(values=c('khaki3','#40960b','#db70e5','brown4','#0a2045'))+labs(x="-log10(P-value) for PC",y="-log10(P-value) for Methylation",title="P-value Change for Min SNPs for each gene in PC")
plot_2=ggplot()+geom_point(data=meth_pc_snpg,aes(x=beta.y,y=beta.x,color=group,alpha=group,size=group))+scale_size_manual(values=c(3,3,7,7,7))+scale_alpha_manual(values = c(0.3,0.3,1,1,1))+scale_color_manual(values=c('khaki3','#40960b','#db70e5','brown4','#0a2045'))+labs(x="beta for PC",y="beta for Methylation",title="Beta Change for Min SNPs for each gene in PC")
print(plot)
print(plot_2)

plot_3=ggplot()+geom_point(data=meth_pc_snp,aes(x=-log10(`p-value.y`),y=-log10(`p-value.x`),color=group),size=4)+scale_color_manual(values=c('khaki3','#40960b','#db70e5','brown4','#0a2045'))+labs(x="-log10(P-value) for PC",y="-log10(P-value) for Methylation",title="P-value Change for Min SNPs for each gene in PC")
plot_4=ggplot()+geom_point(data=meth_pc_snp,aes(x=beta.y,y=beta.x,color=group,alpha=group,size=group))+scale_size_manual(values=c(3,3,7,7,7))+scale_alpha_manual(values = c(0.3,0.3,1,1,1))+scale_color_manual(values=c('khaki3','#40960b','#db70e5','brown4','#0a2045'))+labs(x="beta for PC",y="beta for Methylation",title="Beta Change for Min SNPs for each gene in PC")
print(plot_3)
print(plot_4)

#fwrite(min_pc,"min_pc.txt",sep=" ")
