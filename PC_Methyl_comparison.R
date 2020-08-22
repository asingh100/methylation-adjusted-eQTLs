library(dplyr)
library(ggplot2)
Full = fread("pc_meth_plot_2.5kb.txt")
Full = filter(Full,group == "common")
new = fread("eqtl_inc_notsigPC.txt")
inc = fread("eqtl_inc_sigPC.txt")
dec = fread("eqtl_dec_sigPC.txt")
new_pc = filter(Full,snpg%in%new$snpg)
new_pc$group = "Significant Only In Methylation"
dec_pc = filter(Full,snpg%in%dec$snpg)
dec_pc$group = "Decreased Significance in Methylation"
inc_pc = filter(Full,snpg%in%inc$snpg)
inc_pc$group = "Increased Significance in Methylation"
combined_pc = rbind(new_pc,dec_pc)
combined_pc = rbind(combined_pc,inc_pc)
plot_1 = ggplot(data=combined_pc,aes(x=-log10(V5.x),y=-log10(V5.y),color=group))+scale_color_manual(values=c("black","brown3","darkkhaki"))+geom_point(size=5)+labs(x="-log10(p-value) for PC Adjusted eQTLS",y="-log10(p-value) for Methylation Adjusted eQTLS",title="p-value comparison for PC vs Methylation adjustment")
print(plot_1)

#Full_Meth = fread("output_10_peer_2.5kb.txt")
#Full_Meth = filter(Full_Meth,FDR<0.05)
Sig_Meth = filter(meth,FDR<0.05)
fwrite(Sig_Meth,"output_10_peer_FDR_0.05.txt",sep=" ")
condition1 <- read.csv("condition1_expression_tmm_normalization_exp_v7_10peer_all_dosage_BY_BH_sig_0.05.out", sep="", stringsAsFactors=FALSE)
condition1$snpg = paste0(condition1$SNP,"_",condition1$gene)
common = data.frame(snpg=union(Sig_Meth$snpg,condition1$snpg))
meth_unique = Sig_Meth$snpg[!Sig_Meth$snpg%in%condition1$snpg]
pc_unique = condition1$snpg[!condition1$snpg%in%Sig_Meth$snpg]
meth_unique=meth_unique[meth_unique%in%eqtl_new$snpg]
meth_common = merge(common,Sig_Meth,by="snpg",all=T)
pc_common = merge(common,condition1,by="snpg",all=T)
table = merge(pc_common,meth_common,by="snpg")
table$group = 'common'
table$group[table$snpg%in%meth_unique] = "Methylation Unique"
table$group[table$snpg%in%pc_unique] = "PC Unique"
table = select(table,"p.value","p-value","snpg","group")
colnames(table) = c("V5.x","V5.y","snpg","group")
pc_n = table[is.na(table$V5.x),]
pc_n_1 = filter(pc,snpg%in%pc_n$snpg)
pc_n = merge(pc_n,pc_n_1,by="snpg")
pc_n$V5.x = pc_n$`p-value`
table$V5.x[table$group=="Methylation Unique"] = pc_n$V5.x[pc_n$group =="Methylation Unique"]
meth_n = table[is.na(table$V5.y),]
meth_n_1 = filter(meth,snpg%in%meth_n$snpg)
meth_n = merge(meth_n,meth_n_1,by="snpg")
meth_n$V5.y = meth_n$`p-value`
table$V5.y[table$group=="PC Unique"] = meth_n$V5.y[meth_n$group =="PC Unique"]
table_2=na.omit(table)

method = lm(-log10(V5.y)~-log10(V5.x),table_2,group=="common")

library(ggplot2)
png("pc_meth_10peer_sig.png")
plot=ggplot()+geom_point(data=table_2, aes(x=-log10(V5.x), y=-log10(V5.y), color=group),size=5)+scale_color_manual(values=c("black","brown4","khaki3"))+labs(x="-log10(p) for PC", y =" -log10(p) for Methylation", title = "PC eQTLs vs Methylation eQTLs for 10 peers")
print(plot)
