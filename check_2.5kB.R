library(dplyr)
library(data.table)
library(methods)
kb = fread("matching_positions_unedited_2.txt.txt")
pos = fread("yz_hepatocyte_snp.position")
meth = fread("MethylationEPIC2.txt")
colnames(meth)[1] = "gene"
colnames(pos)[1] = "SNP"
#pos$pos_2500 = pos$pos +2500
#pos$neg_2500 = pos$pos - 2500
combined = merge(meth,kb,by="gene",sort = T)
combined = merge(pos,combined,by="SNP",sort=T)
combined = select(combined,-c("beta","t-stat","p-value","FDR","END"))
#write.table(combined,"overlapping_with_chr_pos.txt",quote=F,row.names=F)
#combined = fread("overlapping_with_chr_pos.txt")
#combined = select(combined,-c("chrCHR"))
#combined_2 = data.frame(SNP=combined$SNP,gene = combined$gene)
combined$dist= abs(combined$MAPINFO - combined$pos)
combined = filter(combined,dist<=2500)
write.table(combined,"overlapping_with_chr_pos_2.txt",quote=F,row.names=F)
#write.table(combined_2,"overlapping.txt",quote=F,row.names=F)
#overlap = combined$dist
#index= which(overlap>2500)
#write.table(cbind(index,overlap[index]),"check_pos_values.txt",quote=F)

