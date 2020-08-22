#######################################
#eQTL estimation
#######################################
library(LAMatrix)
library(dplyr)
#library(gsubfn)
library(methods)
useModel = modelLOCAL
library(data.table)

print("file1: dosage, file2: snp location, file3: expression: file4: gene location, file5: covariate, file6: output, file7: local")

args = commandArgs(trailingOnly=TRUE)
SNP_file_name=args[1]
snps_Location_File_name=args[2]

expression_file=args[3]
gene_location_file_name=args[4]

covariates_file_name=args[5]




#covariates_file_name="hgcr_covarite_with3pc.txt"
output_file_name_tra=tempfile()
output_File_name_cis=args[6]
LOCAL_file_name = args[7]
pvOutputThreshold_cis=1
pvOutputThreshold_tra=0
errorCovariance=numeric()
cisDist=1e6


cvrt=SlicedData$new()
table = read.table(covariates_file_name,header = T,check.names = F)
cvrt$CreateFromMatrix(as.matrix(table[,-1]))
print(str(cvrt))
#########################################
#read gene expression
###########################################
gene=SlicedData$new();
gene$fileDelimiter=" "
gene$fileOmitCharacters="NA"
gene$fileSkipColumns=1
gene$fileSkipRows=1
gene$fileSliceSize=2000
gene$LoadFile(expression_file)

#exp = read.table(expression_file)
#test=read.table("test.txt")
#dim(test)
############################################
#read snp
############################################
#snp_file = read.table(SNP_file_name,header = T)
#local_file = read.table(LOCAL_file_name,header = T)
#snp_pos = read.table("snppos.txt",sep = " ",header = T,fill=T)
#row.names(snp_file) = snp_file[,1]
#row.names(snp_pos) = snp_pos[,1]
#snp_file = snp_file[,colnames(local_file)]
#snp_file = snp_file[row.names(local_file),]
#snp_pos$snpid = row.names(local_file)
#snp_file = lapply(snp_file,gsub,pattern="-",replacement=".")
#snp_cols = read.table("snpcol.txt",header=T)
#colnames(snp_file) = colnames(snp_cols)
#print(colnames(local_file))
#print(colnames(snp_file))
#snp_file = select(data = snp_file, -c(MP.2,MP.25,MP.12,MP.67,MP.41,MP.20,MP.54))
#snp_file = snp_file[,colnames(local_file)]
#snp_file = subset(snp_file,row.names(local_file)%in%snp_file)
#row.names(snp_file) = snp_file[,1]
#snp_file = snp_file[match(row.names(snp_file), row.names(local_file)),]
#print(row.names(snp_file))
#snp_file = snp_file[,-1]
#snp_pos = snp_pos[,-1]
#write.table(snp_file,"snp.dosage",quote=F)
#write.table(snp_pos,"snppos.txt",quote=F)
snps = fread(SNP_file_name,sep="\t")
snp=SlicedData$new()
snp$fileDelimiter=" "
snp$fileOmitCharacters="NA"
snp$fileSkipColumns=1
snp$fileSkipRows=1
snp$fileSliceSize=2000
snp$LoadFile(SNP_file_name)
############################################
#read covariate
############################################

#cvrt=SlicedData$new()
#cvrt$fileDelimiter=" "
#cvrt$fileOmitCharacters="NA"
#cvrt$fileSkipColumns=1
#cvrt$fileSkipRows=1
#cvrt$fileSliceSize=2000
#cvrt$LoadFile(covariates_file_name)

############################################
  #read covariate
############################################
   
     
locals = SlicedData$new();
locals$fileDelimiter = " ";      # the TAB character
locals$fileOmitCharacters = "NA"; # denote missing values;
locals$fileSkipRows = 1;          # one row of column labels
locals$fileSkipColumns = 1;       # one column of row labels
locals$fileSliceSize = 2000;      # read file in slices of 2,000 rows
locals$LoadFile(LOCAL_file_name);

#############################################\
#read snp location and gene position
#################################################


snpspos=read.table(snps_Location_File_name,header=T,stringsAsFactors=F,fill=T)
genepos=read.table(gene_location_file_name,header = T,stringsAsFactors = FALSE)
head(snpspos)
head(genepos)
colnames(genepos) = c("geneid","chr","left","right")
genepos = na.omit(genepos)
print(nrow(genepos))
#write.table(snpspos,snps_Location_File_name,col.names = FALSE,row.names = FALSE,quote = FALSE)
print(ncol(snp))
print(ncol(gene))
#############################################
#compute eQTL
#############################################
me=LAMatrix_main(snps=snp,
                    gene=gene,
                    cvrt=cvrt,
                    local=locals,
                    output_file_name = output_file_name_tra,
                    pvOutputThreshold = pvOutputThreshold_tra,
                    useModel = useModel,
                    errorCovariance = errorCovariance,
                    verbose = TRUE,
                    output_file_name.cis = output_File_name_cis,
                    pvOutputThreshold.cis = pvOutputThreshold_cis,
                    snpspos = snpspos,
                    genepos = genepos,
                    cisDist = cisDist,
                    pvalue.hist = "qqplot",
                    min.pv.by.genesnp = FALSE,
                    noFDRsaveMemory = FALSE)

#cis_eqtl_output = "all_cis_eqtl_matrixEL_local_gtex_muscle.txt"
#write.table(me$cis$eqtls,cis_eqtl_output,row.names = FALSE,col.names = FALSE,quote = FALSE)
#trans_eqtl_output = "trans_eqt.txt"
#write.table(me$tra$eqtls,trans_eqtl_output,row.names = FALSE,col.names = FALSE,quote = FALSE)
#write.table(me,"pvalues_5_peer.txt")
png("all_eqtl_qq_cis_matrixEL_eig_10_peer.png",res=300,width = 17.35, height = 23.35, units = "cm")
plot(me)
#dev.off()
