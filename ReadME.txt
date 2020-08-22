Instructions:
- Import dosage and genotype information into MatrixeQTL as the program describes and replace the gene file positions with CpG site information (for gene expression put in methylation values for CpG sites and for gene position put in CpG site position). Change the cis distance to 2.5E3 to measure CpG sites within a 2.5kB window of a SNP, you can change this depending on what window size you are using.
- Use the check_2.5kB.R program to eliminate all CpG sites that are not in the 2.5kB window as matrixeQTL gives some sites that are a little bit outside of the window
- Use CpG_freq_correl.R program to check the distribution of CpG sites per SNP in your data and to check the average correlation between CpG sites attached to one SNP
- Use weight_setup.R to compute the weighted average of all the CpG sites attached to one SNP to get the methylation value for that SNP to be used in LAMatrix
- Use 2ndIteration.R to run LAMatrix to find the methylation-adjusted eQTLs in your data, use your methylation file from weight_setup.R in the LOCAL_file's place. 
- Finally, use PC_Methyl_comparison.R to compare the results of PC-adjusted eQTL mapping with methylation-adjusted eQTL mapping 

Files and Descriptions:

- 2ndIteration.r:
	- LAMatrix Program to run the methylation-adjusted eQTL analysis

- CpG_freq_correl.R:
	- program to check how many CpG sites there are per SNP and the average correlation between CpG sites attached to one SNP

- PC_Methyl_comparison.R:
	- program to analyze the difference between PC-adjusted and Methylation-adjusted eQTLs

- check_2.5kB.R:
	-program that checks if CpG sites are within 2.5kB of the SNP that Matrix eQTL linked it to and deletes all SNPs that are greater than 2.5 kB away in distance

-esnp_comparison.R:
	- program to compare the effect sizes between methylation-adjusted and PC-adjusted eQTLs

- weight_setup.R:
	- program to gather all of the methylation sites within 2.5kB of a SNP and compute the weighted average of these sites to get methylation as a SNP-based covariate