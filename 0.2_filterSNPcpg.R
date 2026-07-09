### Filtering out CpGs that overlay single nucleotide polymorphisms (SNPs) on the genome

setwd("/path/to/project/folder/")

cpgLoci <- read.csv('CpG.bed.gz', sep = "\t", header = FALSE) # all CpG positions on the chromosomes across the genome
snpLoci <- read.csv('0.1_snps_maf_0.01_chr_pos.txt', header=FALSE, sep='_') # SNP positions
cpgLoci_f <- cpgLoci[!(cpgLoci[,2] %in% snpLoci[,2]), ]

write.csv(cpgLoci_f, '0.2_CpG_filt.csv')
