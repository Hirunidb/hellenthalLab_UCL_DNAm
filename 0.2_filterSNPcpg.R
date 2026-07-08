#setwd("/SAN/ghlab/epigen/Hiruni")
setwd("/Users/hirunidb/Library/CloudStorage/OneDrive-Personal/Projects/Hellenthal-SilverLab_LIDo/hellenthalSilverLab_DNAm_pvt/")

cpgLoci <- read.csv('hannon-mill/loyfer_WGBS_rawData/CpG.bed.gz', sep = "\t", header = FALSE)
snpLoci <- read.csv('metadata/0.1_snps_maf_0.01_chr_pos.txt', header=FALSE, sep='_') # SNP positions
cpgLoci_f <- cpgLoci[!(cpgLoci[,2] %in% snpLoci[,2]), ]

write.csv(cpgLoci_f, 'metadata/0.2_CpG_filt.csv')
