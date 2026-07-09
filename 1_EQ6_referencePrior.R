### reference (proportion) compared against the references (proportion) 
### i.e, how likely is a reference read to come from a certain cell type, determined using the other available reference reads from the given cell type and region
### A measure of how well the cell types can be distinguished within the reference reads

setwd("/path/to/project/folder/")

cpgLoci <- read.csv('metadata/0.2_CpG_filt.csv', row.names = 1) # SNP filtered
betaMat <- read.csv('output/0_loyfer_methylProp27.csv', row.names=1) # methylation proportions of 27 sample
strLen <- 1e4 # Nanopore read length
chr <- paste0('chr', seq.int(1,22))

prod_allchr <- lapply(chr ,function(n){
  message(paste0('Running ', n))
  # separate chromosome
  idx <- which(cpgLoci[,1] == n) 
  mat <- betaMat[idx,]
  loci <- cpgLoci[idx,]
  
  # remove missing CpG
  message('Filtering NA')
  idx <- which(rowSums(is.na(mat)) == 0)
  mat <- mat[idx,]
  loci <- loci[idx,]
  
  # subset 10kb sections
  message('Position calculation')
  strEnd <- loci[nrow(loci), 2]
  start <- seq.int(loci[1,2], strEnd, strLen) # all start positions for 10kb sections
  start <- start[start <= strEnd+1 - strLen] # removing start position if the length is less that 10kb
  allPos <- loci[loci[,2] < start[length(start)]+strLen, 2] # relevant cpg positions
  interval <- findInterval(allPos, start) 
  splitInt <- split(seq_along(interval), interval) # 10kb separations on cpg loci
  
  message('10kb window calculation')
  prodList <- lapply(splitInt, function(rows){
    if (length(rows) <= 10) return(NULL)
    refmat <- as.matrix(mat[rows, , drop=FALSE]) # 10kb proportions

    res <- matrix(0, nrow=ncol(refmat), ncol=ncol(refmat)) 
    for (x in seq_len(ncol(refmat))) {
      D <- refmat[,x]
      prods <- colSums(dnorm(D, refmat, 0.05, log = TRUE))
      
      prods <- exp(prods - max(prods))
      prods[patient == patient[x]] <- NA # remove same indv
      s <- sum(prods, na.rm=TRUE) # normalize, colSums == 1
      if (s>0) {
        res[, x] <- prods/s # if s==0 -> Nan
      }
    }
    res
  })
  Filter(Negate(is.null), prodList)
})

refPriorWeights <- simplify2array(unlist(prod_allchr, recursive=FALSE))

saveRDS(prodArr, 'output/7.1_referenceProdArr_prop_norm.rds')



