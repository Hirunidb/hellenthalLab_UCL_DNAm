### Simulation of methylation data of a blood sample from one healthy donor
### using methylation atlas data from Looyfer et.al., 2023

setwd("/path/to/project/folder/")

cpgLoci <- read.csv('metadata/0.2_CpG_filt.csv', row.names = 1) # snp filtered
betaMat <- read.csv('output/0_loyfer_methylProp27.csv', row.names=1) # meth proportions of 27 sample

# cell type composition of blood sample
colanno <- read.csv('metadata/colanno27.csv', row.names=1)
patient <- 49
patientData <- colanno[, 'PatientID']
ctype <- colanno[, 'stype']
ctProb <- c(0.12, 0.13, 0.09, 0.12, 0.5, 0.04) # in the order:  T, NK, Monocytes, Granulocytes, B 
readDepth <- 30
set.seed(15)
readFrac <- rmultinom(n=1, size=readDepth, prob = ctProb)
ctLabels <- paste0(rep(ctFracMat$cellType, readFrac), '_', seq(readDepth))
ref_cols <- paste0(ctype, '_', seq_along(length(ctype)))

strLen <- 1e4
chr <- paste0('chr', seq.int(1,22))
const <- 1e-10

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
  mseList <- lapply(splitInt, function(rows){
    if (length(rows) <= 10) return(NULL)
    refmat <- as.matrix(mat[rows, , drop=FALSE])
    
    # test data
    patientMat <- refmat[, patientData==patient] # subset where all samples come from one individual
    
    ## simulate the target reads per cell types
    testMat <- do.call(cbind, lapply(seq(ncol(patientMat)), function(i) {
      nReads <- readFrac[i]
      matrix(rbinom(nReads*length(rows), 1, patientMat[,i]), nrow=length(rows))
    }))
    colnames(testMat) <- ctLabels
    

    # reference data - proportions
    refmat[refmat == 0] <- const # ref data
    refmat[refmat == 1] <- 1-const
    A <- log(refmat)
    a <- log(1 - refmat)

    
    # comparison
    res <- matrix(0, nrow=ncol(refmat), ncol=ncol(testMat), dimnames = list(ref_cols, ctLabels)) 
    for (x in seq_len(ncol(testMat))) {
      # test data
      D <- testMat[,x]
      d <- 1 - D
      
      prods <- colSums((A*D) + (a*d)) ## prod( (A^ D) (1-A ^ 1-D) )
      prods <- exp(prods - max(prods))
      prods[patientData == patient] <- NA # remove same indv
      s <- sum(prods, na.rm = TRUE) # normalize, colSums == 1
      if (s>0) {
        res[, x] <- prods/s # if s==0 -> Nan
      }
    }
    res
  })
  Filter(Negate(is.null), mseList)
})

prod_allchrArr <- simplify2array(unlist(prod_allchr, recursive=FALSE))
saveRDS(prod_allchrArr, 'output/2_prodArr_p49.rds')

write.csv(data.frame(
  'cellType'= ctype[patientData==patient], 
  'ctFrac'=ctProb, 
  'readFrac'=readFrac
), 'metadata/2_targetSampleComposition.csv')



