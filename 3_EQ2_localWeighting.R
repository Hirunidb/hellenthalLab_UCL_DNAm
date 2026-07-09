
setwd("/path/to/project/folder")

# reference simulation
refWeights_norm <- 'output/1_referenceProdArr_prop_norm.rds'
refWeights_norm_15kb <- 'output/1_referenceProdArr_prop_norm_15kb.rds'
refWeights_norm_20kb <- 'output/1_referenceProdArr_prop_norm_20kb.rds'

# target data simulation
## 3 patients, 3 read lengths
prodArr49_10 <- 'output/2_prodArr_p49.rds'
prodArr51_10 <- 'output/2_prodArr_p51.rds'
prodArr54_10 <- 'output/2_prodArr_p54.rds'

prodArr49_15 <- 'output/2_prodArr_p49_15kb.rds'
prodArr51_15 <- 'output/2_prodArr_p51_15kb.rds'
prodArr54_15 <- 'output/2_prodArr_p54_15kb.rds'

prodArr49_20 <- 'output/2_prodArr_p49_20kb.rds'
prodArr51_20 <- 'output/2_prodArr_p51_20kb.rds'
prodArr54_20 <- 'output/2_prodArr_p54_20kb.rds'

kb10 <- c(prodArr49_10, prodArr51_10, prodArr54_10)
kb15 <- c(prodArr49_15, prodArr51_15, prodArr54_15)
kb20 <- c(prodArr49_20, prodArr51_20, prodArr54_20)
test <- list(kb10, kb15, kb20)
ref <- c(refWeights_norm, refWeights_norm_15kb, refWeights_norm_20kb)
strLens <- c('10', '15', '20')
pId <- c('49', '51', '54')

# classifications
colanno <- read.csv('metadata/colanno27.csv', row.names=1)
sctype_ref_mod <- c('T-CD4', 'T-CD8', 'T-CD4', unique(colanno[, 'stype'])[4:9])
sctype_ref1 <- rep(sctype_ref_mod, each=3)
sctype_ref2 <- c(rep('Lymphoid', 18), rep('Myeloid', 6), rep('Lymphoid', 3))
clsType <- list(sctype_ref1, sctype_ref2)
clsTypeName <- c('6CT', 'LymMyl')

summaryDF <- data.frame(
  'classification_type' = character(),
  'fileName'= character(),
  'patientID'=character(),
  'window_length_kb' = integer(),
  'total_nWindows' = integer(),
  'prior_null_nWindows'= integer(),
  'final_nWindows'= integer()
)

for (cls in seq(length(clsType))) {
  sctype_ref <- unlist(clsType[cls])
  scts_ref <- unique(sctype_ref)
  
  for (i in seq(length(ref))) {
    message(paste0('Running length', i, '_', clsTypeName[cls]))
    ###### EQ2.1
    refWeights <- readRDS(ref[i])
    refPriorWeights <- array(0, c(dim(refWeights)[2], length(scts_ref), dim(refWeights)[3]), 
                             dimnames=list(NULL, scts_ref, NULL))
    
    for (j in seq(dim(refWeights)[3])) {
      refMat <- refWeights[,,j]
      weightMat <- matrix(0, ncol(refMat), length(scts_ref), dimnames = list(NULL, scts_ref))
      
      # weight calculation per region, per ref, per CT 
      # NOTE: there are NO NAs left after this step
      for (k in seq(ncol(refMat))){
        for (ct in scts_ref) {
          weightMat[k,ct] <- mean(refMat[sctype_ref==ct, k], na.rm=TRUE)
        }
      }
      
      rowsums <- rowSums(weightMat, na.rm = TRUE)
      weightMat_norm <- weightMat
      nonzero <- rowsums != 0
      weightMat_norm[nonzero,] <- weightMat[nonzero,] / rowsums[nonzero]
      # 
      colsums <- colSums(weightMat_norm, na.rm=TRUE)
      weightMat <- weightMat_norm
      nonzero <- colsums != 0
      weightMat[,nonzero] <- weightMat_norm[,nonzero] / rep(colsums[nonzero], each=nrow(weightMat))
      
      refPriorWeights[,,j] <- weightMat
    }
    
    info <- apply(refPriorWeights, 3, colSums)
    ind <- which(apply(info, 2, function(x) all(x==0)))
    refPriorWeights <- refPriorWeights[,,-ind] # 190 removed
    saveRDS(ind, paste0('output/3_allZeroRemovedSliceInd_', strLens[i], 'kb.rds'))
    
    ######EQ2
    testArrs <- unlist(test[i])
    for (m in seq(length(testArrs))) {
      message(paste('Running', pId[m], strLens[i], sep='_'))
      prodArr <- readRDS(testArrs[m])
      prodArr <- prodArr[,,-ind]
      
      ref_groups <- colnames(refPriorWeights)
      samples <- colnames(prodArr)
      eq2_subCT <- array(0, c(length(ref_groups), dim(prodArr)[2:3]), dimnames=list(ref_groups, samples, NULL))
      
      for (n in seq(dim(prodArr)[3])) {
        weightMat <- refPriorWeights[,,n]
        testMat <- prodArr[,,n]
        
        for (z in seq(ncol(testMat))) {
          vals <- colSums(testMat[,z] * weightMat, na.rm = TRUE) # EQ2
          s <- sum(vals, na.rm=T)
          if (s>0) {
            eq2_subCT[,z,n] <- vals / s
          }
        }
      }
      filename <- paste0('output/3_weightedProdEQ2/weightedProd_', clsTypeName[cls],
      '_p', pId[m], '_', strLens[i], 'kb.rds')
      message(paste0('Saving_', filename))
      saveRDS(eq2_subCT, filename)
      
      
      summary <- c(clsTypeName[cls], stringr::str_split_i(basename(filename), '\\.', 1), pId[m], strLens[i],
        dim(refWeights)[3], length(ind), dim(refPriorWeights)[3])
      summaryDF <- rbind(summaryDF, summary)
    }
  }

}

write.csv(summaryDF, 'output/3_weightedProdEQ2/weightedProd_metadata.csv')
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
