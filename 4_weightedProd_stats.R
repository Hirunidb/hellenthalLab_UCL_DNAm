### Analysis of the model performance
### 1: Classification accuracy of the sample
### 2: Change in accuracy across the genome

setwd("/path/to/project/folder")

library(reshape2)
library(ggplot2)
library(dplyr)

# repeated for;
## I.  both 6CT and 2CT classification
## II. 10kb, 15kb, 20kb
eq2_subCT1 <- readRDS('output/3_weightedProdEQ2/weightedProd_6CT_p49_10kb.rds')
eq2_subCT2 <- readRDS('output/3_weightedProdEQ2/weightedProd_6CT_p51_10kb.rds')
eq2_subCT3 <- readRDS('output/3_weightedProdEQ2/weightedProd_6CT_p54_10kb.rds')

resultList <- list(eq2_subCT1, eq2_subCT2, eq2_subCT3)
ctLabels <- colnames(eq2_subCT1)
targetCT <- stringr::str_split_i(ctLabels, '_', 1)
#targetCT <-  c(rep('Lymphoid',12), rep('Myeloid', 17), 'Lymphoid')
targetCTu <- unique(targetCT)
pId <- c('49','51','54')

clsAccuracy_sample_6CT <- list()
clsAccuracy_genome_6CT <- list()
z <- 1
y <- 1


for (id in seq(length(resultList))) {
  message(paste0('Running_', pId[id]))
  
  # TP + FP  
  eq2_subCT <- resultList[[id]]
  maxInd_max <- apply(eq2_subCT, c(2,3), max)
  maxMask <- maxInd_max < (1/length(targetCTu)) # filter out any probabilities below random (all are 0 columns)
  
  maxInd_all <- apply(eq2_subCT, c(2,3), which.max)
  maxInd_all[maxMask] <- NA
  
  # ground truth (TP+FN)
  nReadsNA <- length(maxInd_all) - sum(maxMask)
  GT_all <- matrix(match(targetCT, targetCTu), nrow=dim(eq2_subCT)[2], ncol = dim(eq2_subCT)[3])
  GT_all[maxMask] <- NA
  
  # ratio between the probability of best and second best alternative
  maxRatio <- apply(eq2_subCT, c(2,3), function(x) {
    vals <- sort(x, decreasing = T)[1:2]
    if (any(vals==0)){
      vals <- vals+1e-10
    }
    vals[2]/vals[1]
  })
  maxRatio[maxMask] <- NA
  
  
  threshold <- seq(0.1, 1, 0.1)
  
  # 1: Classification accuracy - confusion matrix
  for (ct in targetCTu){
    message(ct)
    threshAcc_LM <- matrix(0, length(threshold), 8,
                           dimnames =list(threshold, c('pId', 'Class', 'Threshold', 'TP', 'FP', 'TN', 'FN', 'Dropout')))
    threshAcc_LM[,'Class'] <- ct
    ind <- which(targetCTu == ct)
    threshAcc_LM[,'pId'] <- pId[id]
    
     
    for (i in seq(length(threshold))) {
      maxInd <- maxInd_all
      GT <- GT_all
      maxInd[maxRatio > threshold[i]] <- NA
      GT[maxRatio > threshold[i]] <- NA
      threshAcc_LM[i,'Threshold'] <- threshold[i]
      predict_positive <- maxInd == ind
      TP <- sum(predict_positive & (GT==ind), na.rm=T)
      FP <- sum(predict_positive, na.rm=T) - TP
      predict_negative <- maxInd != ind
      FN <- sum(predict_negative & (GT==ind), na.rm=T)
      TN <- sum(predict_negative, na.rm=T) - FN
      
      threshAcc_LM[i, 'TP'] <- TP
      threshAcc_LM[i, 'FP'] <- FP
      threshAcc_LM[i, 'TN'] <- TN
      threshAcc_LM[i, 'FN'] <- FN
      threshAcc_LM[i, 'Dropout'] <- sum(is.na(maxInd))
    }
    clsAccuracy_sample_6CT[[z]] <- threshAcc_LM
    z <- z+1
  }


  # 2: Accuracy across the genome
  ### overall
  acc_threshold <- seq(0, 1, 0.1)
  threshAcc_LM <- matrix(0, length(threshold), length(acc_threshold)+5,
                         dimnames =list(threshold, c('pId', 'Class', 'Threshold', 'Loss', 'Total_nslices_left', 'Null', acc_threshold[-11])))
  threshAcc_LM[,'pId'] <- pId[id]
  threshAcc_LM[,'Class'] <- 'All'

  for (i in seq(length(threshold))){
    i <- 10
    threshAcc_LM[i,'Threshold'] <- threshold[i]
    a <- maxRatio <= threshold[i]
    gt <- apply(a, 2, sum, na.rm=T)
    # 0 - means none of the reads from the region passed the threshold, hence to be removed
    b <- maxInd == GT

    null_slices <- which(gt == 0)
    threshAcc_LM[i,'Loss'] <- length(null_slices)
    threshAcc_LM[i,'Total_nslices_left'] <- length(gt) - length(null_slices)
    if (length(null_slices)>0) {
      a <- a[,-null_slices]
      gt <- gt[-null_slices]
      b <- b[,-null_slices]
    }
    tp <- apply(a&b, 2, sum, na.rm=T)
    tp_prop <- tp/gt

    threshAcc_LM[i,'Null'] <- sum(tp_prop == 0, na.rm=T)
    for (j in seq(length(acc_threshold)-1)){
      con <- (tp_prop > acc_threshold[j]) & (tp_prop <= acc_threshold[j+1])
      threshAcc_LM[i, j+6] <- sum(con, na.rm=T)
    }
  }
  clsAccuracy_genome_6CT[[y]] <- threshAcc_LM
  y <- y+1

  ### CT specific
  for (ct in targetCTu){
    message(ct)
    threshAcc_LM <- matrix(0, length(threshold), length(acc_threshold)+5,
                           dimnames =list(threshold, c('pId', 'Class', 'Threshold', 'Loss', 'Total_nslices_left', 'Null', acc_threshold[-11])))
    threshAcc_LM[,'pId'] <- pId[id]
    threshAcc_LM[,'Class'] <- ct
    ind <- which(targetCT==ct)

    for (i in seq(length(threshold))){
      i <- 10
      threshAcc_LM[i,'Threshold'] <- threshold[i]
      a <- maxRatio[ind, ,drop=F] <= threshold[i]
      gt <- apply(a, 2, sum, na.rm=T)
      b <- maxInd[ind, , drop=F] == GT[ind, , drop=F]

      null_slices <- which(gt == 0)
      threshAcc_LM[i,'Loss'] <- length(null_slices)
      threshAcc_LM[i,'Total_nslices_left'] <- length(gt) - length(null_slices)
      if (length(null_slices)>0) {
        a <- a[,-null_slices]
        gt <- gt[-null_slices]
        b <- b[,-null_slices]
      }
      if (ct=='B'){
        tp_prop <- gt[a&b] / gt[a&b]
      } else{
        tp_prop <- apply(a&b, 2, sum, na.rm=T) / gt
      }

      threshAcc_LM[i,'Null'] <- sum(tp_prop == 0, na.rm=T)
      for (j in seq(length(acc_threshold)-1)){
        con <- (tp_prop > acc_threshold[j]) & (tp_prop <= acc_threshold[j+1])
        threshAcc_LM[i, j+6] <- sum(con, na.rm=T)
      }
    }
    clsAccuracy_genome_6CT[[y]] <- threshAcc_LM
    y <- y+1
  }

}


clsAccuracy_sample_6CT_CM <- do.call(rbind, clsAccuracy_sample_6CT)
clsAccuracy_genome_6CT_all <- do.call(rbind, clsAccuracy_genome_6CT)

saveRDS(clsAccuracy_sample_6CT_CM, 'output/4_weightedProdEQ2_6CT_confMat_10kb.rds')
saveRDS(clsAccuracy_genome_6CT_all, 'output/4_weightedProdEQ2_genome_6CT_stats_10kb.rds')
