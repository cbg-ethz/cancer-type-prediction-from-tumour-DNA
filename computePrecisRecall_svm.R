rm(list=ls())
setwd("~/Documents/Code/ctDNAdiagnostic")

library(LiblineaR)
library(parallel)
#library(snow)       # package for parallel execution on Windows OS
#library(snowfall)   # package for parallel execution on Windows OS

load("mcmat_expBnrCNA2017.RData")       # gene alteration data matrix
load("testIdx_global.RData")            # a set of indices as the final unseen test set
load("testIdx_50bootstrapSets.RData")   # a list of 50 sets of indices for test cases
load("geneRank_mcmat_svmrfe.RData")     # gene ranking results from feature selection via SVM-RFE

####### Function to compute class precision, recall and accuracy #######
getPrecisionRecall <- function(u, indexSets, testIdx.global, dataMat, rankingList) {
  allTestCases <- c(indexSets[[u]], testIdx.global)
  lastColIdx <- ncol(dataMat)
  labels <- dataMat[,lastColIdx]
  x.train <- dataMat[-allTestCases,-lastColIdx]
  x.test <- dataMat[indexSets[[u]],-lastColIdx]
  y.train <- labels[-allTestCases]
  y.test <- labels[indexSets[[u]]]
  
  featureRankedList <- rankingList[[u]]
  selectedFeatures <- featureRankedList[1:900]   # number of genes needed to give highest accuracy
    
  # find the best parameter for svm
  tryCosts <- c(0.01,0.02,0.03)
  bestCost <- NA
  bestAcc <- 0
  for (ct in tryCosts) {
    training.acc <- LiblineaR(data=x.train[,selectedFeatures], target=as.factor(y.train), 
                              type=1, cost=ct, bias=TRUE, cross=10)
    if (training.acc > bestAcc) {
      bestCost <- ct
      bestAcc <- training.acc
    }
  }
    
  bestModel <- LiblineaR(data=x.train[,selectedFeatures], target=as.factor(y.train), 
                         type=1, cost=bestCost, bias=TRUE)
  pr <- FALSE
  linSVM.pred <- predict(bestModel, x.test[,selectedFeatures], proba=pr, decisionValues=TRUE)
  
  precisionVec <- vector(mode="numeric", length=28)
  recallVec <- vector(mode="numeric", length=28)
  classAccu <- vector(mode="numeric", length=28)
  
  for (cancerClass in 1:28) {
    truePos <- sum((linSVM.pred$predictions == cancerClass) & (y.test == cancerClass))
    falsePos <- sum((linSVM.pred$predictions == cancerClass) & (y.test != cancerClass))
    falseNeg <- sum((linSVM.pred$predictions != cancerClass) & (y.test == cancerClass))
    
    precisionVec[cancerClass] <- truePos / (truePos + falsePos)
    recallVec[cancerClass] <- truePos / (truePos + falseNeg)
    classAccu[cancerClass] <- truePos / sum(y.test == cancerClass)
  }
  
  res <- cbind(precisionVec, recallVec, classAccu)
  colnames(res) <- c("Precision", "Recall", "Class accuracy")
  
  return(res)  
}
########################################################################


ptm <- proc.time()
  allDataSetsRes <- mclapply(1:50, getPrecisionRecall, indexSets, testIdx, redMutMat, ranking.res, mc.cores=8)
proc.time() - ptm

## For Windows OS:
#sfInit(parallel=TRUE, cpus=8, type="SOCK")
#sfLibrary(LiblineaR)  # load library on workers
##sfExport("")  # export objects to workers
#ptm <- proc.time()
#  allDataSetsRes <- sfLapply(1:50, getPrecisionRecall, indexSets, testIdx, redMutMat, ranking.res)
#proc.time() - ptm
#sfStop()  # stop cluster (for 'snowfall')

save(allDataSetsRes, file="precisRecall_svm.RData")

# Plot
precisionMat <- matrix(data=0, nrow=28, ncol=50)
recallMat <- matrix(data=0, nrow=28, ncol=50)
classAccuMat <- matrix(data=0, nrow=28, ncol=50)
rownames(precisionMat) <- as.character(1:28)
rownames(recallMat) <- as.character(1:28)
rownames(classAccuMat) <- as.character(1:28)
for (v in 1:50) {
  precisionMat[,v] <- allDataSetsRes[[v]][,1]
  recallMat[,v] <- allDataSetsRes[[v]][,2]
  classAccuMat[,v] <- allDataSetsRes[[v]][,3]
}

classPrecis <- rowMeans(precisionMat, na.rm=TRUE)
precis.sErr <- 1.96 * (apply(precisionMat, 1, sd, na.rm=TRUE) / sqrt(50))
classRecall <- rowMeans(recallMat, na.rm=TRUE)
recall.sErr <- 1.96 * (apply(recallMat, 1, sd, na.rm=TRUE) / sqrt(50))
classAccu <- rowMeans(classAccuMat, na.rm=TRUE)

#par(mfrow=c(1,2))
precisBars <- barplot(classPrecis, ylim=c(0,1), 
                      main="Precision of SVM in cancer class prediction", 
                      xlab="Cancer classes", col="darkgrey", cex.names=0.7)
arrows(precisBars, classPrecis-precis.sErr, 
       precisBars, classPrecis+precis.sErr, 
       length=0.03, angle=90, code=3, lwd=2, col="black")
mtext(paste0("(Av. precision = ", round(mean(classPrecis), digits=3), ")"), side=3)

recallBars <- barplot(classRecall, ylim=c(0,1), 
                      main="Recall of SVM in cancer class prediction", 
                      xlab="Cancer classes", col="darkgrey", cex.names=0.7)
arrows(recallBars, classRecall-recall.sErr, 
       recallBars, classRecall+recall.sErr, 
       length=0.03, angle=90, code=3, lwd=2, col="black")
mtext(paste0("(Av. recall = ", round(mean(classRecall), digits=3), ")"), side=3)
#par(mfrow=c(1,1))
