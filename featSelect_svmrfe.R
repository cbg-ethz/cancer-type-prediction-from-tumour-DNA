rm(list=ls())
setwd("~/Documents/Code/ctDNAdiagnostic")

library(LiblineaR)
library(parallel)
#library(snow)       # package for parallel execution on Windows OS
#library(snowfall)   # package for parallel execution on Windows OS

load("mcmat_expBnrCNA2017.RData")       # gene alteration data matrix
load("testIdx_50bootstrapSets.RData")   # a list of 50 sets of indices for test cases
load("testIdx_global.RData")            # a set of indices as the final unseen test set

################### Function for parallel mSVM-RFE ###################
svmrfeRankFeatures <- function(i, indexSets, testIdx.global, dataMat) {
  testCases <- c(indexSets[[i]], testIdx.global)
  lastColIdx <- ncol(dataMat)
  x <- dataMat[-testCases, -lastColIdx]
  y <- dataMat[-testCases, lastColIdx]
  
  p <- ncol(x)
  
  survivingFeaturesIndexes <- seq(1:p)
  featureRankedList <- vector(mode="integer", length=p)
  rankedFeatureIndex <- p
  
  while (length(survivingFeaturesIndexes) > 1) {
    
    # Remarks: the optimum parameter for svm seems to increase as the 
    # number of genes decreases. Hence, we do a simple parameter search 
    # when the number of genes gets smaller.
    if (rankedFeatureIndex > 500) {
      bestCost <- 0.02
      
    } else {
      # find the best parameter for svm
      tryCosts <- c(0.02,0.07,0.2,0.7)
      bestCost <- NA
      bestAcc <- 0
      for (ct in tryCosts) {
        acc <- LiblineaR(data=x[,survivingFeaturesIndexes], target=as.factor(y),
                         type=1, cost=ct, bias=TRUE, cross=10)
        if (acc > bestAcc) {
          bestCost <- ct
          bestAcc <- acc
        }
      }
    }
    
    bestSVM <- LiblineaR(data=x[,survivingFeaturesIndexes], target=as.factor(y), 
                         type=1, cost=bestCost, bias=TRUE)
    
    # compute the weight vector
    biasIdx <- ncol(bestSVM$W)
    W.mat <- bestSVM$W[,-biasIdx]
    
    # compute ranking criteria
    multiclassWeights <- W.mat * W.mat
    rankingCriteria <- colMeans(multiclassWeights)
    
    # rank the features
    ranking <- sort(rankingCriteria, index.return=TRUE)$ix  # ascending order
    
    # update feature ranked list
    #featureRankedList[rankedFeatureIndex] <- survivingFeaturesIndexes[ranking[1]]
    #rankedFeatureIndex <- rankedFeatureIndex - 1
    numToRmv <- ceiling(0.005 * length(survivingFeaturesIndexes))
    for (i in 1:numToRmv) {
      featureRankedList[rankedFeatureIndex] <- survivingFeaturesIndexes[ranking[i]]
      rankedFeatureIndex <- rankedFeatureIndex - 1
    }
    
    # eliminate the feature with smallest ranking criterion
    #survivingFeaturesIndexes <- survivingFeaturesIndexes[-ranking[1]]
    survivingFeaturesIndexes <- survivingFeaturesIndexes[-ranking[1:numToRmv]]    
  }
  
  featureRankedList[1] <- survivingFeaturesIndexes
  
  return(featureRankedList)
}
################# Function for computing accuracies ##################
getError <- function(j, rankingList, numToUse, dataMat, indexSets, testIdx.global){
  
  allTestCases <- c(indexSets[[j]], testIdx.global)
  lastColIdx <- ncol(dataMat)
  labels <- dataMat[,lastColIdx]
  x.train <- dataMat[-allTestCases,-lastColIdx]
  x.test <- dataMat[indexSets[[j]],-lastColIdx]
  y.train <- labels[-allTestCases]
  y.test <- labels[indexSets[[j]]]
  
  featureRankedList <- rankingList[[j]]
  progressive.acc <- vector(mode="numeric", length=length(numToUse))
    
  for (k in 1:length(numToUse)) {
    selectedFeatures <- featureRankedList[1:numToUse[k]]
    
    # pull out CNA genes
    #if (sum(selectedFeatures > 7693)) {
    #  cnaIdx <- which(selectedFeatures > 7693)
    #  
    #  for (u in 1:length(cnaIdx)) {
    #    # read last character of each CNA gene
    #    lastCharPos <- nchar(colnames(dataMat)[selectedFeatures[cnaIdx[u]]])
    #    lastChar <- substr(colnames(dataMat)[selectedFeatures[cnaIdx[u]]], lastCharPos, lastCharPos)
    #    
    #    # if 'p', include next gene; if 'n', include previous gene
    #    if (lastChar == "p") {
    #      gene2add <- selectedFeatures[cnaIdx[u]] + 1
    #    } else if (lastChar == "n") {
    #      gene2add <- selectedFeatures[cnaIdx[u]] - 1
    #    }
    #    selectedFeatures <- c(selectedFeatures, gene2add)
    #  }
    #}
    #selectedFeatures <- sort(unique(selectedFeatures))
    
    if (numToUse[k] > 500) {
      bestCost <- 0.02
    } else {      
      # find the best parameter for svm
      tryCosts <- c(0.02,0.07,0.2,0.7)
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
    }
    
    bestModel <- LiblineaR(data=x.train[,selectedFeatures], target=as.factor(y.train),
                           type=1, cost=bestCost, bias=TRUE)
    pr <- FALSE
    linSVM.pred <- predict(bestModel, x.test[,selectedFeatures], proba=pr, 
                           decisionValues=TRUE)
    progressive.acc[k] <- mean(linSVM.pred$predictions==y.test)
  }
  
  return(progressive.acc)
}
########################## End of functions ##########################


## Parallelised feature ranking with mSVM-RFE
#############################################
ptm <- proc.time()
  ranking.res <- mclapply(1:50, svmrfeRankFeatures, indexSets, testIdx, redMutMat, mc.cores=8)
proc.time() - ptm

# For Windows OS:
#sfInit(parallel=TRUE, cpus=8, type="SOCK")
#sfLibrary(LiblineaR)  # load library on workers
##sfExport("")  # export objects to workers
#ptm <- proc.time()
#  ranking.res <- sfLapply(1:50, svmrfeRankFeatures, indexSets, testIdx, redMutMat)
#proc.time() - ptm
#sfStop()  # stop cluster (for 'snowfall')

save(ranking.res, file="geneRank_mcmat_svmrfe.RData")

## Compute accuracy of ranked genes on test data set
####################################################
numOfFeatures <- c(10, 20, 50, seq(from=70,to=160,by=30), seq(from=200,to=7600,by=50), 7673)
#numOfFeatures <- c(10, 20, 50, seq(from=70,to=160,by=30), seq(from=200,to=8200,by=50))

ptm <- proc.time()
  output <- mclapply(1:50, getError, ranking.res, numOfFeatures, redMutMat, indexSets, testIdx, mc.cores=8)
proc.time() - ptm

# For Windows OS:
#sfInit(parallel=TRUE, cpus=8, type="SOCK")
#sfLibrary(LiblineaR)  # load library on workers
##sfExport("")  # export objects to workers
#ptm <- proc.time()
#  output <- sfLapply(1:50, getError, ranking.res, numOfFeatures, redMutMat, indexSets, testIdx)
#proc.time() - ptm
#sfStop()  # stop cluster (for 'snowfall')

test.acc <- matrix(unlist(output, use.names=FALSE), ncol=length(numOfFeatures), byrow=TRUE)
save(test.acc, file="accuracyTrace_mcmat_svmrfe.RData")

## Plot results
###############
av.acc <- colMeans(test.acc)
sdev <- apply(test.acc, 2, sd)
sErr <- sdev / sqrt(50)
ci95 <- 1.96 * sErr
max(av.acc)
max.acc.index <- which(av.acc==max(av.acc))

plot(x=numOfFeatures, y=av.acc, main="Performance of SVM-RFE", xlab="Number of genes used", ylab="Overall accuracy", 
     type="b", pch=20, ylim=range(c(av.acc-ci95, av.acc+ci95)))
arrows(numOfFeatures, av.acc-ci95, numOfFeatures, av.acc+ci95, length=0.03, angle=90, code=3, col="darkgrey")
abline(v=numOfFeatures[max.acc.index], col="red", lty=3, lwd=2)
