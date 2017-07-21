rm(list=ls())

library(ranger)
library(foreach)
library(parallel)
library(doParallel)

load("mcmat_expBnrCNA2017.RData")       # gene alteration data matrix
load("testIdx_50bootstrapSets.RData")   # a list of 50 sets of indices for test cases
load("testIdx_global.RData")            # a set of indices as the final unseen test set

########## Funtion returns OOB accuracy and variable importance ranking ##########
getVarImpt <- function(i, testIdxSets, testIdxGlobal, dataMat) {
  allTestCases <- c(testIdxSets[[i]], testIdxGlobal)
  trainData <- dataMat[-allTestCases, ]
  lastColIdx <- ncol(dataMat)
  splitNum <- 9 * round( sqrt(lastColIdx-1) )
  # remarks: splitNum is determined from earlier optimisation of random forest on the data set
  # in our case, it should be '2' for somatic point mutation matrix, '9' for spm + cna matrix
  
  rdmWald.fit <- ranger(classLabels ~ ., data=as.data.frame(trainData), num.trees=500, mtry=splitNum, 
                        importance="permutation", num.threads=1, verbose=FALSE, classification=TRUE)
  
  trng.acc <- 1 - rdmWald.fit$prediction.error
  rankedFeatIdx <- order(rdmWald.fit$variable.importance, decreasing=TRUE)
  
  final.res <- list(oob.acc=trng.acc, ranked.feat=rankedFeatIdx)
  return(final.res)
}
##################################################################################
## Function grows a random forest on the training data using selected features, ##
## and returns the accuracy on test data #########################################
multiRdmWald <- function(i, dataMat, geneRanking, numToUse, testIdxGlobal, testIdxSets, nthSet) {
  allTestCases <- c(testIdxSets[[nthSet]], testIdxGlobal)
  lastColIdx <- ncol(dataMat)
  selectedColms <- c(geneRanking[1:numToUse[i]], lastColIdx)
  trainData <- dataMat[-allTestCases, selectedColms]
  testData <- dataMat[testIdxSets[[nthSet]], selectedColms]
  yLabels.test <- dataMat[testIdxSets[[nthSet]], lastColIdx]
  
  rdmWald.fit <- ranger(classLabels ~ ., data=as.data.frame(trainData), num.trees=500, num.threads=1, 
                        verbose=FALSE, classification=TRUE)
  
  rdmWald.pred <- predict(rdmWald.fit, as.data.frame(testData), type="response", num.threads=1)
  overallAccu <- mean(rdmWald.pred$predictions==yLabels.test)
  
  return(overallAccu)
}
##################################################################################


## Compute variable importance of the features
##############################################
clus1 <- makeForkCluster(8)
registerDoParallel(clus1)
ptm <- proc.time()
varImptRanking <- foreach(u=1:50, .inorder=TRUE, .packages="ranger"
                          ) %dopar% getVarImpt(u, indexSets, testIdx, redMutMat)
stopCluster(clus1)

elapsedTime <- proc.time() - ptm
cat("Time taken: ", elapsedTime[3], " s\n", sep="")
save(varImptRanking, file="varImptRes_mcmat_rf.RData")

## Compute accuracy of ranked genes on test data set
####################################################
rankedFeatSets <- matrix(data=0, nrow=50, ncol=length(varImptRanking[[1]]$ranked.feat))
for (i in 1:50) {
  rankedFeatSets[i,] <- varImptRanking[[i]]$ranked.feat
}

numOfGenes <- c(10,20,50,70, seq(from=100,to=1800,by=50), seq(from=2000,to=3000,by=100), seq(from=3200,to=7400,by=200), 7673)
#numOfGenes <- c(10,20,50,70, seq(from=100,to=1800,by=50), seq(from=2000,to=3000,by=100), seq(from=3200,to=8400,by=200))

test.acc <- matrix(nrow=50, ncol=length(numOfGenes))
forestClus <- makeForkCluster(8)
registerDoParallel(forestClus)
ptm <- proc.time()
for(u in 1:50) {
  rankedFeatIdx <- rankedFeatSets[u,]
  test.acc[u,] <- foreach(v=1:length(numOfGenes), .combine="c", .inorder=TRUE,
                          .packages="ranger") %dopar% multiRdmWald(
                            v, redMutMat, rankedFeatIdx, numOfGenes, testIdx, indexSets, u)
}
stopCluster(forestClus)

elapsedTime <- proc.time() - ptm
cat("Time taken: ", elapsedTime[3], " s\n", sep="")
save(test.acc, file="accuracyTrace_mcmat_rf.RData")
