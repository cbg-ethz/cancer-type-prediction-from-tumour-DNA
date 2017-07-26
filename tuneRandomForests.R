rm(list=ls())
setwd("~/Documents/Code/ctDNAdiagnostic")

library(ranger)
load("mcmat_expBnrCNA2017.RData")       # gene alteration data matrix
load("testIdx_50bootstrapSets.RData")   # a list of 50 sets of indices for test cases
load("testIdx_global.RData")            # a set of indices as the final unseen test set

## Optimise mtry parameter
##########################
# Compute accuracy for a subset of the training sets, using different mtry values
overallAccy <- matrix(data=0.0, nrow=10, ncol=10)
colnames(overallAccy) <- paste0(1:10, "rootP")
choiceOfSet <- 16:25
for (v in 1:length(choiceOfSet)) {
  allTestCases <- c(indexSets[[choiceOfSet[v]]], testIdx)
  lastColIdx <- ncol(redMutMat)
  trainData <- redMutMat[-allTestCases, ]
  testData <- redMutMat[indexSets[[choiceOfSet[v]]], ]
  yLabels.test <- redMutMat[indexSets[[choiceOfSet[v]]], lastColIdx]
  
  for (w in 1:ncol(overallAccy)) {
    rdmWald.fit <- ranger(classLabels ~ ., data=as.data.frame(trainData), num.trees=500, mtry=(w*round(sqrt(lastColIdx-1))), 
                          num.threads=1, verbose=FALSE, classification=TRUE)
    rdmWald.pred <- predict(rdmWald.fit, as.data.frame(testData), type="response", num.threads=1)
    overallAccy[v,w] <- mean(rdmWald.pred$predictions==yLabels.test)
  }  
  print(v)
}

# Visualisation
accu <- colMeans(overallAccy)
ci95 <- 1.96 * (apply(overallAccy, 2, sd) / sqrt(length(choiceOfSet)))
plot(x=1:10, y=accu, type="b", pch=20, main="", xlab="mtry value", ylab="Overall accuracy", 
     ylim=range(c(accu-ci95, accu+ci95)))
arrows(1:10, accu-ci95, 1:10, accu+ci95, length=0.03, angle=90, code=3, col="darkgrey")


## Check if the default number of trees used (500) is sufficient for best mtry value
####################################################################################
bestMtry <- which(accu==max(accu))
overallAccy.1kTrees <- vector(mode="numeric", length=10)
for (v in 1:length(choiceOfSet)) {
  allTestCases <- c(indexSets[[choiceOfSet[v]]], testIdx)
  lastColIdx <- ncol(redMutMat)
  trainData <- redMutMat[-allTestCases, ]
  testData <- redMutMat[indexSets[[choiceOfSet[v]]], ]
  yLabels.test <- redMutMat[indexSets[[choiceOfSet[v]]], lastColIdx]
  
  rdmWald.fit <- ranger(classLabels ~ ., data=as.data.frame(trainData), 
                        num.trees=1000, mtry=(bestMtry*round(sqrt(lastColIdx-1))), 
                        num.threads=1, verbose=FALSE, classification=TRUE)
  rdmWald.pred <- predict(rdmWald.fit, as.data.frame(testData), type="response", num.threads=1)
  overallAccy.1kTrees[v] <- mean(rdmWald.pred$predictions==yLabels.test)
  print(v)
}
avAccu.1kTrees <- mean(overallAccy.1kTrees)
ci95.1kTrees <- 1.96 * (sd(overallAccy.1kTrees) / sqrt(length(choiceOfSet)))
dodge <- 0.3
points(x=(bestMtry+dodge), y=avAccu.1kTrees, pch=20, col="red")
arrows((bestMtry+dodge), avAccu.1kTrees-ci95.1kTrees, (bestMtry+dodge), avAccu.1kTrees+ci95.1kTrees, 
       length=0.03, angle=90, code=3, col="darkgrey")

## Closing remarks: 
# These parameters should be used only for gene ranking, but not when computing accuracy trace! 
# In the latter case, the default mtry value is good enough for small number of genes.
