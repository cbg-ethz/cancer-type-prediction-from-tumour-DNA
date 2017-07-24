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

################# Function for computing accuracies ##################
compareAccu <- function(u, rankingList, numToUse, dataMat, indexSets, testIdx.global){
  
  allTestCases <- c(indexSets[[u]], testIdx.global)
  lastColIdx <- ncol(dataMat)
  labels <- dataMat[,lastColIdx]
  x.train <- dataMat[-allTestCases,-lastColIdx]
  x.test <- dataMat[indexSets[[u]],-lastColIdx]
  y.train <- labels[-allTestCases]
  y.test <- labels[indexSets[[u]]]
  
  featureRankedList <- rankingList[[u]]
  progressive.acc <- vector(mode="numeric", length=3*length(numToUse))
  
  for (v in 1:length(numToUse)) {
    featureIdxVec <- 1:numToUse[v]
    
    for (w in 1:3) {
      if (w != 1) featureIdxVec <- featureIdxVec + numToUse[v]
      selectedFeatures <- featureRankedList[featureIdxVec]
      
      if (numToUse[v] > 500) {
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
      progressive.acc[(3*v) - 3 + w] <- mean(linSVM.pred$predictions==y.test)
    }    
  }
  
  return(progressive.acc)
}
########################## End of functions ##########################

## Compute accuracies of top 1st-, 2nd-, 3rd- tier predictors
#############################################################
numOfFeatures <- c(10, 20, 50, seq(from=70,to=160,by=30), seq(from=200,to=2500,by=50))

ptm <- proc.time()
  output <- mclapply(1:50, compareAccu, ranking.res, numOfFeatures, redMutMat, indexSets, testIdx, mc.cores=8)
proc.time() - ptm

## For Windows OS:
#sfInit(parallel=TRUE, cpus=8, type="SOCK")
#sfLibrary(LiblineaR)  # load library on workers
##sfExport("")  # export objects to workers
#ptm <- proc.time()
#  output <- sfLapply(1:50, compareAccu, ranking.res, numOfFeatures, redMutMat, indexSets, testIdx)
#proc.time() - ptm
#sfStop()  # stop cluster (for 'snowfall')

test.acc <- matrix(unlist(output, use.names=FALSE), ncol=3*length(numOfFeatures), byrow=TRUE)
save(test.acc, file="topFeatPerform_mcmat_svmrfe.RData")

## Visualisation of performances of top feature sets
####################################################
library(ggplot2)
library(grid)
library(scales)
load("topFeatPerform_mcmat_svmrfe.RData")

# Re-structure the data into data frame, in long format
sErr <- apply(test.acc, 2, sd) / sqrt(50)
featSetPerform <- data.frame(numOfGenes = rep(numOfFeatures, each=3), 
                             avAccu = colMeans(test.acc), 
                             lwr95CI = colMeans(test.acc) - (1.96 * sErr),
                             upr95CI = colMeans(test.acc) + (1.96 * sErr),
                             geneSets = rep(c("Top gene set","2nd best set","3rd best set"), times=length(numOfFeatures)), 
                             stringsAsFactors=FALSE)

# Plot
colBrewPalette <- c("turquoise4", "#003366", "firebrick4")
ggplot(data=featSetPerform, aes(x=numOfGenes, y=avAccu, group=geneSets, colour=geneSets, shape=geneSets)) + 
  geom_ribbon(data=featSetPerform.mcmat, aes(ymin=lwr95CI, ymax=upr95CI, fill=geneSets, linetype=NA), alpha=0.8) + 
  geom_line(data=featSetPerform.mcmat, linetype="solid", size=0.3) + 
  scale_colour_manual(name="Predictor sets", values=colBrewPalette, breaks=c("Top gene set","2nd best set","3rd best set")) + 
  scale_fill_manual(name="Predictor sets", values=c("paleturquoise", "skyblue1", "salmon"), 
                    breaks=c("Top gene set","2nd best set","3rd best set")) + 
  xlab("Number of genes used") + ylab("Overall accuracy") + ggtitle("") + 
  scale_x_continuous(breaks=seq(from=0,to=2500,by=500)) + 
  scale_y_continuous(breaks=seq(from=0.1,to=0.9,by=0.2),limits=c(0.1,0.9),labels=percent,oob=squish) + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), aspect.ratio=1, 
                     axis.text=element_text(size=8), axis.title=element_blank(), plot.title=element_blank(), 
                     #axis.title.x=element_text(size=24,vjust=-1), axis.title.y=element_text(size=24,vjust=2),
                     legend.title=element_blank(), legend.text=element_text(size=8), 
                     legend.position=c(0.75,0.25), plot.margin=unit(c(0.3,0.3,0.7,0.7), "lines"))
