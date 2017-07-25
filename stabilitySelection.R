rm(list=ls())
setwd("~/Documents/Code/ctDNAdiagnostic")
load("geneRank_mcmat_svmrfe.RData")     # gene ranking results from feature selection via SVM-RFE

numOfRuns <- 50
numOfGenes <- 50
topGenes <- matrix(nrow=numOfRuns, ncol=numOfGenes)

# Compute gene frequency
for (i in 1:numOfRuns) {
  topGenes[i,] <- ranking.res[[i]][1:numOfGenes]
}
gene.tab <- table(topGenes)
gene.freq <- sort(gene.tab, decreasing=TRUE) / numOfRuns
gene.num <- 1:length(gene.freq)

# Choosing the top 50 genes
indicator <- (gene.freq >= 0.3)
(numToSelect <- sum(indicator))
shortlistedGenes <- as.numeric(names(gene.freq[1:numToSelect]))
sort(gene.tab[match(as.character(shortlistedGenes), names(gene.tab))], decreasing=TRUE)
selectedGenes <- as.numeric(names(sort(gene.tab[match(as.character(shortlistedGenes), 
                                                      names(gene.tab))], decreasing=TRUE)))[1:50]

# For CNA genes, we pull out their corresponding amp/del genes since we 
# would also get those information when testing for copy number changes
cnaIdx <- which(selectedGenes > 7673)
load("mcmat_expBnrCNA2017.RData")       # load gene alteration data matrix
for (u in 1:length(cnaIdx)) {
  # read last character of each CNA gene
  lastCharPos <- nchar(colnames(redMutMat)[selectedGenes[cnaIdx[u]]])
  lastChar <- substr(colnames(redMutMat)[selectedGenes[cnaIdx[u]]], lastCharPos, lastCharPos)
  
  # if 'p', include next gene; if 'n', include previous gene
  if (lastChar == "p") {
    gene2add <- selectedGenes[cnaIdx[u]] + 1
  } else if (lastChar == "n") {
    gene2add <- selectedGenes[cnaIdx[u]] - 1
  }
  selectedGenes <- c(selectedGenes, gene2add)
}
selectedGenes <- sort(unique(selectedGenes))
matrix(colnames(redMutMat)[selectedGenes[19:82]], nrow=2, byrow=FALSE)  # double check
sum(selectedGenes <= 7673)      # number of genes associated with somatic point mutations
sum(selectedGenes > 7673) / 2   # number of genes associated with CNAs


## Check performance of most frequently selected genes
######################################################
library(LiblineaR)
dim(redMutMat)
lastColIdx <- ncol(redMutMat)
labels <- redMutMat[,lastColIdx]

load("testIdx_global.RData")            # load the set of unseen indices for the final test
x.train <- redMutMat[-testIdx,-lastColIdx]
x.test <- redMutMat[testIdx,-lastColIdx]
y.train <- labels[-testIdx]
y.test <- labels[testIdx]

tryCosts <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
bestCost <- NA
bestAcc <- 0
for (co in tryCosts) {
  acc <- LiblineaR(data=x.train[,selectedGenes], target=y.train, 
                   type=1, cost=co, bias=TRUE, cross=10, verbose=FALSE)
  cat("Results for C=",co," : ",acc," accuracy.\n",sep="")
  if (acc > bestAcc) {
    bestCost <- co
    bestAcc <- acc
  }
}
cat("Best cost is:", bestCost, "\n")
cat("Best accuracy is:", bestAcc, "\n")

bestLinMod <- LiblineaR(data=x.train[,selectedGenes], target=y.train, type=1, 
                        cost=bestCost, bias=TRUE, verbose=FALSE)
pr <- FALSE
linMod.pred <- predict(bestLinMod, x.test[,selectedGenes], proba=pr, decisionValues=TRUE)
mean(linMod.pred$predictions==y.test)


## Compute precision and recall
###############################
precisionVec <- vector(mode="numeric", length=28)
recallVec <- vector(mode="numeric", length=28)
for (cancerClass in 1:28) {
  truePos <- sum((linMod.pred$predictions == cancerClass) & (y.test == cancerClass))
  falsePos <- sum((linMod.pred$predictions == cancerClass) & (y.test != cancerClass))
  falseNeg <- sum((linMod.pred$predictions != cancerClass) & (y.test == cancerClass))
  
  precisionVec[cancerClass] <- truePos / (truePos + falsePos)
  recallVec[cancerClass] <- truePos / (truePos + falseNeg)
}
precisionVec[which(is.na(precisionVec))] <- 0
recallVec[which(is.na(recallVec))] <- 0
barplot(precisionVec, ylim=c(0,1), main="Precision of SVM in cancer class prediction", 
        xlab="Cancer classes", col="darkgrey", cex.names=0.7)
mtext(paste0("(Av. precision = ", round(mean(precisionVec), digits=3), ")"), side=3)
barplot(recallVec, ylim=c(0,1), main="Recall of SVM in cancer class prediction", 
        xlab="Cancer classes", col="darkgrey", cex.names=0.7)
mtext(paste0("(Av. recall = ", round(mean(recallVec), digits=3), ")"), side=3)
