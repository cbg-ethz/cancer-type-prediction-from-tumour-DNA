rm(list=ls())
setwd("~/Documents/Code/ctDNAdiagnostic")
load("mmat_raw.RData")
load("mcmat_raw.RData")

## Remove uninformative rows and columns that are all zeros
###########################################################
# For point mutation matrix ...
numOfMutGenes.mmat <- rowSums(mmat.raw[,-ncol(mmat.raw)])
range(numOfMutGenes.mmat)
hist(numOfMutGenes.mmat, breaks=seq(0,(max(numOfMutGenes.mmat)+50),by=50), freq=TRUE, xlim=c(1,1000), 
     main="Distribution of patients with different number of point mutated genes", 
     xlab="Number of mutated genes", ylab="Number of patients")
sum(numOfMutGenes.mmat==0)  # no. of patients with zero mutated genes
combMutMat[(numOfMutGenes.mmat==0), lastColIdx]  # class labels of patients with zero mutated gene

numOfPatientWMut.mmat <- colSums(mmat.raw[,-ncol(mmat.raw)])
range(numOfPatientWMut.mmat)
sum(numOfPatientWMut.mmat==0)  # no. of genes that are not mutated in any patient

# For point mutation + CNA combined matrix ...
numOfMutGenes.mcmat <- rowSums(combMutMat[,-ncol(combMutMat)])
range(numOfMutGenes.mcmat)
hist(numOfMutGenes.mcmat, breaks=seq(0,(max(numOfMutGenes.mcmat)+50),by=50), freq=TRUE, xlim=c(1,4000), 
     main="Distribution of patients with different number of mutated genes", 
     xlab="Number of mutated genes", ylab="Number of patients")
sum(numOfMutGenes.mcmat==0)  # no. of patients with zero mutated genes
combMutMat[(numOfMutGenes.mcmat==0), lastColIdx]  # class labels of patients with zero mutated gene

numOfPatientWMut.mcmat <- colSums(combMutMat[,-ncol(combMutMat)])
range(numOfPatientWMut.mcmat)
sum(numOfPatientWMut.mcmat==0)  # no. of genes that are not mutated in any patient

# Check that no. of patients with zero mutated genes in M-C matrix are found in those in M matrix
which(numOfMutGenes.mcmat==0) %in% which(numOfMutGenes.mmat==0)

redMutMat <- mmat.raw[-which(numOfMutGenes.mmat==0), -which(numOfPatientWMut.mmat==0)]
dim(redMutMat)
save(redMutMat, file="mmat_final.RData")

## Check for duplicated columns
###############################
# Genes from the same chromosome regions that have their copy number altered will 
# have exactly the same 0s and 1s for all patients. We merge such genes into 
# a single feature, with a name like gene1-gene2.
redMutMat <- combMutMat[-which(numOfMutGenes.mmat==0), -which(numOfPatientWMut.mcmat==0)]
colInd.dupl <- duplicated(t(redMutMat))  # indicator of duplicated columns
sum(colInd.dupl)

colInd.allRep <- duplicated(t(redMutMat)) | duplicated(t(redMutMat), fromLast=TRUE)  # indicator of all repeated columns
sum(colInd.allRep)      # no. of columns that are repeated
colInd.uniqueRep <- colInd.allRep - colInd.dupl
sum(colInd.uniqueRep)   # no. of unique rows that are repeated
colIdx.uniqueRep <- which(colInd.uniqueRep==TRUE)

redMutMat2 <- redMutMat
for (z in 1:length(colIdx.uniqueRep)) {
  identCol <- which(apply(redMutMat2, 2, identical, redMutMat[,colIdx.uniqueRep[z]]))
  if (length(identCol) != 1) {
    #all(redMutMat2[, identCol[1]]==redMutMat2[, identCol[2]])
    colnames(redMutMat2)[identCol[1]] <- paste(names(which(apply(redMutMat2, 2, identical, 
                                                                 redMutMat[,colIdx.uniqueRep[z]]))), collapse="-")
    identCol <- identCol[-1]
    redMutMat2 <- redMutMat2[, -identCol]
    print(z)
  }
}
dim(redMutMat2)
sum(duplicated(t(redMutMat2)))    # verify that no. of duplicated columns is 0
redMutMat <- redMutMat2
save(redMutMat, file="mcmat_expBnrCNA2017.RData")

## Separate data into training and test sets
############################################
dim(redMutMat)
lastColIdx <- ncol(redMutMat)
labels <- redMutMat[,lastColIdx]
table(labels)

# Randomly draw out 1/4 of data from each class and set them aside as (global) test data
cancerTypes <- length(unique(labels))
endRange <- cumsum(table(labels))
temp <- endRange + 1
startRange <- c(1, temp[1:(length(temp)-1)])
drawSize <- round(table(labels)/4, digits=0)

testIdxList <- vector(mode="list", length=cancerTypes)
for (t in 1:cancerTypes) {
  testIdxList[[t]] <- sample(startRange[t]:endRange[t], drawSize[t], replace=FALSE)
}

testIdx <- unlist(testIdxList)
length(testIdx)
save(testIdx, file="testIdx_global.RData")

## Generate 50 bootstrapped sets of indices for training
########################################################
trainIdx <- 1:nrow(redMutMat)
trainIdx <- trainIdx[-testIdx]

which(match(endRange, trainIdx, nomatch=-1) == -1)
endRange[11] <- 3442
endRange[12] <- 3496
endRange[13] <- 3776
endRange[15] <- 4128
endRange[16] <- 4359
endRange[19] <- 5114
endRange[28] <- 6638
match(endRange, trainIdx)
endRange.new <- match(endRange, trainIdx)

temp <- endRange.new + 1
startRange.new <- c(1, temp[1:(length(temp)-1)])

drawSize.new <- vector(mode="numeric", length=cancerTypes)
drawSize.new[1] <- round(95/5)
for (u in 2:cancerTypes) 
  drawSize.new[u] <- round((endRange.new[u] - endRange.new[u-1]) / 5)

indexSets <- vector(mode="list", length=50)
for (v in 1:50){
  testIdxList <- vector(mode="list", length=cancerTypes)
  for (w in 1:cancerTypes) {
    testIdxList[[w]] <- sample(startRange.new[w]:endRange.new[w], drawSize.new[w], replace=FALSE)
  }
  temp2 <- unlist(testIdxList)
  indexSets[[v]] <- trainIdx[temp2]
}

# Double check
length(indexSets[[7]])
length(indexSets[[49]])
intersect(testIdx, indexSets[[23]])

save(indexSets, file="testIdx_50bootstrapSets.RData")
