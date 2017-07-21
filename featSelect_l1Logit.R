rm(list=ls())
setwd("~/Documents/Code/ctDNAdiagnostic")

library(Matrix)
library(glmnet)
library(parallel)
#library(snow)       # package for parallel execution on Windows OS
#library(snowfall)   # package for parallel execution on Windows OS

load("mcmat_expBnrCNA2017.RData")       # gene alteration data matrix
load("testIdx_50bootstrapSets.RData")   # a list of 50 sets of indices for test cases
load("testIdx_global.RData")            # a set of indices as the final unseen test set

#######################################################################
## Function returns the index of features selected at each lambda value
getFeatures <- function(weights.list, dms) {
  featureList <- vector(mode="list",length=dms)
  
  for (k in 1:dms) {
    b1 <- weights.list$"1"[,k]
    b2 <- weights.list$"2"[,k]
    b3 <- weights.list$"3"[,k]
    b4 <- weights.list$"4"[,k]
    b5 <- weights.list$"5"[,k]
    b6 <- weights.list$"6"[,k]
    b7 <- weights.list$"7"[,k]
    b8 <- weights.list$"8"[,k]
    b9 <- weights.list$"9"[,k]
    b10 <- weights.list$"10"[,k]
    b11 <- weights.list$"11"[,k]
    b12 <- weights.list$"12"[,k]
    b13 <- weights.list$"13"[,k]
    b14 <- weights.list$"14"[,k]
    b15 <- weights.list$"15"[,k]
    b16 <- weights.list$"16"[,k]
    b17 <- weights.list$"17"[,k]
    b18 <- weights.list$"18"[,k]
    b19 <- weights.list$"19"[,k]
    b20 <- weights.list$"20"[,k]
    b21 <- weights.list$"21"[,k]
    b22 <- weights.list$"22"[,k]
    b23 <- weights.list$"23"[,k]
    b24 <- weights.list$"24"[,k]
    b25 <- weights.list$"25"[,k]
    b26 <- weights.list$"26"[,k]
    b27 <- weights.list$"27"[,k]
    b28 <- weights.list$"28"[,k]
    betaMat <- cbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,
                     b16,b17,b18,b19,b20,b21,b22,b23,b24,b25,b26,b27,b28)
    b.combined <- rowSums(betaMat)
    featureList[[k]] <- which(b.combined!=0)
  }
  return(featureList)  
}
#############################################################################
## Function returns no. of genes and test accuracy given a sequence of lambda
getL1logitRes <- function(i, indexSets, testIdx.global, dataMat, lambdaSeq) {
  allTestCases <- c(indexSets[[i]], testIdx.global)
  acc.vec <- vector(mode="numeric", length=length(lambdaSeq))
  
  lastColIdx <- ncol(dataMat)
  labels <- dataMat[,lastColIdx]
  x.train <- dataMat[-allTestCases,-lastColIdx]
  x.test <- dataMat[indexSets[[i]],-lastColIdx]
  y.train <- labels[-allTestCases]
  y.test <- labels[indexSets[[i]]]
  
  logit.fit <- glmnet(x.train, y.train, family="multinomial", standardize=FALSE,
                      alpha=1, lambda=lambdaSeq, maxit=150000, 
                      type.multinomial="ungrouped")
  logit.pred <- predict(logit.fit, newx=x.test, type="class")
  
  oDim <- length(logit.fit$lambda)
  for (j in 1:oDim) {
    acc.vec[j] <- mean(logit.pred[,j]==y.test)
  }
  # sometimes convergence for the smaller lambda values is not reached after 'maxit' iterations
  if (oDim < length(lambdaSeq)) {
    for (j in (oDim+1):length(lambdaSeq)) {
      acc.vec[j] <- NA
    }
  }
  
  geneList <- getFeatures(logit.fit$beta, oDim)
  
  all.res <- list(gene.num=logit.fit$df, test.acc=acc.vec, gene.sel=geneList)  
  return(all.res)  
}
#######################################################################


lambda.max <- 0.048
lambda.min <- 0.01 * 0.001
step.size <- (log(lambda.min)-log(lambda.max)) / 71
lambda.log <- seq(from=log(lambda.max), to=log(lambda.min), by=step.size)
length(lambda.log)   # to get finer intervals near the peak
lambda.log <- c(lambda.log[1:40], 
                seq(from=lambda.log[41], to=lambda.log[50], by=((lambda.log[50]-lambda.log[41])/18)), 
                lambda.log[51:72])
(lambda.seq <- exp(lambda.log))

ptm <- proc.time()
  l1Logit.res <- mclapply(1:50, getL1logitRes, indexSets, testIdx, redMutMat, lambda.seq, mc.cores=8)
proc.time() - ptm

# For Windows OS:
#sfInit(parallel=TRUE, cpus=8, type="SOCK")
#sfLibrary(glmnet)  # load library on workers
#sfExport("getFeatures")  # export objects to workers
#ptm <- proc.time()
#  l1Logit.res <- sfLapply(1:50, getL1logitRes, indexSets, testIdx, redMutMat, lambda.seq)
#proc.time() - ptm
#sfStop()  # stop cluster (for 'snowfall')

save(l1Logit.res, file="geneSelect_mcmat_l1logit.RData")

## Plot results
###############
names(l1Logit.res[[1]])
test.acc <- matrix(nrow=length(indexSets), ncol=length(lambda.seq))
gene.num <- matrix(nrow=length(indexSets), ncol=length(lambda.seq))

for (setNum in 1:length(indexSets)) {
  test.acc[setNum,] <- l1Logit.res[[setNum]]$test.acc
  if (length(l1Logit.res[[setNum]]$gene.num)==81) {
    gene.num[setNum,] <- l1Logit.res[[setNum]]$gene.num
  } else {
    gene.num[setNum,] <- c(l1Logit.res[[setNum]]$gene.num, 
                           rep(NA, times=(81-length(l1Logit.res[[setNum]]$gene.num))))
  }
}

av.acc <- colMeans(test.acc, na.rm=TRUE)
av.genes <- round(colMeans(gene.num, na.rm=TRUE))
sdev <- apply(test.acc, 2, sd, na.rm=TRUE)
sErr <- sdev / sqrt(50)
ci95 <- 1.96 * sErr
max(av.acc)
max.acc.index <- which(av.acc==max(av.acc))

#par(mfrow=c(1,2))
plot(log(lambda.seq), av.acc, main="Choice of lambda on accuracy of classifier", 
     xlab=expression("ln (" ~ lambda ~ ")"), ylab="Overall accuracy", type="b", pch=20, 
     ylim=range(c(av.acc-ci95, av.acc+ci95)))
arrows(log(lambda.seq), av.acc-ci95, log(lambda.seq), av.acc+ci95, length=0.03, angle=90, code=3, col="darkgrey")
abline(v=log(lambda.seq[max.acc.index]), col="red", lty=3, lwd=2)

plot(av.genes, av.acc, main="Performance of L1-logistic regression", 
     xlab="Number of genes selected (averaged)", ylab="Overall accuracy", 
     type="b", pch=20, ylim=range(c(av.acc-ci95, av.acc+ci95)))
arrows(av.genes, av.acc-ci95, av.genes, av.acc+ci95, length=0.03, angle=90, code=3, col="darkgrey")
abline(v=av.genes[max.acc.index], col="red", lty=3, lwd=2)
#par(mfrow=c(1,1))
