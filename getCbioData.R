rm(list=ls())
setwd("~/Documents/Code/ctDNAdiagnostic")

## Retrieve mutation data
#########################
#install.packages("cgdsr")
library(cgdsr)
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
studies <- getCancerStudies(mycgds)
cancerTypes <- 28
myCaseList <- vector(mode="character", length=cancerTypes)
myGeneticProfile <- vector(mode="character", length=cancerTypes)

## Select the relevant studies
## Remarks: the indices for 'getCaseLists' and 'getGeneticProfiles' below may change 
## as cBioPortal make more studies available. It is important to check before download.
# Class 1: Bladder Urothelial Carcinoma (TCGA, Nature 2014)
blca <- studies[which(studies[,2]=="Bladder Urothelial Carcinoma (TCGA, Nature 2014)"), 1]
myCaseList[1] <- getCaseLists(mycgds,blca)[8,1]
myGeneticProfile[1] <- getGeneticProfiles(mycgds,blca)[7,1]  # 7 for pt mutations; 2 for CNAs
# Class 2: Breast Invasive Carcinoma (TCGA, Cell 2015)
brca <- studies[which(studies[,2]=="Breast Invasive Carcinoma (TCGA, Cell 2015)"), 1]
myCaseList[2] <- getCaseLists(mycgds,brca)[23,1]
myGeneticProfile[2] <- getGeneticProfiles(mycgds,brca)[10,1]  # 10 for pt mutations; 3 for CNAs
# Class 3: Colorectal Adenocarcinoma (TCGA, Nature 2012)
coadread <- studies[which(studies[,2]=="Colorectal Adenocarcinoma (TCGA, Nature 2012)"), 1]
myCaseList[3] <- getCaseLists(mycgds,coadread)[16,1]
myGeneticProfile[3] <- getGeneticProfiles(mycgds,coadread)[8,1]  # 8 for pt mutations; 1 for CNAs
# Class 4: Glioblastoma (TCGA, Cell 2013)
gbm <- studies[which(studies[,2]=="Glioblastoma (TCGA, Cell 2013)"), 1]
myCaseList[4] <- getCaseLists(mycgds,gbm)[7,1]
myGeneticProfile[4] <- getGeneticProfiles(mycgds,gbm)[7,1]  # 7 for pt mutations; 3 for CNAs
# Class 5: Head and Neck Squamous Cell Carcinoma (TCGA, Nature 2015)
hnsc <- studies[which(studies[,2]=="Head and Neck Squamous Cell Carcinoma (TCGA, Nature 2015)"), 1]
myCaseList[5] <- getCaseLists(mycgds,hnsc)[12,1]
myGeneticProfile[5] <- getGeneticProfiles(mycgds,hnsc)[8,1]  # 8 for pt mutations; 3 for CNAs
# Class 6: Kidney Renal Clear Cell Carcinoma (TCGA, Nature 2013)
kidrc <- studies[which(studies[,2]=="Kidney Renal Clear Cell Carcinoma (TCGA, Nature 2013)"), 1]
myCaseList[6] <- getCaseLists(mycgds,kidrc)[14,1]
myGeneticProfile[6] <- getGeneticProfiles(mycgds,kidrc)[10,1]  # 10 for pt mutations; 3 for CNAs
# Class 7: Acute Myeloid Leukemia (TCGA, NEJM 2013)
laml <- studies[which(studies[,2]=="Acute Myeloid Leukemia (TCGA, NEJM 2013)"), 1]
myCaseList[7] <- getCaseLists(mycgds,laml)[9,1]
myGeneticProfile[7] <- getGeneticProfiles(mycgds,laml)[7,1]  # 7 for pt mutations; 1 for CNAs
# Class 8: Lung Adenocarcinoma (TCGA, Nature 2014)
luad <- studies[which(studies[,2]=="Lung Adenocarcinoma (TCGA, Nature 2014)"), 1]
myCaseList[8] <- getCaseLists(mycgds,luad)[21,1]
myGeneticProfile[8] <- getGeneticProfiles(mycgds,luad)[10,1]  # 10 for pt mutations; 3 for CNAs
# Class 9: Lung Squamous Cell Carcinoma (TCGA, Nature 2012)
lusq <- studies[which(studies[,2]=="Lung Squamous Cell Carcinoma (TCGA, Nature 2012)"), 1]
myCaseList[9] <- getCaseLists(mycgds,lusq)[17,1]
myGeneticProfile[9] <- getGeneticProfiles(mycgds,lusq)[8,1]  # 8 for pt mutations; 1 for CNAs
# Class 10: Ovarian Serous Cystadenocarcinoma (TCGA, Nature 2011)
ov <- studies[which(studies[,2]=="Ovarian Serous Cystadenocarcinoma (TCGA, Nature 2011)"), 1]
myCaseList[10] <- getCaseLists(mycgds,ov)[18,1]
myGeneticProfile[10] <- getGeneticProfiles(mycgds,ov)[6,1]  # 6 for pt mutations; 1 for CNAs
# Class 11: Uterine Corpus Endometrioid Carcinoma (TCGA, Nature 2013)
ucec <- studies[which(studies[,2]=="Uterine Corpus Endometrial Carcinoma (TCGA, Nature 2013)"), 1]
myCaseList[11] <- getCaseLists(mycgds,ucec)[15,1]
myGeneticProfile[11] <- getGeneticProfiles(mycgds,ucec)[8,1]  # 8 for pt mutations; 2 for CNAs
# Class 12: Adenoid Cystic Carcinoma (MSKCC, Nat Genet 2013)
acyc <- studies[which(studies[,2]=="Adenoid Cystic Carcinoma (MSKCC, Nat Genet 2013)"), 1]
myCaseList[12] <- getCaseLists(mycgds,acyc)[4,1]
myGeneticProfile[12] <- getGeneticProfiles(mycgds,acyc)[2,1]  # 2 for pt mutations; 1 for CNAs
# Class 13: Brain Lower Grade Glioma (TCGA, Provisional)
blgg <- studies[which(studies[,2]=="Brain Lower Grade Glioma (TCGA, Provisional)"), 1]
myCaseList[13] <- getCaseLists(mycgds,blgg)[8,1]
myGeneticProfile[13] <- getGeneticProfiles(mycgds,blgg)[9,1]  # 9 for pt mutations; 3 for CNAs
# Class 14: Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (TCGA, Provisional)
csccea <- studies[which(studies[,2]=="Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (TCGA, Provisional)"), 1]
myCaseList[14] <- getCaseLists(mycgds,csccea)[8,1]
myGeneticProfile[14] <- getGeneticProfiles(mycgds,csccea)[8,1]  # 8 for pt mutations; 3 for CNAs
# Class 15: Kidney Renal Papillary Cell Carcinoma (TCGA, Provisional)
kidrp <- studies[which(studies[,2]=="Kidney Renal Papillary Cell Carcinoma (TCGA, Provisional)"), 1]
myCaseList[15] <- getCaseLists(mycgds,kidrp)[9,1]
myGeneticProfile[15] <- getGeneticProfiles(mycgds,kidrp)[10,1]  # 10 for mutations; 3 for CNAs
# Class 16: Liver Hepatocellular Carcinoma (AMC, Hepatology 2014)
livhc <- studies[which(studies[,2]=="Liver Hepatocellular Carcinoma (AMC, Hepatology 2014)"), 1]
myCaseList[16] <- getCaseLists(mycgds,livhc)[10,1]
myGeneticProfile[16] <- getGeneticProfiles(mycgds,livhc)[2,1]  # 2 for pt mutations; 1 for CNAs
# Class 17: Pancreatic Adenocarcinoma (TCGA, Provisional)
pancr <- studies[which(studies[,2]=="Pancreatic Adenocarcinoma (TCGA, Provisional)"), 1]
myCaseList[17] <- getCaseLists(mycgds,pancr)[8,1]
myGeneticProfile[17] <- getGeneticProfiles(mycgds,pancr)[8,1]  # 8 for pt mutations; 3 for CNAs
# Class 18: Prostate Adenocarcinoma (TCGA, Cell 2015)
prost <- studies[which(studies[,2]=="Prostate Adenocarcinoma (TCGA, Cell 2015)"), 1]
myCaseList[18] <- getCaseLists(mycgds,prost)[8,1]
myGeneticProfile[18] <- getGeneticProfiles(mycgds,prost)[7,1]  # 7 for pt mutations; 2 for CNAs
# Class 19: Skin Cutaneous Melanoma (TCGA, Provisional)
skcm <- studies[which(studies[,2]=="Skin Cutaneous Melanoma (TCGA, Provisional)"), 1]
myCaseList[19] <- getCaseLists(mycgds,skcm)[8,1]
myGeneticProfile[19] <- getGeneticProfiles(mycgds,skcm)[8,1]  # 8 for pt mutations; 3 for CNAs
# Class 20: Stomach Adenocarcinoma (TCGA, Nature 2014)
stom <- studies[which(studies[,2]=="Stomach Adenocarcinoma (TCGA, Nature 2014)"), 1]
myCaseList[20] <- getCaseLists(mycgds,stom)[11,1]
myGeneticProfile[20] <- getGeneticProfiles(mycgds,stom)[5,1]  # 5 for pt mutations; 1 for CNAs
# Class 21: Papillary Thyroid Carcinoma (TCGA, Cell 2014)
thyr <- studies[which(studies[,2]=="Papillary Thyroid Carcinoma (TCGA, Cell 2014)"), 1]
myCaseList[21] <- getCaseLists(mycgds,thyr)[8,1]
myGeneticProfile[21] <- getGeneticProfiles(mycgds,thyr)[8,1]  # 8 for pt mutations; 3 for CNAs
# Class 22: Adrenocortical Carcinoma (TCGA, Provisional)
adncc <- studies[which(studies[,2]=="Adrenocortical Carcinoma (TCGA, Provisional)"), 1]
myCaseList[22] <- getCaseLists(mycgds,adncc)[8,1]
myGeneticProfile[22] <- getGeneticProfiles(mycgds,adncc)[8,1]  # 8 for pt mutations; 3 for CNAs
# Class 23: Kidney Chromophobe (TCGA, Cancer Cell 2014)
kidchm <- studies[which(studies[,2]=="Kidney Chromophobe (TCGA, Cancer Cell 2014)"), 1]
myCaseList[23] <- getCaseLists(mycgds,kidchm)[7,1]
myGeneticProfile[23] <- getGeneticProfiles(mycgds,kidchm)[6,1]  # 6 for pt mutations; 1 for CNAs
# Class 24: Pheochromocytoma & Paraganglioma (TCGA, Provisional)
pcpg <- studies[which(studies[,2]=="Pheochromocytoma and Paraganglioma (TCGA, Provisional)"), 1]
myCaseList[24] <- getCaseLists(mycgds,pcpg)[8,1]
myGeneticProfile[24] <- getGeneticProfiles(mycgds,pcpg)[8,1]  # 8 for pt mutations; 3 for CNAs
# Class 25: Sarcoma (TCGA, Provisional)
sarc <- studies[which(studies[,2]=="Sarcoma (TCGA, Provisional)"), 1]
myCaseList[25] <- getCaseLists(mycgds,sarc)[8,1]
myGeneticProfile[25] <- getGeneticProfiles(mycgds,sarc)[8,1]  # 8 for pt mutations; 3 for CNAs
# Class 26: Testicular Germ Cell Cancer (TCGA, Provisional)
ttgc <- studies[which(studies[,2]=="Testicular Germ Cell Cancer (TCGA, Provisional)"), 1]
myCaseList[26] <- getCaseLists(mycgds,ttgc)[8,1]
myGeneticProfile[26] <- getGeneticProfiles(mycgds,ttgc)[8,1]  # 8 for pt mutations; 3 for CNAs
# Class 27: Uterine Carcinosarcoma (TCGA, Provisional)
ucs <- studies[which(studies[,2]=="Uterine Carcinosarcoma (TCGA, Provisional)"), 1]
myCaseList[27] <- getCaseLists(mycgds,ucs)[8,1]
myGeneticProfile[27] <- getGeneticProfiles(mycgds,ucs)[8,1]  # 8 for pt mutations; 3 for CNAs
# Class 28: Uveal Melanoma (TCGA, Provisional)
uvm <- studies[which(studies[,2]=="Uveal Melanoma (TCGA, Provisional)"), 1]
myCaseList[28] <- getCaseLists(mycgds,uvm)[7,1]
myGeneticProfile[28] <- getGeneticProfiles(mycgds,uvm)[6,1]  # 6 for pt mutations; 1 for CNAs


# Check that all are in order
myCaseList
myGeneticProfile

# Divide the long gene list into smaller batches
load("mGeneList_final.RData")   # pre-compiled list of somatic point mutated genes
load("cGeneList_final.RData")   # pre-compiled list of copy number altered genes
batchSize <- 90
batchNum <- (length(geneList) %/% batchSize) + 1
genesBatch <- vector(mode="list", length=batchNum)
for (i in 1:batchNum) {
  if (i==1) genesBatch[[1]] <- geneList[1:batchSize]
  else if (i==batchNum) genesBatch[[batchNum]] <- geneList[(batchNum*batchSize-batchSize+1):(length(geneList))]
  else genesBatch[[i]] <- geneList[(i*batchSize-batchSize+1):(i*batchSize)]
}

# Construct data matrix
#install.packages("RCurl")
library(RCurl)
curlSetOpt(timeout=200)
for (j in 1:cancerTypes) {
  # inner loop loops through all gene batches to construct data matrix for a cancer type
  for (k in 1:batchNum) {
    batchMat <- getProfileData(mycgds,genesBatch[[k]],myGeneticProfile[j],myCaseList[j])
    
    # for somatic point mutation matrix, we convert it to binary form in the next line
    # for CNA matrix, the next line is not needed as the raw data will be in the form of +2, +1, 0, -1, -2
    binBatchMat <- apply(batchMat, 1:2, function(x) if (is.na(x)|x=="NaN") return(0) else return(1))
    
    if (k==1) {
      binMat <- binBatchMat    
    } else {
      binMat <- cbind(binMat,binBatchMat)
    }
  }
  
  # outer loop loops through all cancer types and compiles the data matrix into combMutMat
  print(j)  # track progress
  label <- j
  if (j==1) {
    checkGenes <- colnames(binMat)
    combMutMat <- cbind(binMat,label)
  } else {
    if (all(checkGenes==colnames(binMat))) {
      combMutMat <- rbind(combMutMat,cbind(binMat,label))
    } else {
      print("Error: Column names do not match. Please check.")
      break
    }    
  }  
}
dim(combMutMat)
mmat.raw <- as.matrix(combMutMat)
save(mmat.raw, file="mmat_raw.RData")
#cmat.raw <- as.matrix(combMutMat)
#save(cmat.raw, file="cmat_raw.RData")


## For CNA matrix, create two variables for each gene, 
## one for presence/absence of amplification, the other for presence/absence of deletion
########################################################################################
load("mmat_raw.RData")
load("cmat_raw.RData")
yLabels <- cmat.raw[, ncol(cmat.raw)]
cmat.red <- cmat.raw[, -ncol(cmat.raw)]        # remove label column before conversion to binary format
cmat.red[is.na(cmat.red)] <- 0                 # convert NAs (including NaN) to '0'

numOfGenes <- ncol(cmat.raw) - 1
cmat.exp <- cmat.red[, rep(1:numOfGenes, each=2)]
cmat.exp[, seq(from=1,to=(numOfGenes*2-1),by=2)] <- (cmat.red > 0) + 0    # amplification
cmat.exp[, seq(from=2,to=(numOfGenes*2),by=2)] <- (cmat.red < 0) + 0      # deletion

# Rename column names by adding a suffix to the gene symbols
colnames(cmat.exp)[seq(from=1,to=(numOfGenes*2-1),by=2)] <- paste0(colnames(cmat.red), ".p")    # amplification
colnames(cmat.exp)[seq(from=2,to=(numOfGenes*2),by=2)] <- paste0(colnames(cmat.red), ".n")      # deletion

# Concatenate copy number matrix to somatic point mutation matrix
combMutMat <- cbind(mmat.raw[,-ncol(mmat.raw)], cmat.exp, yLabels)
dim(combMutMat)
save(combMutMat, file="mcmat_raw.RData")
