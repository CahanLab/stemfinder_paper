### SingleCellNet in R
### From https://pcahan1.github.io/singleCellNet/

#Setup 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(rgl)
library(singleCellNet)

# Load training and query data as Seurat objects 

  #Training data 
TM = readRDS("../../LLi_data/LLibatch_epi_results_withvelocity_10082021.rds")
expTM = TM@assays$RNA@counts 
stTM = TM@meta.data 
stTM = droplevels(stTM)
stTM$barcodes = rownames(stTM)

  #Query data
adata = readRDS("~/Downloads/olbrecht_tumorepicellsonly_results_113021.rds")

#Load mouse-human ortholog table for species conversion if needed

#oTab = utils_loadObject("~/Documents/human_mouse_genes_Jul_24_2018.rda")
#expQuery = adata@assays$RNA@counts
#aa = csRenameOrth(expTM, expQuery, oTab)
#expQueryOrth = aa[['expTrain']] #this is correct if converting to mouse; reverse it if converting training to human
#expTrainOrth = aa[['expQuery']]

#If no species conversion needed, proceed: 

expQueryOrth = adata@assays$RNA@counts
expTrainOrth = expTM
commonGenes = intersect(rownames(expQueryOrth), rownames(expTrainOrth)) 
length(commonGenes)
expTM2 = expTrainOrth[,rownames(stTM)]

#Train classifier

#split training data
stList = splitCommon(stTM, ncells=200, dLevel = "gini_potency_bin")
stTrain = stList[[1]]
expTrain = expTrainOrth[commonGenes,rownames(stTrain)]

#train classifier
system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 10, nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = "gini_potency_bin", colName_samp = "barcodes"))

# Assess it with held out data

#validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=200, dLevel="gini_potency_bin") 
stTest = stTestList[[1]]
expTest = expTrainOrth[,rownames(stTest)]

#predict
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)

# PR curves
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "barcodes", classTrain = "gini_potency_bin", classQuery = "gini_potency_bin", nRand = 50)
plot_PRs(tm_heldoutassessment)

# Apply to query data
nqRand = 50
system.time(crQuery<-scn_predict(class_info[['cnProc']], expQueryOrth, nrand=nqRand))

# Heat map
stQuery = adata@meta.data
sgrp = as.vector(stQuery$Phenotype)
names(sgrp) = as.vector(rownames(stQuery))
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)

sc_hmClass(crQuery, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)

# Attribution plot
stQuery$barcode = rownames(stQuery)
plot_attr(crQuery, stQuery, nrand=nqRand, sid="barcode", dLevel="Phenotype")

# Classify cells 
  ## cell is classified as category with highest classification score or higher than a classification score threshold of your choosing.
  ## returns annotations in metadata column 
stQuery <- get_cate(classRes = crQuery, sampTab = stQuery, dLevel = "seurat_clusters", sid = "barcode", nrand = nqRand)
head(stQuery)


### function get_cate
#' @export
assign_cate <- function (classRes, sampTab, cThresh = 0) 
{
  topCats <- rownames(classRes)[apply(classRes, 2, which.max)]
  sampTab <- cbind(sampTab, category = topCats)
  sampTab
}

#' @export
get_cate <- function (classRes, sampTab, dLevel, sid, nrand, cThresh=0, keepRand = FALSE) 
{
  if(is.data.frame(classRes)){
    classRes = as.matrix(classRes)
  }
  if(nrand == 0){
    stTmp = sampTab[,c(sid, dLevel)]
  }else{
    stTmp <- addToST(sampTab, nrand = nrand, sid = sid, dLevels = dLevel)
  }
  #stTmp <- assign_cate(classRes, stTmp)
  colnames(stTmp)[2] <- "group"
  
  topCat_score=c()
  topCats = c()
  
  for(i in 1:ncol(classRes)){
    tmp =  max(classRes[,i])
    topCat_score=c(topCat_score, tmp)
    if(tmp < cThresh){
      tmp2 = "rand"
      topCats = c(topCats, tmp2)
    }else{
      tmp2 = names(classRes[,i][classRes[,i] == tmp])[1]
      topCats = c(topCats, tmp2)
    }
    
  }
  
  if(keepRand){
    sampTab <- cbind(stTmp, category = topCats, scn_score = topCat_score)
    return(sampTab)
  }
  
  sampTab <- cbind(sampTab, category = topCats[1:nrow(sampTab)], scn_score = topCat_score[1:nrow(sampTab)])
  return(sampTab) 
}
