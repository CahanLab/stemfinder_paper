# CCAT

#Setup
library(SCENT)
library(gprofiler2)
setwd("~/Dropbox (CahanLab)/Kathleen.Noller/Stochasticity/")

#Load PPI and gene references
data(net17Jan16)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
library(org.Mm.eg.db)
ms = org.Mm.eg.db
library(org.Dr.eg.db) #zebrafish
dr = org.Dr.eg.db

#Load data - takes normalized and log transformed matrix
expDat = as.matrix(adata@assays$RNA@data)

  #Map genes to human homologs: gprofiler2
my.symbols <- rownames(expDat)
my.symbols.human = gorth(query = my.symbols, source_organism = "mmusculus", target_organism = "hsapiens", mthreshold=1, filter_na=T)
#my.symbols.human = gorth(query = my.symbols, source_organism = "celegans", target_organism = "hsapiens", mthreshold=1, filter_na=T)
genes_conversion = data.frame("mouse_geneid" = my.symbols, "human_geneid" = rep(NA, length(my.symbols)), "human_entrez" = rep(NA, length(my.symbols)))

  #convert to ENTREZ gene IDs
entrez = select(hs, keys = my.symbols.human$ortholog_name, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
head(entrez)
#entrez = AnnotationDbi::select(org.Hs.eg.db, keys=as.character(my.symbols.human$ortholog_name), columns = c('ENTREZID','SYMBOL'), keytype = 'SYMBOL')

for (i in 1:nrow(genes_conversion)){
  mouse_geneid = genes_conversion$mouse_geneid[i]
  if (length(my.symbols.human[my.symbols.human$input ==mouse_geneid,]$ortholog_name >0)){
    human_geneid = unique(my.symbols.human[my.symbols.human$input ==mouse_geneid,]$ortholog_name)[1]
    genes_conversion$human_entrez[i] = unique(entrez[entrez$SYMBOL == human_geneid,]$ENTREZID)
  } else{
    human_geneid = NA 
  }
  genes_conversion$human_geneid[i] = human_geneid
}

rownames(expDat) = genes_conversion$human_entrez
dim(expDat)
expDat_sub = expDat[!is.na(rownames(expDat)),]
dim(expDat_sub)

#Run CCAT
ccat.v <- CompCCAT(exp = expDat_sub, ppiA = net17Jan16.m)

#Add to adata
adata$ccat = ccat.v
adata$ccat_normandinvert = 1 - (adata$ccat)/max(adata$ccat) #reformat CytoTrace score


#########################
#CCAT function from github: https://github.com/aet21/SCENT/blob/master/R/CompCCAT.R

CompCCAT <- function(exp.m, ppiA.m){
  
  if(max(exp.m) > 100){ ### check if data has been log2-transformed or not and if not, then log2-transform with a pseudcount of 1, as for CCAT we can allow 0s.
    exp.m <- log2(exp.m+1);
  }
  # get input class of data matrix
  classMATRIX <- class(exp.m);
  
  # set common gene IDs
  commonEID.v <- intersect(rownames(ppiA.m),rownames(exp.m));
  
  # check the consistency of gene IDs and size of overlap
  #if ( length(commonEID.v) < 5000 ){
  #  stop(paste("The overlap of common genes between PPI and expression matrix is only ",length(commonEID.v)," so check that gene identifiers are correct. We don't recommend running CCAT on less than 5000 overlapping genes.",sep=""));
  #}
  
  # compute degrees
  k.v <- rowSums(ppiA.m[match(commonEID.v, rownames(ppiA.m)),]);
  
  if(classMATRIX=="matrix"){ ## ordinary matrix
    ccat.v <- as.vector(cor(exp.m[match(commonEID.v,rownames(exp.m)),],k.v));
  }
  else if (classMATRIX=="dgCMatrix"){
    ccat.v <- as.vector(corSparse(exp.m[match(commonEID.v,rownames(exp.m)),],Matrix(matrix(k.v,ncol=1))));
  }
  
  # subsample data set
  #  if (parallelMode) {
  #      chunk <- round(length(col_names)/subsamplesize)
  #      subsamples <- split(1:length(col_names), sample(factor(1:length(col_names) %% chunk)))
  #    }else if (length(col_names) > 10000) {
  #      chunk <- round(length(col_names)/1000)
  #      subsamples <- split(1:length(col_names), sample(factor(1:length(col_names) %% chunk)))
  #    }else{
  #    subsamples <- NULL
  #    }
  
  return(ccat.v);   
  
} ## EOF
