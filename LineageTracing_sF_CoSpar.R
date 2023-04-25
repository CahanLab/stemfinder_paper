### CoSpar lineage tracing dataset
## Mar 17, 2023 KNoller

#Setup
library(Seurat)
library(stemFinder)
setwd("~/Dropbox (CahanLab)/Kathleen.Noller/Stochasticity/data_from_pc/")
load("clone.rda") #matrix of cell barcodes x clone IDs
load("times.rda") #time points for each cell barcode: 2, 4, or 6 days
adata = readRDS("adata_cospar_seurat_030823.rds")

#Limit analysis to only clones that contain cells from time point 6day 
clones_incl = c()
for (c in colnames(clone_bin)){
  times_c = names(table(adata@meta.data[names(which(clone_bin[,c] == TRUE)),]$time_info))
  if("6" %in% times_c){
    clones_incl = c(clones_incl, c)
  }
}

#Count how many different lineages a given clone gives rise to
lin_counts = data.frame("clone_ID" = colnames(clone_bin), "Num_lineages" = rep("", 500), "Mean_sFinvert_undiff" = rep("", 500), "Min_sFinvert_undiff" = rep("", 500), "Mean_CCATinvert_undiff" = rep("",500), "Min_CCATinvert_undiff" = rep("", 500), "Mean_Cytoinvert_undiff" = rep("",500), "Min_Cytoinvert_undiff" = rep("", 500))
for (c in colnames(clone_bin)){
  #count number of lineages
  lins = unique(names(table(adata@meta.data[names(which(clone_bin[,c] == TRUE)),]$state_info)))
  if("undiff" %in% times_c){
    lins = length(lins) - 1
  }else {lins = length(lins)}
  lin_counts[lin_counts$clone_ID == c,]$Num_lineages = lins
  
  #retrieve stemFinder inverted score for time point 2 undifferentiated cells in that lineage
    ## 2 ways: minimum sF score per undiff cells vs. avg sF score for undiff cells
  bars = names(which(clone_bin[,c] == TRUE))
  min_time = min(adata@meta.data[bars,]$time_info)
  avg_sF = mean(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$stemFinder_invert)
  min_sF = min(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$stemFinder_invert)
  lin_counts[lin_counts$clone_ID == c,]$Min_sFinvert_undiff = min_sF
  lin_counts[lin_counts$clone_ID == c,]$Mean_sFinvert_undiff = avg_sF
  
  avg_CCAT = mean(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$ccat_normandinvert)
  min_CCAT = min(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$ccat_normandinvert)
  lin_counts[lin_counts$clone_ID == c,]$Min_CCATinvert_undiff = min_CCAT
  lin_counts[lin_counts$clone_ID == c,]$Mean_CCATinvert_undiff = avg_CCAT

  avg_Cyto = mean(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$CytoTRACE_normandinvert)
  min_Cyto = min(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$CytoTRACE_normandinvert)
  lin_counts[lin_counts$clone_ID == c,]$Min_Cytoinvert_undiff = min_Cyto
  lin_counts[lin_counts$clone_ID == c,]$Mean_Cytoinvert_undiff = avg_Cyto
}
  

#Plots

  ## Num lineages vs. sF inverted score
lin_counts$Num_lineages = as.numeric(lin_counts$Num_lineages)
#ggplot(lin_counts, aes(x = Min_sFinvert_undiff, y = Num_lineages)) + geom_point()
cor.test(x = as.numeric(lin_counts$Min_sFinvert_undiff), y = lin_counts$Num_lineages, method = 'pearson')

ggplot(lin_counts[,c(2,4)], aes(y = as.numeric(Min_sFinvert_undiff), x = as.numeric(Num_lineages))) + geom_boxplot(aes(fill = Num_lineages)) + guides(fill = 'none', color = 'none') + geom_point(aes(color = as.numeric(Min_sFinvert_undiff))) + xlab("Number of downstream lineages") + ylab("Minimum inverted stemFinder score") + ggtitle("Inverted stemFinder vs. number of lineages")

  #For given lineage, dot plot of sF score
lineage = 'Monocyte'
marker = 'Mmp8'

      #find barcodes of undifferentiated cells that give rise to lineage of interest
clones_interest = c()
cells_interest = c()
for (c in colnames(clone_bin)){
  ds_names = (names(table(adata@meta.data[names(which(clone_bin[,c] == TRUE)),]$state_info)))
  if(lineage %in% ds_names){
    clones_interest = c(clones_interest, c)
    sub = adata@meta.data[names(which(clone_bin[,c] == TRUE)),]
    if(2 %in% sub$time_info){ #make sure original cell starts at time 2
      sub2 = sub[sub$state_info %in% c('undiff',lineage),]
      cells_interest = c(cells_interest, rownames(sub2))
    }
  }
}
adata_sub = subset(adata, cells = cells_interest)

toplot = data.frame("marker"= adata_sub@assays$RNA@data[marker,], "sF_invert" = adata_sub$stemFinder_invert)
ggplot(toplot, aes(x = marker, y = sF_invert)) + geom_point()

cor.test(x = toplot$marker, y = toplot$sF_invert, method = 'pearson')
FeaturePlot(adata, features = c(marker, 'stemFinder_invert'))
