## Automated potency quantification


#Goal: robustness testing. Test various values of K for each dataset, 10 iterations for each value
#Only need Gini outputs (scaled zero, CC genes only)

library(Seurat)
library(pryr)
library(dplyr)
library(MASS)
library(ggplot2)
setwd("~/Documents/Stochasticity/")
load("~/Documents/Stochasticity/s_genes_mouse.rda")
load("~/Documents/Stochasticity/g2m_genes_mouse.rda")
load("~/Documents/Stochasticity/s_genes_human.rda")
load("~/Documents/Stochasticity/g2m_genes_human.rda")
load("~/Documents/Stochasticity/s_genes_celeg.rda")
load("~/Documents/Stochasticity/g2m_genes_celeg.rda")
setwd("~/Documents/Stochasticity")
set.seed(123)
list_all = list(spear_all_gini = c(), auc_gini = c(), spear_pheno_gini = c(), dataset = c(), iteration = c(), markerlist = c())
d=1 #set your counter once before running everything

#get info about each dataset
filenames = list.files(pattern = ".rds", full.names = T)

for (f in filenames[c(26:length(filenames))]){ 
  id = gsub(f, pattern = ".rds", replacement = "")
  id = gsub(id, pattern = "./", replacement = "")
  
  adata = readRDS(f)
  
  if (grepl("mouse", id, ignore.case = T) == T){
    org = 'mouse'
    s_genes = s_genes_mouse[s_genes_mouse %in% rownames(adata)]
    g2m_genes = g2m_genes_mouse[g2m_genes_mouse %in% rownames(adata)]
    g1_genes = g1_mouse[g1_mouse %in% rownames(adata)]
    g0_genes = go_mouse[go_mouse %in% rownames(adata)]
  } else if (grepl("hum", id, ignore.case = T) == T | grepl("zebra", id, ignore.case = T) == T){
    org = 'human'
    s_genes = s_genes_human[s_genes_human %in% rownames(adata)]
    g2m_genes = g2m_genes_human[g2m_genes_human %in% rownames(adata)]
    g1_genes = g1_human[g1_human %in% rownames(adata)]
    g0_genes = go_human[go_human %in% rownames(adata)]
  }else if (grepl("celeg", id, ignore.case = T) == T){
    org = 'celeg'
    s_genes = s_genes_celeg[s_genes_celeg %in% rownames(adata)]
    g2m_genes = s_genes_celeg[s_genes_celeg %in% rownames(adata)]
    g1_genes = g1_celeg[g1_celeg %in% rownames(adata)]
    g0_genes = go_celeg[go_celeg %in% rownames(adata)]
  }
  
  markers_std = c(s_genes, g2m_genes)
  markers_list = list(markers_std, s_genes, g2m_genes, g1_genes, g0_genes)
  names(markers_list) = c('Full_Regev','S','G2M','G1','G0')
  
  #VariableFeatures(adata) = VariableFeatures(adata)[!VariableFeatures(adata) %in% c(markers_std, g1_genes)]
  #pcs = 35 #not great that we have to use this for every dataset
  #adata = CellCycleScoring(adata, s.features = s_genes, g2m.features = g2m_genes)
  #adata = ScaleData(adata, features = rownames(adata), vars.to.regress = c('S.Score','G2M.Score'))
  #expDat = as.matrix(adata@assays$RNA@scale.data)
  
  #k = round(sqrt(ncol(adata)))
  #adata = FindNeighbors(adata, dims = 1:pcs, k.param = k)
  #knn = adata@graphs$RNA_nn
  
      #run Gini 
  for (m in 1:length(markers_list)){
    markers_name = names(markers_list[m])
    markers = markers_list[[m]]
    adata = gene_set_score(adata, gene.set = markers)

    if('Ground_truth' %in% colnames(adata@meta.data)){
      list_all = plot_and_store_gt(adata, d, id, list_all, markerlist = markers_name)

    }else {
      list_all = plot_and_store_order(adata, d, id, list_all, markerlist = markers_name)
    }
    d = d+1 
  }
  save(list_all, file = "list_all_test_CCgenesetscore_andcompetitors_inprog.rda")
}

######################

#Plotting
    #want: box plot for each k value, where points = iterations
    #y axis: performance value of choice, x axis = k values
    #mark performance at ideal_k with horizontal dotted line

df = data.frame(dataset = list_all$dataset, k = list_all$k, k_ratio = list_all$k_ratio, auc_gini = list_all$auc_gini, spear_all_gini = list_all$spear_all_gini, spear_pheno_gini = list_all$spear_pheno_gini, iteration = list_all$iteration)
df$k_ratio_round = round(df$k_ratio, digits = 2)
df$k_ratio_factor = as.factor(df$k_ratio)
#df$k_ratio_fraction = as.factor(as.character(fractions(list_all$k_ratio)))
df$k_ratio_isone = rep(0, nrow(df))
df[df$k_ratio == 1,]$k_ratio_isone = 1
df$UMI = rep('UMI-based data', nrow(df))
df[df$dataset %in% unique(df$dataset)[c(1,4,9,14,19,23,24,27,28,29,30)],]$UMI = 'Counts or TPM/FPKM-based data'

  #remove datasets that you're no longer using
toremove = which(df$dataset == unique(df$dataset)[33])
df = df[-toremove,]

df$neg_deviation_from_auc = rep(NA, nrow(df))
df$neg_deviation_from_spearss = rep(NA, nrow(df))
df$neg_deviation_from_spearp = rep(NA, nrow(df))

for (d in df$dataset){
  df[df$dataset == d,]$neg_deviation_from_auc = df[df$dataset == d,]$auc_gini - df[df$k_ratio == 1 & df$dataset == d, ]$auc_gini
  df[df$dataset == d,]$neg_deviation_from_spearss = df[df$dataset == d,]$spear_all_gini - df[df$k_ratio == 1 & df$dataset == d, ]$spear_all_gini
  df[df$dataset == d,]$neg_deviation_from_spearp = df[df$dataset == d,]$spear_pheno_gini - df[df$k_ratio == 1 & df$dataset == d, ]$spear_pheno_gini
}

  #separate plot for each dataset
xlabels = sort(unique(df$k_ratio_round))
#xlabels[seq(2, length(xlabels), 2)] <- '' #hide every other label
xlabels[!xlabels %in% c('0','0.25','0.5','1','1.5','2','2.5','3','3.5','4')] = ""
ggplot(df[df$dataset %in% datasets[5:6],], aes(x = k_ratio_factor, y = neg_deviation_from_auc)) + geom_boxplot(aes(color = k_ratio_isone)) + geom_point(position = position_dodge(width=0.75),aes(group = k_ratio_factor, color = k_ratio_isone)) + facet_wrap(~dataset) + theme_classic() + theme(legend.position = 'none') + scale_x_discrete(labels = xlabels) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.ticks = element_blank()) + xlab("K/Kideal") + ylab("-(Deviation from AUC at Kideal)") 

  #for all datasets
xlabels = sort(unique(df$k_ratio_factor))
xlabels[!xlabels %in% c('0','0.25','0.5','1','1.5','2','2.5','3','3.5','4')] = ""
ggplot(df, aes(x = k_ratio_factor, y = neg_deviation_from_auc)) + geom_boxplot() + geom_point(position = position_dodge(width=0.75),aes(group = k_ratio_factor, color = dataset)) + ggtitle("Robustness of Gini method to Changes in K") + scale_fill_discrete(name = "Dataset") + geom_hline(yintercept = 0, linetype = 'dashed', col = 'black') + scale_x_discrete(labels = xlabels) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.ticks = element_blank()) + xlab("K/Kideal") + ylab("-(Deviation from AUC at Kideal)")  + theme_classic()
  
  #for all datasets with bins
df$k_bins = rep("0 - 0.5", nrow(df))
df[df$k_ratio >= 0.5 & df$k_ratio < 1,]$k_bins = "0.5 - 1"
df[df$k_ratio >= 1 & df$k_ratio < 1.5,]$k_bins = "1 - 1.5"
df[df$k_ratio >= 1.5 & df$k_ratio < 2,]$k_bins = "1.5 - 2"
df[df$k_ratio >= 2 & df$k_ratio < 2.5,]$k_bins = "2 - 2.5"
df[df$k_ratio >= 2.5 & df$k_ratio < 3,]$k_bins = "2.5 - 3"
df[df$k_ratio >= 3 & df$k_ratio <= 3.5,]$k_bins = "3 - 3.5"
df[df$k_ratio >= 3.5 & df$k_ratio <= 4,]$k_bins = "3.5 - 4"
ggplot(df, aes(x = k_bins, y = neg_deviation_from_auc)) + geom_boxplot() + geom_point(position = position_dodge(width=0.75),aes(group = k_bins, color = UMI)) + theme_classic() + xlab("K/Kideal") + ylab("-(Deviation from AUC at Kideal)") + ggtitle("Robustness of Gini method to Changes in K") + scale_fill_discrete(name = "Dataset")  
#################
#################
#################
#################

#ANOVA one way: tells you if any groups are signif diff
one.way.auc = aov(auc_gini ~ k_bins, data = df)
summary(one.way.auc)

#Tukey's HSD: tells you which groups are signif diff
tukey.auc = TukeyHSD(one.way.auc)
tukey.auc 



####################
#FUNCTIONS ABOVE


####
## Compute Stochasticity using Gini index

run_gini <- function(adata, nn = knn, k = k, thresh = 0, markers = markers){
  
  gini_agg = data.frame("cell" = colnames(nn), "gini_index_agg" = rep(NA, ncol(nn)))
  
  for (cell in rownames(nn)){
    neigh = names(nn[cell,][nn[cell,]== 1]) #barcodes of the cell + its nearest neighbors
    exp = expDat[markers,neigh] > thresh #binarize gene expression matrix
    exp_match = exp == exp[,cell] #this will be T if expression pattern of a neighboring cell matches expression pattern of cell of interest, and vice versa
    n_match = apply(exp_match, 1, sum) - 1 #number of neighboring cells with the same expression pattern of gene g as the cell of interest
    p_g = n_match/(k - 1) #probability of expression pattern matching the cell of interest in this neighborhood: calculates separate p_g for each gene
    gini_g = p_g * (1 - p_g) #gini index for gene g for cell of interest
    
    gini_agg[gini_agg$cell == cell,]$gini_index_agg = sum(gini_g) #report mean gini index for each cell
  }
  adata@meta.data$gini_cellcycle = gini_agg$gini_index_agg #add to Seurat metadata
  adata@meta.data$gini_cellcycle_normandinvert = 1 - (adata@meta.data$gini_cellcycle)/max(adata@meta.data$gini_cellcycle)
  
  return(adata)
}


#########
######### store Spearman and AUC values in a list

plot_and_store_gt <- function(adata, d = d, id = id, list_all = list_all, markerlist = markers_name){
  
  #Spearman all cells
  spear_all_geneset = cor.test(x = adata$gene.set.score, y = adata$Ground_truth, method = "spearman", exact = F)$estimate
  spear_all_cyto = cor.test(x = adata$CytoTRACE_normandinvert, y = adata$Ground_truth, method = "spearman", exact = F)$estimate
  spear_all_ccat = cor.test(x = adata$ccat_normandinvert, y = adata$Ground_truth, method = "spearman", exact = F)$estimate
  
  #AUC
  mostpotent = adata@meta.data[adata$Ground_truth == max(adata$Ground_truth),]
  leastpotent = adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),]
  category = c(rep(1, nrow(mostpotent)), rep(0, nrow(leastpotent))) #gold standard
  
  prediction.geneset = c(mostpotent$gene.set.score, leastpotent$gene.set.score) #gene set score results
  auc_geneset = auc_probability(category, prediction.geneset)
  
  prediction.cyto = c(mostpotent$CytoTRACE_normandinvert, leastpotent$CytoTRACE_normandinvert) #calculator results
  auc_cyto = auc_probability(category, prediction.cyto)
  
  prediction.ccat = c(mostpotent$ccat_normandinvert, leastpotent$ccat_normandinvert) #calculator results
  auc_ccat = auc_probability(category, prediction.ccat)
  
  
  #Spearman pheno 
  clusters = as.character(unique(adata$Phenotype))
  meanpotency = data.frame("cluster" = clusters, "gene_set_score" = rep(NA, length(clusters)), "ccat_potency" = rep(NA, length(clusters)), "cyto_potency" = rep(NA, length(clusters)), "ground_truth" = rep(NA, length(clusters)))
  for (i in clusters){
    meanpotency$gene_set_score[which(clusters == i)] = mean(adata@meta.data$gene.set.score[adata$Phenotype == i])
    meanpotency$ground_truth[which(clusters == i)] = mean(adata@meta.data$Ground_truth[adata$Phenotype == i])
  }
  spear_pheno_geneset = cor.test(x = meanpotency$gene_set_score, y = meanpotency$ground_truth)$estimate
  
  #Store in list
  list_all$dataset[d] = id
  list_all$spear_all_geneset[d] = spear_all_geneset
  list_all$auc_geneset[d] = auc_geneset
  list_all$spear_pheno_geneset[d] = spear_pheno_geneset
  list_all$markerlist[d] = markers_name
  
  print(list_all)
  return(list_all)
}




plot_and_store_order <- function(adata, d = d, id = id, list_all = list_all, markerlist = markers_name){
  
  #Spearman all cells
  spear_all_gini = cor.test(x = adata$gini_cellcycle_normandinvert, y = adata$Order, method = "spearman", exact = F)$estimate
  print(paste("Spearman, all cells, Gini:", spear_all_gini, sep = " "))
  
  #AUC
  mostpotent = adata@meta.data[adata$Order == max(adata$Order),]
  leastpotent = adata@meta.data[adata$Order == min(adata$Order),]
  category = c(rep(1, nrow(mostpotent)), rep(0, nrow(leastpotent))) #gold standard
  prediction.gini = c(mostpotent$gini_cellcycle_normandinvert, leastpotent$gini_cellcycle_normandinvert) #calculator results
  auc_gini = auc_probability(category, prediction.gini)
  
  #Spearman pheno 
  clusters = as.character(unique(adata$Phenotype))
  meanpotency = data.frame("cluster" = clusters, "gini_potency" = rep(NA, length(clusters)), "ccat_potency" = rep(NA, length(clusters)), "cyto_potency" = rep(NA, length(clusters)), "ground_truth" = rep(NA, length(clusters)))
  for (i in clusters){
    meanpotency$gini_potency[which(clusters == i)] = mean(adata@meta.data$gini_cellcycle_normandinvert[adata$Phenotype == i])
    meanpotency$ground_truth[which(clusters == i)] = mean(adata@meta.data$Order[adata$Phenotype == i])
  }
  spear_pheno_gini = cor.test(x = meanpotency$gini_potency, y = meanpotency$ground_truth)$estimate
  
  #Store in list
  list_all$dataset[d] = id
  list_all$spear_all_gini[d] = spear_all_gini
  list_all$auc_gini[d] = auc_gini
  list_all$spear_pheno_gini[d] = spear_pheno_gini
  list_all$markerlist[d] = markers_name
  
  print(list_all)
  return(list_all)
}


#AUC probability
auc_probability <- function(labels, scores, N=1e7){
  pos <- sample(scores[labels], N, replace=TRUE)
  neg <- sample(scores[!labels], N, replace=TRUE)
  # sum( (1 + sign(pos - neg))/2)/N # does the same thing
  (sum(pos > neg) + sum(pos == neg)/2) / N # give partial credit for ties
}
