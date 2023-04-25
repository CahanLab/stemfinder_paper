### Find TFs in dataset with high stochasticity in one phenotype-defined cluster & low stochasticity in another 
### TFs overlap with top DE genes
### Inputs:
  # TFs is list of species-specific TF names
  # lowthresh_tf and highthresh_tf are quantile-based thresholds for stochasticity-based TF selection
  # min.pct and logfc.threshold are for DE analysis 
  # cell_cycle_genes is list of S and G2M phase genes

find_sparse_tfs <- function(adata, TFs, lowthresh_tf = 0.4, highthresh_tf = 0.6, min.pct = 0.3, logfc.threshold = 0.5, cell_cycle_genes = c(s_genes, g2m_genes)){
  
  #Retrieve expression data and metadata
  markers = TFs[TFs %in% rownames(adata)]
  expDat = as.matrix(adata@assays$RNA@scale.data)
  sampTab = adata@meta.data
  idents = names(which(table(adata$Phenotype) != 1)) #do not select phenotype categories with only 1 cell
  
  #Compute Gini index of all TFs in dataset
  gini_agg = data.frame("cluster" = rep(idents, length(markers)), "gini_sum" = rep(NA, length(idents)*length(markers)), "TF" = rep("", length(idents)*length(markers)))
  
  for (i in idents){
    gini_agg[gini_agg$cluster ==i,]$TF = markers
    
    neigh = rownames(sampTab[sampTab$Phenotype == i,]) 
    exp = expDat[markers,neigh] > 0 
    
    gini_allsingles = c()
    
    for (c in neigh){
      gini_singles = data.frame("cell" = rep(c, length(markers)), "gini" = rep(NA, length(markers)), "TF" = markers)
      
      exp_match = exp == exp[,c] 
      n_match = apply(exp_match, 1, sum) - 1
      p_g = n_match/(length(neigh)-1) 
      gini_g = p_g * (1 - p_g) 
      
      gini_singles$gini = gini_g
      gini_allsingles = rbind(gini_allsingles, gini_singles)
    }
    
    for (m in markers){
      gini_agg[gini_agg$cluster == i & gini_agg$TF == m,]$gini_sum = sum(gini_allsingles[gini_allsingles$TF == m,]$gini)
    }
  }
  
  #Select TFs with high stochasticity in one cluster & low stochasticity of expression in another
  finals_down = gini_agg[gini_agg$gini_sum < quantile(gini_agg$gini_sum, lowthresh_tf),]$TF
  finals_up = gini_agg[gini_agg$gini_sum > quantile(gini_agg$gini_sum, highthresh_tf),]$TF
  finals_stoch = unique(finals_down[finals_down %in% finals_up])
  
  if(length(finals_stoch)==0){
    stop(print("No sporadically expressed TFs found with given thresholds. Please select more permissive thresholds."))
  }
  
  #Select TFs in top DE genes
  Idents(adata) = 'Phenotype'
  markers_de = FindAllMarkers(adata, min.pct = min.pct, logfc.threshold = logfc.threshold, only.pos=T, test.use="wilcox")
  markers_de = markers_de[markers_de$p_val_adj < 0.05,]$gene
  finals = finals_stoch[finals_stoch %in% markers_de]
  print("Sporadically expressed TFs overlapping with DE genes:")
  print(finals)
  if(length(finals) == 0 & length(finals_stoch > 0)){
    stop(print("No overlapping TFs found with given thresholds for differential gene expression. Please select more permissive thresholds."))
  }
  
  #Randomly sample cell cycle genes of equal length, n = 5 iterations
  finals_withcc = list()
  
  for (iter in seq(1:5)){
    if (length(finals) <= length(cell_cycle_genes)){
      cc_random = sample(cell_cycle_genes, size = length(finals), replace = F)
    }else{
      cc_random = cell_cycle_genes
    }
    finals_withcc[[iter]] = c(finals, cc_random)
  }
  
  finals_withcc[[6]] = finals 
  
  return(finals_withcc)
}


