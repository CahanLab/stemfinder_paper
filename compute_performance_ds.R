compute_performance_ds <- function(adata, d = d, id = id, list_all = list_all, iter = iter, ds_ratio, ds_type){
  
  #Single-cell Spearman correlation
  spear_all_gini = cor.test(x = adata$stemFinder_invert, y = adata$Ground_truth, method = "spearman", exact = F)$estimate

  #AUC
  mostpotent = adata@meta.data[adata$Ground_truth == max(adata$Ground_truth),]
  leastpotent = adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),]
  category = c(rep(1, nrow(mostpotent)), rep(0, nrow(leastpotent))) #gold standard
  
  prediction.gini = c(mostpotent$stemFinder_invert, leastpotent$stemFinder_invert) #calculator results
  auc_gini = auc_probability(category, prediction.gini)

  #Phenotypic Spearman correlation
  clusters = as.character(unique(adata$Phenotype))
  meanpotency = data.frame("cluster" = clusters, "gini_potency" = rep(NA, length(clusters)), "ccat_potency" = rep(NA, length(clusters)), "cyto_potency" = rep(NA, length(clusters)), "ground_truth" = rep(NA, length(clusters)))
  for (i in clusters){
    meanpotency$gini_potency[which(clusters == i)] = mean(adata@meta.data$stemFinder_invert[adata$Phenotype == i])
    meanpotency$ground_truth[which(clusters == i)] = mean(adata@meta.data$Ground_truth[adata$Phenotype == i])
  }
  spear_pheno_gini = cor.test(x = meanpotency$gini_potency, y = meanpotency$ground_truth)$estimate

  #Percent highly potent cell recovery
  num.stem.gt = nrow(adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),])
  pct.recov_gini = (nrow(adata@meta.data[adata$stemFinder_invert < quantile(adata$stemFinder_invert, (1/length(unique(adata$Ground_truth)))) & adata$Ground_truth == min(adata$Ground_truth),])/num.stem.gt) * 100

  #Store in list
  list_all$dataset[d] = id
  list_all$iteration[d] = iter
  list_all$ds_ratio[d] = ds_ratio
  list_all$ds_type[d] = ds_type

  list_all$spear_all_sF[d] = spear_all_gini
  list_all$auc_sF[d] = auc_gini
  list_all$spear_pheno_sF[d] = spear_pheno_gini
  list_all$pct.recov_sF[d] = pct.recov_gini
  
  return(list_all)
}

