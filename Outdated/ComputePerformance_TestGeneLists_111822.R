## Automated potency quantification


#Setup 

library(Seurat)
library(dplyr)
setwd("~/Documents/ValidationandGeneLists/")

load("s_genes_mouse.rda")
load("g2m_genes_mouse.rda")
load("s_genes_human.rda")
load("g2m_genes_human.rda")
load("s_genes_celeg.rda")
load("g2m_genes_celeg.rda")
load("G1_GSEAMsigDB_celeg.rda")
load("G1_GSEAMsigDB_human.rda")
load("G1_GSEAMsigDB_mouse.rda")
load("GO_cellcycle_celeg.rda")
load("GO_cellcycle_mouse.rda")
load("GO_cellcycle.rda")
load("mmTFs.rda")
load("hsTFs.rda")
load("ceTFs.rda")

load("pcs_all.rda") 
    ##Numeric vector of the number of PCs that should be inputted to FindNeighbors() 
    ##Values are chosen based on elbow plot inflection point for each dataset

set.seed(123)
d=1 

list_all = list(spear_all_gini = c(), auc_gini = c(), spear_pheno_gini = c(), dataset = c(), iteration = c(), genelist = c(), regression = c(), pct.recov_gini = c(), pct.recov_cyto = c(), pct.recov_ccat = c(), auc_cyto = c(), auc_ccat = c(), spear_all_cyto = c(), spear_all_ccat = c(), spear_pheno_cyto = c(), spear_pheno_ccat = c())


#Load validation datasets, select marker genes, and run stemFinder

filenames = list.files(pattern = ".rds", full.names = T)

for (f in filenames[1:(length(filenames))]){
  
  id = gsub(f, pattern = ".rds", replacement = "")
  id = gsub(id, pattern = "./", replacement = "")
  
  adata = readRDS(f)
  
  if (grepl("mouse", id, ignore.case = T) == T){
    org = 'mouse'
    s_genes = s_genes_mouse[s_genes_mouse %in% rownames(adata)]
    g2m_genes = g2m_genes_mouse[g2m_genes_mouse %in% rownames(adata)]
    g0_genes = go_cc_mouse[go_cc_mouse %in% rownames(adata)]
    TFs = mmTFs
  } else if (grepl("hum", id, ignore.case = T) == T | grepl("zebra", id, ignore.case = T) == T){
    org = 'human'
    s_genes = s_genes_human[s_genes_human %in% rownames(adata)]
    g2m_genes = g2m_genes_human[g2m_genes_human %in% rownames(adata)]
    g0_genes = go_cc[go_cc %in% rownames(adata)]
    TFs = hsTFs
  }else if (grepl("celeg", id, ignore.case = T) == T){
    org = 'celeg'
    s_genes = s_genes_celeg[s_genes_celeg %in% rownames(adata)]
    g2m_genes = s_genes_celeg[s_genes_celeg %in% rownames(adata)]
    g0_genes = go_celeg[go_celeg %in% rownames(adata)]
    TFs = ceTFs
  }
  
  #Compute sporadically expressed TFs
  finals_withcc = find_sparse_tfs(adata, TFs)
  
  #Compile all gene lists 
  full_regev = c(s_genes, g2m_genes)
  markers_all = c(list(full_regev, g0_genes), finals_withcc)
  names(markers_all) = c('Full_Regev','G0','CCandTF_1','CCandTF_2','CCandTF_3','CCandTF_4','CCandTF_5','SporadTFonly')
    
  #Set input parameters
  k = round(sqrt(ncol(adata)))
  pcs = pcs_all[which(filenames==f)] 
  
  #Run stemFinder for n = 5 iterations per gene list
  
  niter = 5
  
  for(m in 1:length(markers_all)){
    markers = markers_all[[m]]
    genelist_name = names(markers_all)[m]
    
    for(iter in 1:niter){
      
      #Compute HVG, scale, PCA, KNN 
      adata = FindVariableFeatures(adata, selection.method = 'vst', nfeatures = 2500)
      VariableFeatures(adata) = VariableFeatures(adata)[!VariableFeatures(adata) %in% markers]
      adata = ScaleData(adata, features = rownames(adata))
      adata = RunPCA(adata)
      adata = FindNeighbors(adata, dims = 1:pcs, k.param = k)
      knn = adata@graphs$RNA_nn
      
      #Run stemFinder
      adata = run_stemFinder(adata, k = k, nn = knn, thresh = 0, markers = markers)
      print("stemFinder finished running successfully!")
      
      #Save anndata files for standard cell cycle gene list iterations
      if(genelist_name == 'Full_Regev' & iter == 1){ 
        saveRDS(adata, file = f)
      }
      
      #Compute performance
      if('Ground_truth' %in% colnames(adata@meta.data)){
        list_all = compute_performance(adata, d, id, list_all, iter = iter, genelist = genelist_name, regression = "Not_Regressed")
      }else {
        adata$Ground_truth = adata$Order
        list_all = compute_performance(adata, d, id, list_all, iter = iter, genelist = genelist_name, regression = "Not_Regressed")
      }
      d = d+1
    }
  }
  save(list_all, file = "list_all_automatedperformance.rda")
} 
