## Automated potency quantification & test of various gene lists in stemFinder

  #wd_string is a string describing the location of the folder containing all query datasets and gene lists
        ##validation dataset filename must contain string "mouse," "hum," "zebra," or "celeg" corresponding to its species
  #pcs_all is a numeric vector of the number of PCs that should be inputted to FindNeighbors(). Values are chosen based on elbow plot inflection point for each dataset
  #niter is the number of iterations to process validation data, run stemFinder and compute performance 

  #requires: Seurat, dplyr, matrixStats

automated_performance_genelists <- function(wd_string = "/home/ubuntu/data/ValidationandGeneLists", pcs_all = pcs_all, niter = 3){

  library(Seurat)
  library(dplyr)
  library(matrixStats)
  setwd(wd_string)
  
  #Regev lab CC gene lists
  load("s_genes_mouse.rda")
  load("g2m_genes_mouse.rda")
  load("s_genes_human.rda")
  load("g2m_genes_human.rda")
  load("s_genes_celeg.rda")
  load("g2m_genes_celeg.rda")
  
  #Gene ontology CC gene lists
  load("GO_cellcycle.rda")
  load("GO_cellcycle_mouse.rda")
  load("GO_cellcycle_celeg.rda")
  
  #KEGG CC gene lists
  load("kegg_cellcyclegenes_human.rda")
  load("kegg_cellcyclegenes_mouse.rda")
  load("kegg_cellcyclegenes_celeg.rda")

  set.seed(123)
  d=1 

  list_all = list(dataset = c(), iteration = c(), genelist = c(), spear_all_sF = c(), auc_sF = c(), spear_pheno_sF = c(), pct.recov_sF = c())

  filenames = list.files(pattern = ".rds", full.names = T)
  
  for (f in filenames[1:(length(filenames))]){
    id = gsub(f, pattern = ".rds", replacement = "")
    id = gsub(id, pattern = "./", replacement = "")
    adata = readRDS(f)
    
    # Cell cycle gene lists
    if (grepl("mouse", id, ignore.case = T) == T){
      s_genes = s_genes_mouse[s_genes_mouse %in% rownames(adata)]
      g2m_genes = g2m_genes_mouse[g2m_genes_mouse %in% rownames(adata)]
      GO_genes = go_cc_mouse[go_cc_mouse %in% rownames(adata)]
      KEGG_genes = keggcc_mouse[keggcc_mouse %in% rownames(adata)]
    } else if (grepl("hum", id, ignore.case = T) == T | grepl("zebra", id, ignore.case = T) == T){
      s_genes = s_genes_human[s_genes_human %in% rownames(adata)]
      g2m_genes = g2m_genes_human[g2m_genes_human %in% rownames(adata)]
      GO_genes = go_cc[go_cc %in% rownames(adata)]
      KEGG_genes = kegggenes[kegggenes %in% rownames(adata)]
    }else if (grepl("celeg", id, ignore.case = T) == T){
      s_genes = s_genes_celeg[s_genes_celeg %in% rownames(adata)]
      g2m_genes = s_genes_celeg[s_genes_celeg %in% rownames(adata)]
      GO_genes = go_celeg[go_celeg %in% rownames(adata)]
      KEGG_genes = keggcc_celeg[keggcc_celeg %in% rownames(adata)]
    }
    
    full_regev = c(s_genes, g2m_genes)
    
    # Random gene lists
    len = length(full_regev) # same length as default CC Gene list
    m_and_v <- data.frame(mean = rowMeans(x = as.matrix(adata@assays$RNA@data)), variance = rowVars(x = as.matrix(adata@assays$RNA@data)), row.names = rownames(adata)) 
    m_and_v_cc <- m_and_v[full_regev,]
    m_and_v_other <- m_and_v[!rownames(m_and_v) %in% full_regev,]
    med_mean <- median(m_and_v_cc$mean)
    med_var <- median(m_and_v_cc$variance)
        ##bin genes by whether mean & variance are < or >= median CC gene expression 
    lowm_lowv <- nrow(m_and_v_cc[m_and_v_cc$mean < med_mean & m_and_v_cc$variance < med_var,])
    lowm_highv <- nrow(m_and_v_cc[m_and_v_cc$mean < med_mean & m_and_v_cc$variance >= med_var,])
    highm_lowv <- nrow(m_and_v_cc[m_and_v_cc$mean >= med_mean & m_and_v_cc$variance < med_var,])
    highm_highv <- nrow(m_and_v_cc[m_and_v_cc$mean >= med_mean & m_and_v_cc$variance >= med_var,])
   
        ##select random gene lists with similar binned mean & variances as CC list
    markers_all = list(full_regev, GO_genes, KEGG_genes)
    markers_rand = c()
    for(rep in 1:5){
      if(nrow(m_and_v_other[m_and_v_other$mean < med_mean & m_and_v_other$mean < med_var,]) > 0){
        rand_lowm_lowv <- sample(x = rownames(m_and_v_other[m_and_v_other$mean < med_mean & m_and_v_other$mean < med_var,]), size = lowm_lowv, replace = F)
        markers_rand = c(markers_rand, rand_lowm_lowv)
      }
      if(nrow(m_and_v_other[m_and_v_other$mean < med_mean & m_and_v_other$mean >= med_var,]) > 0){
        rand_lowm_highv <- sample(x = rownames(m_and_v_other[m_and_v_other$mean < med_mean & m_and_v_other$mean >= med_var,]), size = lowm_highv, replace = F)
        markers_rand = c(markers_rand, rand_lowm_highv)
      }
      if(nrow(m_and_v_other[m_and_v_other$mean >= med_mean & m_and_v_other$mean < med_var,]) > 0){
        rand_highm_lowv <- sample(x = rownames(m_and_v_other[m_and_v_other$mean >= med_mean & m_and_v_other$mean < med_var,]), size = highm_lowv, replace = F)
        markers_rand = c(markers_rand, rand_highm_lowv)
      }
      if(nrow(m_and_v_other[m_and_v_other$mean >= med_mean & m_and_v_other$mean >= med_var,]) > 0){
        rand_highm_highv <- sample(x = rownames(m_and_v_other[m_and_v_other$mean >= med_mean & m_and_v_other$mean >= med_var,]), size = highm_highv, replace = F)
        markers_rand = c(markers_rand, rand_highm_highv)
      }
      markers_all[[rep + 3]] = markers_rand
      markers_rand = c()
    }
    names(markers_all) = c('Full_Regev','GeneOntology','KEGG',rep('Random', 5))
    
    #Select other input parameters, run HVG
    k = round(sqrt(ncol(adata)))
    pcs = pcs_all[which(filenames==f)] # from elbow plots examined beforehand
    adata = FindVariableFeatures(adata, selection.method = 'vst', nfeatures = 2500)
    
    #Adjust HVG, scale, PCA, KNN and run stemFinder 
    for(m in 1:length(markers_all)){
      markers = markers_all[[m]]
      genelist_name = names(markers_all)[m]
      
      for(iter in 1:niter){
        
        VariableFeatures(adata) = VariableFeatures(adata)[!VariableFeatures(adata) %in% markers]
        adata = ScaleData(adata, features = rownames(adata))
        adata = RunPCA(adata)
        adata = FindNeighbors(adata, dims = 1:pcs, k.param = k)
        knn = adata@graphs$RNA_nn
        
        adata = run_stemFinder(adata, k = k, nn = knn, thresh = 0, markers = markers)
        print("stemFinder finished running successfully!")
        
        if(genelist_name == 'Full_Regev' & iter == 1){ 
          saveRDS(adata, file = f)
        }
        
        if(!('Ground_truth' %in% colnames(adata@meta.data))){
          adata$Ground_truth = adata$Order
        }
        list_all = compute_performance_genelists(adata, d, id, list_all, iter = iter, genelist = genelist_name)
        
        d = d+1
      }
      save(list_all, file = "list_all_genelist.rda")
    }
  }
  return(list_all)
}

