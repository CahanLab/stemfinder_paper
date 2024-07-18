## Robustness testing - K input parameter to KNN 

#wd_string is a string describing the location of the folder containing all query datasets and gene lists
##validation dataset filename must contain string "mouse," "hum," or "celeg" corresponding to its species
##gene lists must include S and G2M lists for each species in validation datasets
#pcs_all is a numeric vector of the number of PCs that should be inputted to FindNeighbors(). Values are chosen based on elbow plot inflection point for each dataset
#niter is the number of iterations to process validation data, run stemFinder and compute performance 

#requires: Seurat, dplyr

automated_performance_K <- function(wd_string = "~/data", pcs_all = pcs_all, niter = 3){
    
    library(Seurat)
    library(dplyr)
    setwd(wd_string)

    load("s_genes_mouse.rda")
    load("g2m_genes_mouse.rda")
    load("s_genes_human.rda")
    load("g2m_genes_human.rda")
    load("s_genes_celeg.rda")
    load("g2m_genes_celeg.rda")
    
    set.seed(123)
    d=1 
    
    list_all = list(dataset = c(), iteration = c(), k = c(), kideal = c(), kratio = c(), spear_all_sF = c(), auc_sF = c(), spear_pheno_sF = c(), pct.recov_sF = c())
    
    filenames = list.files(pattern = ".rds", full.names = T)
    
    for (f in filenames[1:(length(filenames))]){
        id = gsub(f, pattern = ".rds", replacement = "")
        id = gsub(id, pattern = "./", replacement = "")
        adata = readRDS(f)
        
        if (grepl("mouse", id, ignore.case = T) == T){
            s_genes = s_genes_mouse[s_genes_mouse %in% rownames(adata)]
            g2m_genes = g2m_genes_mouse[g2m_genes_mouse %in% rownames(adata)]
        } else if (grepl("hum", id, ignore.case = T) == T | grepl("zebra", id, ignore.case = T) == T){
            s_genes = s_genes_human[s_genes_human %in% rownames(adata)]
            g2m_genes = g2m_genes_human[g2m_genes_human %in% rownames(adata)]
        }else if (grepl("celeg", id, ignore.case = T) == T){
            s_genes = s_genes_celeg[s_genes_celeg %in% rownames(adata)]
            g2m_genes = s_genes_celeg[s_genes_celeg %in% rownames(adata)]
        }
        #Default inputs of CC genes and PCs
        markers = c(s_genes, g2m_genes)
        pcs = pcs_all[which(filenames==f)] #PCs from previously examined elbow plot

        #Choose K values ranging from K/Kideal = 0.25 to 4
        k_ideal = round(sqrt(ncol(adata))) 
        k_ratio = seq(0.25, 4, by = 0.25)
        k_values = round(k_ratio * k_ideal)
        
        #HVG, scale, PCA
        adata = FindVariableFeatures(adata, selection.method = 'vst', nfeatures = 2500)
        VariableFeatures(adata) = VariableFeatures(adata)[!VariableFeatures(adata) %in% markers]
        adata = ScaleData(adata, features = rownames(adata))
        adata = RunPCA(adata)
        
        #Compute KNN and run stemFinder across n = niter iterations per K
        for(k in 1:length(k_values)){
            for(iter in 1:niter){
                
                adata = FindNeighbors(adata, dims = 1:pcs, k.param = k_values[k])
                knn = adata@graphs$RNA_nn
                
                adata = run_stemFinder(adata, k = k_values[k], nn = knn, thresh = 0, markers = markers)
                print("stemFinder finished running successfully!")
                print(paste("Iteration: ", iter))
                
                if(!('Ground_truth' %in% colnames(adata@meta.data))){
                    adata$Ground_truth = adata$Order
                }
                
                list_all = compute_performance_K(adata, d, id, list_all, iter = iter, k = k_values[k], k_ratio = k_ratio[k], k_ideal = k_ideal)
                d = d+1
            }
            save(list_all, file = "list_all_k.rda")
        }
    }
    return(list_all)
}

