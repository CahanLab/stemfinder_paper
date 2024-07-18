## Robustness testing - Downsampling of potent population OR proportional downsampling of all celltypes

#wd_string is a string describing the location of the folder containing all query datasets and gene lists
##validation dataset filename must contain string "mouse," "hum," "zebra," or "celeg" corresponding to its species
##gene lists must include S and G2M lists for each species in validation datasets
#pcs_all is a numeric vector of the number of PCs that should be inputted to FindNeighbors(). Values are chosen based on elbow plot inflection point for each dataset
#niter is the number of iterations to process validation data, run stemFinder and compute performance 

#requires: Seurat, dplyr

automated_performance_ds <- function(wd_string = "~/data", pcs_all = pcs_all, niter = 3){
    
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
    
    list_all = list(dataset = c(), iteration = c(), ds_ratio = c(), ds_type = c(), spear_all_sF = c(), auc_sF = c(), spear_pheno_sF = c(), pct.recov_sF = c())
    
    filenames = list.files(pattern = ".rds", full.names = T)
    
    for (f in 1:(length(filenames))){
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
        
        if(!('Ground_truth' %in% colnames(adata@meta.data))){ #rename metadata column
            adata$Ground_truth = adata$Order
        }
        
        #Default inputs of CC genes, PCs
        markers = c(s_genes, g2m_genes)
        pcs = pcs_all[which(filenames==f)] #PCs from previously examined elbow plot
        
        #Downsampling 
        ds_ratio = rep(seq(0.1, 1, by = 0.1),2)
        names(ds_ratio) = c(rep("All cell types", 10), rep("Highly potent only", 10))
        
        for (iter in 1:niter){
            for (r in 1:length(ds_ratio)){
                ratio = ds_ratio[[r]]
                type = names(ds_ratio[r]) 
                
                #Proportional downsampling of all cell types
                if (type == 'All cell types'){
                    num_all = round(table(adata$Phenotype)*ratio)
                    cells_sub = c()
                    for(l in 1:length(num_all)){
                        cells_sub = c(cells_sub, sample(rownames(adata@meta.data[adata$Phenotype == names(num_all[l]),]), size = num_all[[l]], replace = F)) #barcodes of cells in a given phenotype
                    }
                    adata_sub = subset(adata, cells = cells_sub)}
                
                #Downsampling of highly potent population only, as defined by ground truth potency
                else{cells_potent = rownames(adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),])
                num_potent = length(cells_potent)
                cells_sub = c(rownames(adata@meta.data)[!rownames(adata@meta.data) %in% cells_potent], sample(cells_potent, round(num_potent * ratio), replace = F))
                adata_sub = subset(adata, cells = cells_sub)}
                
                #HVG, scale, PCA, KNN
                adata_sub = FindVariableFeatures(adata_sub, selection.method = 'vst', nfeatures = 2500)
                VariableFeatures(adata_sub) = VariableFeatures(adata_sub)[!VariableFeatures(adata_sub) %in% markers]
                adata_sub = ScaleData(adata_sub, features = rownames(adata_sub))
                if (ncol(adata_sub) < 50){ #ensure PCs not greater than number of cells in dataset
                    adata_sub = RunPCA(adata_sub, npcs = (ncol(adata_sub) - 1))
                }else {adata_sub = RunPCA(adata_sub, npcs = 50)}

                k = round(sqrt(ncol(adata_sub))) #default value 
                if (pcs > ncol(adata_sub)){
                    pcs = ncol(adata_sub) - 1
                }
                adata_sub = FindNeighbors(adata_sub, dims = 1:pcs, k.param = k)
                knn = adata_sub@graphs$RNA_nn
                
                #Run stemFinder and compute performance
                adata_sub = run_stemFinder(adata_sub, k = k, nn = knn, thresh = 0, markers = markers)
                print("stemFinder finished running successfully!")
                
                list_all = compute_performance_ds(adata_sub, d, id, list_all, iter = iter, ds_ratio = ratio, ds_type = type)
                d = d+1
            }
            save(list_all, file = "list_all_downsample.rda")
        }
    }
return(list_all)
}
