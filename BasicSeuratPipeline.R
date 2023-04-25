#Basics of creating a Seurat object


raw = Read10X("~/Downloads/GSE147729_RAW/oldB/")
adata = CreateSeuratObject(counts = raw, min.cells = 5, min.features = 200, project = "GSE166504_liver")
adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^mt-")
adata[["percent.ribo"]] <- PercentageFeatureSet(adata, pattern = "^Rp[sl]")

VlnPlot(adata, features=c("nCount_RNA", "percent.mt", "percent.ribo"), ncol =3)

#Remove mito, ribo, MALAT1 genes
ribogenes <- grep("^Rp[sl]*", rownames(adata), value=TRUE)
mtgenes <- grep("^mt--*", rownames(adata), value=TRUE)
length(ribogenes)
length(mtgenes)

adata <- adata[!(rownames(adata) %in% c(ribogenes, mtgenes, "Malat1")),]
fmax <- quantile(adata$nFeature_RNA, 0.95)
mtmax <- 5
adata <- subset(adata, subset = nFeature_RNA > 200 & nFeature_RNA < fmax & percent.mt < mtmax)

adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2500)

saveRDS(adata, file = "~/Downloads/GSE147729_RAW/oldA/aginghsc_oldB_qc.rds")

### Process all samples together

#Merge and harmonize

#merge files into a single Seurat object
adata1 <- merge(li3, li6, add.cell.ids=c("LI3", "LI6"), project="Wang", merge.data = TRUE)
adata <- merge(adata8, li21, add.cell.ids=c("","LI21"), project="Wang", merge.data = TRUE)

#get highly variable genes for all samples and for the group together
adata.concat <- c(li3, li6, li9, li10, li11, li12, li13, li14, li15, li21)
adata <- FindVariableFeatures(adata, selection.method = "vst",
                              nfeatures = 2000)
hvgindiv <- unique(as.character(sapply(adata.concat, VariableFeatures)))
length(hvgindiv)
hvg <- unique(c(hvgindiv, VariableFeatures(adata))) #### check this is okay
length(hvg)
write.table(hvg, file = "hvg_li_092121.txt", row.names=F, sep = '\t')

# Scale data
adata <- ScaleData(adata, features = rownames(adata))

# Cell cycle genes
cell_cycle_genes = read.table("../regev_lab_cell_cycle_genes.txt")$V1
s_genes = cell_cycle_genes[1:43]
g2m_genes = cell_cycle_genes[43:length(cell_cycle_genes)]
s_genes = s_genes[s_genes %in% rownames(adata)]
length(s_genes)
g2m_genes = g2m_genes[g2m_genes %in% rownames(adata)]
length(g2m_genes)

adata <- CellCycleScoring(adata, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)


rm(adata.concat)

#PCA
adata <- RunPCA(adata) #, features = hvg)
ElbowPlot(adata, ndims=50)
DimHeatmap(adata, dims = 30:40, cells = 500, balanced = TRUE)
N <- 30

#Harmonize
adata <- adata %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, assay.use = "RNA", reduction.save = "harmony", reduction = "pca", dims.use=1:N)

#Downstream analysis - UMAP, neighbors, clusters
adata <- adata %>% 
  RunUMAP(reduction = "harmony", dims = 1:N) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:N, k.param = round(sqrt(ncol(adata)))) %>% 
  FindClusters(resolution = 0.5) 
