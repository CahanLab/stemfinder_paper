### Create a Seurat object from raw counts data - example 
### K Noller

# Setup & load raw data
library(Seurat)

load("s_genes_mouse.rda")
load("g2m_genes_mouse.rda")
s_genes = s_genes_mouse
g2m_genes = g2m_genes_mouse

raw = Read10X("~/Downloads/GSE147729_RAW/oldB/")

# Create Seurat object
adata = CreateSeuratObject(counts = raw, min.cells = 5, min.features = 200, project = "GSE166504_liver")
adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^mt-")
adata[["percent.ribo"]] <- PercentageFeatureSet(adata, pattern = "^Rp[sl]")
VlnPlot(adata, features=c("nCount_RNA", "percent.mt", "percent.ribo"), ncol =3)

# Remove mito, ribo, MALAT1 genes
ribogenes <- grep("^Rp[sl]*", rownames(adata), value=TRUE)
mtgenes <- grep("^mt--*", rownames(adata), value=TRUE)
length(ribogenes)
length(mtgenes)
adata <- adata[!(rownames(adata) %in% c(ribogenes, mtgenes, "Malat1")),]

# Filter by max nFeature and percent.mt
fmax <- quantile(adata$nFeature_RNA, 0.95)
mtmax <- 5
adata <- subset(adata, subset = nFeature_RNA > 200 & nFeature_RNA < fmax & percent.mt < mtmax)

# Normalize, compute HVG
adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2500)
VariableFeatures(adata) = VariableFeatures(adata)[!VariableFeatures(adata) %in% c(s_genes, g2m_genes)] #for stemFinder only: exclude CC genes from HVG

saveRDS(adata, file = "adata_qc.rds")
write.table(hvg, file = "hvg.txt", row.names=F, sep = '\t')

# Scale, CC scoring
adata <- ScaleData(adata, features = rownames(adata))

s_genes = s_genes[s_genes %in% rownames(adata)]
length(s_genes)
g2m_genes = g2m_genes[g2m_genes %in% rownames(adata)]
length(g2m_genes)
adata <- CellCycleScoring(adata, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)

# PCA
adata <- RunPCA(adata) 
ElbowPlot(adata, ndims=50) #use this to choose number of principal components for KNN 
DimHeatmap(adata, dims = 30:40, cells = 500, balanced = TRUE)
pcs <- 30

# KNN 
k <- round(sqrt(ncol(adata))) #standard formula for k 
adata <- FindNeighbors(dims = 1:pcs, k.param = k) 
