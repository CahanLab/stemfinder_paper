### Code for stemFinder manuscript figures

# updated 7/17/24 K.Noller

#Setup
library(stemFinder)
library(ggthemes)
library(RColorBrewer)
library(forcats)
library(stringr)
library(stringi)
library(rstatix)
library(ggpubr)
library(ggrepel)
library(ggthemes)
library(rstatix)
setwd("~/Dropbox (CahanLab)/Kathleen.Noller/Stochasticity/OtherCode/FiguresCode")

#Example code for saving as TIFF 300dpi figures
#tiff("test.tiff", units="in", width=6, height=5, res=300)
##insert ggplot code
#dev.off()

# To download datasets and objects needed below from S3, type the following URL with desired filename into browser: 
# https://cnobjects.s3.amazonaws.com/stemFinder/FiguresCode/filename.ext
# ex: https://cnobjects.s3.amazonaws.com/stemFinder/FiguresCode/compute_performance_K.R 
# file will automatically download

########################################

# Fig 1 - graphical abstract: majority generated with BioRender

    ##1C: Feature plots of BMMC 10X Tabula Muris, scatter plot for Spearman and pct.recov
adata = readRDS("MurineBoneMarrow10X_GSE109774.rds")
cell_cycle_genes = c(s_genes_mouse, g2m_genes_mouse)[c(s_genes_mouse, g2m_genes_mouse) %in% rownames(adata)]
VariableFeatures(adata) = VariableFeatures(adata)[!VariableFeatures(adata) %in% cell_cycle_genes]
adata = RunPCA(adata)
k = round(sqrt(ncol(adata)))
adata = FindNeighbors(adata, k.param = k, dims = 1:32) 
knn = adata@graphs$RNA_nn
adata = run_stemFinder(adata, k = k, nn = knn, thresh = 0, markers = cell_cycle_genes)

figc.1 = FeaturePlot(adata, features = 'stemFinder', cols = c('blue','red')) + ggtitle("stemFinder score") + xlab("UMAP_1") + ylab("UMAP_2")
figc.2 = FeaturePlot(adata, features = 'Ground_truth', cols = c('blue','red')) + ggtitle("Ground truth potency") + xlab("UMAP_1") + ylab("UMAP_2")
      
    ###single-cell spearman correlation representation
figc.3 = ggplot(adata@meta.data, aes(x = Ground_truth, y = stemFinder)) + geom_jitter(aes(color = Ground_truth)) + stat_smooth(method = 'lm', formula = y ~ x, geom = 'smooth') + xlab("Ground truth potency") + ylab("Inverted stemFinder score") + ggtitle("Single-cell Spearman Correlation") + stat_cor(method = 'spearman', label.x = 0.6, label.y = 1)
    
    ###phenotypic Spearman correlation representation 
clusters = as.character(unique(adata$Phenotype))
meanpotency = data.frame("cluster" = clusters, "sf_potency" = rep(NA, length(clusters)), "ground_truth" = rep(NA, length(clusters)))
for (i in clusters){
  meanpotency$sf_potency[which(clusters == i)] = mean(adata$stemFinder[adata$Phenotype == i])
  meanpotency$ground_truth[which(clusters == i)] = mean(adata$Ground_truth[adata$Phenotype == i])
}
figc.4 = ggplot(meanpotency, aes(x = ground_truth, y = sf_potency)) + geom_point(aes(color = cluster)) + theme_classic() + geom_label_repel(aes(label = cluster, color = cluster), box.padding = 0.35, point.padding = 0.3, label.padding = 0.1) + ggtitle("Phenotypic Spearman correlation") + ylab("Mean inverted stemFinder score") + xlab("Ground truth potency") + stat_smooth(method = 'lm', formula = y ~ x, geom = 'smooth', alpha = 0) + stat_cor(method = 'spearman', label.x = 0.6, label.y=1)

    ###pct.recov representation
thresh_cells = rownames(adata@meta.data[adata$stemFinder < quantile(adata$stemFinder, (1/length(unique(adata$Ground_truth)))),])
figc.5 = DimPlot(adata, cells.highlight = thresh_cells, cols.highlight = 'blue') + ggtitle("Recovered potent cells") + xlab("UMAP_1") + ylab("UMAP_2")
  
########################################

# Fig 2 - murine dentate gyrus example

adata = readRDS("DentateGyrus_mouse_phenotypes_10X.rds")

#2A - B: UMAP of cell phenotype (A) and CC phase (B)
fig2a <- DimPlot(adata, group.by = 'Phenotype', cols = c('pink','red','yellow','orange','green','darkgreen','black','grey','blue','darkblue','brown','darkred','purple','magenta','tan')) + ggtitle("Original cell phenotypes")
fig2b <- DimPlot(adata, group.by = 'Phase') + ggtitle("Cell cycle phase")

#2C - D: UMAP of potency values (GT and sF)
fig2c <- FeaturePlot(adata, features = 'Ground_truth', cols = c('blue','red')) + ggtitle("Ground truth potency")
fig2d <- FeaturePlot(adata, features = 'stemFinder', cols = c('blue','red')) + ggtitle("stemFinder score")

#2E: UMAP with highlighted non-selected RGL cells
highlypot.sF <-rownames(adata@meta.data[adata$stemFinder < quantile(adata$stemFinder, (1/length(unique(adata$Ground_truth)))),]) #use same threshold in pct.recov() to get highly potent cells per potency calculator
highlypot.cyto <-rownames(adata@meta.data[adata$CytoTRACE_invert < quantile(adata$CytoTRACE_invert, (1/length(unique(adata$Ground_truth)))),]) 
highlypot.ccat <-rownames(adata@meta.data[adata$ccat_invert < quantile(adata$ccat_invert, (1/length(unique(adata$Ground_truth)))),]) 
truepot <- rownames(adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),]) #cells with min ground truth potency
computedpot <- union(union(highlypot.sF, highlypot.cyto), highlypot.ccat) #cells marked as highly potent by any calculator
length(truepot[!truepot %in% computedpot]) #num highly potent cells missed by all calculators
fig2e <- DimPlot(adata, cells.highlight = truepot[!truepot %in% computedpot], pt.size = 0.5, sizes.highlight = 0.5) + ggtitle("Non-selected RGL cells") + NoLegend()

#2F: Violin plots of potency scores by phenotype
fig2f_sF <- VlnPlot(adata, features = 'stemFinder', group.by = 'Phenotype') + ggtitle("stemFinder")
fig2f_cyto <- VlnPlot(adata, features = 'CytoTRACE_invert', group.by = 'Phenotype') + ggtitle("CytoTRACE")
fig2f_ccat <- VlnPlot(adata, features = 'ccat_invert', group.by = 'Phenotype') + ggtitle("CCAT")

#2G - I: RGL cells: Leiden clusters, potency bin, and CC phase
rgl <- readRDS("dentategyrus_rglonly.rds") #subsetted RGL cells, re-analyzed using BasicSeuratPipeline.R, clustered with Leiden resolution 0.25
fig2g <- DimPlot(rgl, group.by = 'seurat_clusters') + ggtitle("RGL cells, Leiden clusters")  
rgl$highly.pot <- 'Not marked as highly potent'
rgl@meta.data[colnames(rgl)[colnames(rgl) %in% computedpot],]$highly.pot = 'Marked as highly potent'
fig2h <- DimPlot(rgl, group.by = 'highly.pot') + ggtitle("RGL cells, potency designation")
fig2i <- DimPlot(rgl, group.by = 'Phase') + ggtitle("RGL cells, cell cycle phase")

#2J: RGL cells: number of counts and features per cluster
fig2j_counts <- VlnPlot(rgl, features = 'nCount_RNA', group.by = 'seurat_clusters') + ggtitle("Counts") + xlab("RGL subcluster")
fig2j_feat <- VlnPlot(rgl, features = 'nFeature_RNA', group.by = 'seurat_clusters') + ggtitle("Features") + xlab("RGL subcluster")

#2K: heat map of DE genes
Idents(rgl) <- 'seurat_clusters'
de.genes <- FindMarkers(rgl, ident.1 = '0', only.pos = F, min.pct = 0, logfc.threshold = 0)
de.genes$p_val_adj_fdr <- p.adjust(p = de.genes$p_val, method = 'fdr') #FDR adjustment 
de.genes[c('S100a16','Gfap','Sparcl1','Pla2g7','Mycn','Hes6'),] #genes in heat map
fig6k <- DoHeatmap(adata, features = c('S100a16','Gfap','Sparcl1','Pla2g7','Mycn','Hes6'), group.by = 'Phenotype', group.colors = c('pink','red','yellow','orange','green','darkgreen','black','grey','blue','darkblue','brown','darkred','purple','magenta','tan'), angle = 90, size = 3)

#2L: revised phenotypes
rev.pheno <- read.csv("../../ValidationandGeneLists/DentateGyrusPhenotype_RevisedAnnotandGT.csv", header = T, row.names = 1)
adata$Phenotype_revised = rev.pheno$Revised_Phenotype
fig2l <- DimPlot(adata, group.by = 'Phenotype_revised', pt.size = 0.5, cols = c('pink','red','black','yellow','orange','green','darkgreen','black','grey','blue','darkblue','brown','darkred','purple','magenta','tan')) + ggtitle("Revised cell phenotypes")

#2M: GSEA, RGL clust 0 vs. 1
m_df = msigdbr(species = "Mus musculus", category = 'C8') %>%  #C8 is cell type signatures
  dplyr::select(gs_name, gene_symbol)
genelist_rf = de.genes$avg_log2FC
names(genelist_rf) = rownames(de.genes)
genelist_rf = na.omit(genelist_rf)
genelist_rf = sort(genelist_rf, decreasing=TRUE)
gsea_res = GSEA(genelist_rf, TERM2GENE = m_df, pvalueCutoff=0.05, verbose = TRUE)
#dot plot
dot_df = as.data.frame(gsea_res[gsea_res$p.adjust < 0.05,])
dot_df$absNES = abs(dot_df$NES)
dot_df = arrange(dot_df, desc(absNES),p.adjust)  #reorder by NES then by p value
dot_df$type = "upregulated"
dot_df$type[dot_df$NES < 0] = "downregulated"
dot_df$Description = stri_replace_all_fixed(dot_df$Description, "_", " ") 
#plotting terms 1, 2, 6, 11, 16, 35, 23, 43, 45: astrocyte & progenitor terms
fig6m <- ggplot(dot_df[c(1,2,6,11,16,23,35,43,45),], aes(x = NES, y = fct_reorder(Description, NES))) + 
  geom_point(aes(size=2, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) +
  ggtitle("GSEA, RGL cluster 0 vs. 1") +
  scale_y_discrete(labels = function(y) str_wrap(y, width=45)) +
  theme(axis.text.y = element_text(size=11)) +
  guides(size="none")

##########################################

#Fig 3 - Performance comparison for UMI data 

  ## four performance metrics box plots
load("df_performance.rda")
fig3a = ggplot(df_all[df_all$iter == 1 & df_all$genelist == 'Full_Regev' & df_all$method %in% c('stemFinder','CytoTRACE','CCAT'),], aes(x = method, y = spear_all)) + geom_boxplot(aes(group = method)) + geom_point(aes(fill = dataset_forplotting), size=2, shape=21,position = position_dodge(0.1)) + theme_clean() + ggtitle("Single-Cell Spearman Correlation") + ylab("Single-cell Spearman Correlation") + xlab("Potency quantification method") + scale_fill_discrete(name = 'Validation Datasets') + scale_colour_discrete(guide = 'none') + theme(axis.text.x=element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), plot.title = element_text(size=20)) 
fig3b = ggplot(df_all[df_all$iter == 1 & df_all$genelist == 'Full_Regev' & df_all$method %in% c('stemFinder','CytoTRACE','CCAT'),], aes(x = method, y = spear_pheno)) + geom_boxplot(aes(group = method)) + geom_point(aes(fill = dataset_forplotting), size=2, shape=21,position = position_dodge(0.1)) + theme_clean() + ggtitle("Phenotypic Spearman Correlation") + ylab("Phenotypic Spearman Correlation") + xlab("Potency quantification method") + scale_fill_discrete(name = 'Validation Datasets') + scale_colour_discrete(guide = 'none') + theme(axis.text.x=element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), plot.title = element_text(size=20)) 
fig3c = ggplot(df_all[df_all$iter == 1 & df_all$genelist == 'Full_Regev' & df_all$method %in% c('stemFinder','CytoTRACE','CCAT'),], aes(x = method, y = auc)) + geom_boxplot(aes(group = method)) + geom_point(aes(fill = dataset_forplotting), size=2, shape=21,position = position_dodge(0.1)) + theme_clean() + ggtitle("Discrimination Accuracy") + ylab("AUC") + xlab("Potency quantification method") + scale_fill_discrete(name = 'Validation Datasets') + scale_colour_discrete(guide = 'none') + theme(axis.text.x=element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), plot.title = element_text(size=20)) 
fig3d = ggplot(df_all[df_all$iter == 1 & df_all$genelist == 'Full_Regev' & df_all$method %in% c('stemFinder','CytoTRACE','CCAT'),], aes(x = method, y = pct.recov)) + geom_boxplot(aes(group = method)) + geom_point(aes(fill = dataset_forplotting), size=2, shape=21,position = position_dodge(0.1)) + theme_clean() + ggtitle("Percent recovery of highly potent cells") + ylab("Percent recovery") + xlab("Potency quantification method") + scale_fill_discrete(name = 'Validation Datasets') + scale_colour_discrete(guide = 'none') + theme(axis.text.x=element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), plot.title = element_text(size=20)) 

  ## memory & run time
load("runandmem_ubuntus3_periphblood.rda")
fig3e = ggplot(runandmem_copy, aes(x = NumCells, y = MemoryKB)) + geom_line(aes(color = Method), size=2, alpha = 0.5) + geom_point(aes(color = Method), size=3, alpha = 0.5) + theme_clean() + ggtitle("Memory usage") + ylab("Memory (kB)") + xlab("Number of cells") + theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
fig3f = ggplot(runandmem_copy, aes(x = NumCells, y = RunTime.sec.)) + geom_line(aes(color = Method), size=2, alpha = 0.5) + geom_point(aes(color = Method), size=3, alpha = 0.5) + theme_clean() + ggtitle("Run time") + ylab("Run time (sec)") + xlab("Number of cells") + theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))

########################################

#Fig 4 - Cell cycle gene expression in Tabula Muris BMMC 10X dataset

adata = readRDS("MurineBoneMarrow10X_GSE109774.rds")
Idents(adata) = 'Phenotype'

  #4A: GSEA for stem cells vs. all other cells
markers_stem = FindMarkers(adata, ident.1 = "Stem_Progenitors", min.pct = 0, logfc.threshold=0, only.pos=FALSE)
m_df = msigdbr(species = "Mus musculus", category = 'C8') %>% 
  dplyr::select(gs_name, gene_symbol)
genelist_rf = markers_stem$avg_log2FC
names(genelist_rf) = rownames(markers_stem)
genelist_rf = na.omit(genelist_rf)
genelist_rf = sort(genelist_rf, decreasing=TRUE)
gsea_res = GSEA(genelist_rf, TERM2GENE = m_df, pvalueCutoff=0.05, verbose = TRUE)
    #dot plot
dot_df = as.data.frame(gsea_res[gsea_res$p.adjust < 0.05,])
dot_df$absNES = abs(dot_df$NES)
dot_df = arrange(dot_df, desc(absNES),p.adjust)  #reorder by NES then by p value
dot_df$type = "upregulated"
dot_df$type[dot_df$NES < 0] = "downregulated"
dot_df$Description = stri_replace_all_fixed(dot_df$Description, "_", " ") #adds space
    #terms 2, 4, 6, 7, 10, 21, 43, 46, 59 with NES>0 are progenitor cell types
fig4a <- ggplot(dot_df[dot_df$type == 'upregulated',][c(2,4,6,7,10,21,43,46,59),], aes(x = NES, y = fct_reorder(Description, NES))) + 
  geom_point(aes(size=2, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) +
  ggtitle("MSigDB GSEA results") +
  scale_y_discrete(labels = function(y) str_wrap(y, width=45)) +
  theme(axis.text.y = element_text(size=11)) +
  guides(size="none") 

  #4B-D heat maps of binarized scaled CC gene expression
load("G1_MsigDB_mouse.rda") 
markers_cc = g1_mouse #set as needed for figB-D

expDat_sub = as.matrix(adata@assays$RNA@scale.data)[markers_cc[markers_cc %in% rownames(adata)],]
sf = adata@meta.data$stemFinder
names(sf) = rownames(adata@meta.data)
sf = sf[order(sf, decreasing = F)]
toplot = expDat_sub[,order(match(colnames(expDat_sub), names(sf)))]
toplot_binary = (as.matrix(toplot > 0))*1 #binarized scaled expression matrix with cells in order of ascending inverted stemFinder score
rows.cor <- cor(t(expDat_sub), use = "pairwise.complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor)) #dendrogram

bar = as.numeric(as.factor(names(sf))) #colors
newcol <- colorRampPalette(c('blue','orange'))
ncols <- length(sf)
colSide <- newcol(ncols)
dev.off()
heatmap(toplot_binary, Colv=NA, Rowv = as.dendrogram(hclust.row), scale = 'none', col = c('black','red'), ColSideColors = colSide, cexRow = 0.5) # distfun=function (y) dist(y,method = "euclidean"))
      #note: to save as TIFF, just run heatmap code after tiff() line, not object 
      #cexRow is rowname font size

  #alternative heat map: with phenotype annotations
# meta = adata@meta.data
# meta$Phenotype_simplified = meta$Phenotype
# col_annots <- meta[,c('Phenotype','Ground_truth','stemFinder')]
# col_annots = with(col_annots, col_annots[order(stemFinder, decreasing = F),]) 
# toplot2 = toplot_binary[,rownames(col_annots)] #can use toplot if not binarized is desired
# pheatmap(toplot2, scale = 'none', color = colorRampPalette(c('black','red'))(100), cluster_rows = T, show_colnames=F, cluster_cols = F, annotation_col = col_annots)

  #4E: gene set score vs. sF score 
adata = gene_set_score(gene.set = s_genes_mouse, adata) #returns in meta column called gene.set.score
fig4e_s <- ggplot(adata@meta.data, aes(x = stemFinder, y = gene.set.score)) + geom_point(aes(color = stemFinder)) + scale_color_gradient(low = 'blue', high = 'orange', breaks = c(0, max(adata$stemFinder)), labels = c('Min','Max'), name = 'stemFinder score') + geom_smooth(method = 'loess', colour = 'blue') + ylab("S phase gene set score") + xlab("stemFinder score") 
    ## etc. for G2M and G1 gene sets

  #4F-G: histograms of CC phase
toplot = adata@meta.data[,c('Phase','stemFinder')]
fig4f <- ggplot(toplot, aes(stemFinder)) + geom_histogram(binwidth = 0.05, aes(fill = Phase)) + ggtitle("Distribution of CC Phase over stemFinder Trajectory") + theme_classic() + ylab("Cell count") + xlab("stemFinder score")
toplot = adata@meta.data[,c('Phase','Ground_truth')]
fig4g <- ggplot(toplot, aes(Ground_truth)) + geom_histogram(binwidth = .5, aes(fill = Phase)) + ggtitle("Distribution of CC Phase over Ground truth potency") + theme_classic() + ylab("Cell count") + xlab("Ground truth potency")

#######################################################

#Fig 5 - Lineage tracing
adata = readRDS("cospar_hematopoeitic_lineagetracing.rds")
load("clone.rda") #matrix of cell barcodes x clone IDs
load("times.rda") #time points for each cell barcode: 2, 4, or 6 days

    #5A-B: UMAP cell type, time point, sF score
fig5a <- DimPlot(adata, group.by = 'state_info') + ggtitle("Cell type")
fig5b <- DimPlot(adata, group.by = 'time_info') + ggtitle("Time point")
fig5c <- FeaturePlot(adata, features = 'stemFinder') + ggtitle("stemFinder score")

    #5D: box plot, sF vs number lineages

#Limit analysis to only clones that contain cells from time point 6day 
clones_incl = c()
for (c in colnames(clone_bin)){
  times_c = names(table(adata@meta.data[names(which(clone_bin[,c] == TRUE)),]$time_info))
  if("6" %in% times_c){
    clones_incl = c(clones_incl, c)
  }
}

#Count how many different lineages a given clone gives rise to
lin_counts = data.frame("clone_ID" = colnames(clone_bin), "Num_lineages" = rep("", 500), "Mean_sF_undiff" = rep("", 500), "Min_sF_undiff" = rep("", 500), "Mean_CCATinvert_undiff" = rep("",500), "Min_CCATinvert_undiff" = rep("", 500), "Mean_Cytoinvert_undiff" = rep("",500), "Min_Cytoinvert_undiff" = rep("", 500))
for (c in colnames(clone_bin)){
  #count number of lineages
  lins = unique(names(table(adata@meta.data[names(which(clone_bin[,c] == TRUE)),]$state_info)))
  if("undiff" %in% times_c){
    lins = length(lins) - 1
  }else {lins = length(lins)}
  lin_counts[lin_counts$clone_ID == c,]$Num_lineages = lins
  
  #retrieve minimum and mean stemFinder scores for time point 2 undifferentiated cells in that lineage
  bars = names(which(clone_bin[,c] == TRUE))
  min_time = min(adata@meta.data[bars,]$time_info)
  avg_sF = mean(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$stemFinder)
  min_sF = min(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$stemFinder)
  lin_counts[lin_counts$clone_ID == c,]$Min_sF_undiff = min_sF
  lin_counts[lin_counts$clone_ID == c,]$Mean_sF_undiff = avg_sF
  
  avg_CCAT = mean(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$ccat_invert)
  min_CCAT = min(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$ccat_invert)
  lin_counts[lin_counts$clone_ID == c,]$Min_CCATinvert_undiff = min_CCAT
  lin_counts[lin_counts$clone_ID == c,]$Mean_CCATinvert_undiff = avg_CCAT
  
  avg_Cyto = mean(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$CytoTRACE_invert)
  min_Cyto = min(adata@meta.data[rownames(adata@meta.data) %in% bars & adata$time_info == min_time,]$CytoTRACE_invert)
  lin_counts[lin_counts$clone_ID == c,]$Min_Cytoinvert_undiff = min_Cyto
  lin_counts[lin_counts$clone_ID == c,]$Mean_Cytoinvert_undiff = avg_Cyto
}

fig5d_sF <- ggplot(lin_counts[,c(2,4)], aes(y = as.numeric(Min_sF_undiff), x = as.numeric(Num_lineages))) + geom_boxplot(aes(group = Num_lineages)) + guides(fill = 'none', color = 'none') + geom_point(aes(color = as.numeric(Min_sF_undiff))) + xlab("Number of downstream lineages") + ylab("Minimum stemFinder score") + ggtitle("stemFinder vs. number of lineages")
fig5d_cyto <- ggplot(lin_counts[,c(2,8)], aes(y = as.numeric(Min_Cytoinvert_undiff), x = as.numeric(Num_lineages))) + geom_boxplot(aes(group = Num_lineages)) + guides(fill = 'none', color = 'none') + geom_point(aes(color = as.numeric(Min_Cytoinvert_undiff))) + xlab("Number of downstream lineages") + ylab("Minimum inverted CytoTRACE score") + ggtitle("CytoTRACE vs. number of lineages")
fig5d_ccat <- ggplot(lin_counts[,c(2,6)], aes(y = as.numeric(Min_CCATinvert_undiff), x = as.numeric(Num_lineages))) + geom_boxplot(aes(group = Num_lineages)) + guides(fill = 'none', color = 'none') + geom_point(aes(color = as.numeric(Min_CCATinvert_undiff))) + xlab("Number of downstream lineages") + ylab("Minimum inverted CCAT score") + ggtitle("CCAT vs. number of lineages")

    #5E-F: UMAP marker gene expression
fig5e <- FeaturePlot(adata, features = 'Ngp') 
fig5f <- FeaturePlot(adata, features = 'Mmp8') 

##################################################

# SUPPLEMENTAL FIGURES

#Fig S1

  #S1A: box plot, single-cell Spearman for CC gene expression heterogeneity vs. mean CC expression 
load("df_performance.rda")
s1a <- ggplot(df_all[df_all$method %in% c('stemFinder','GeneSetScore') & df_all$genelist %in% c('Full_Regev','Full_Regev_Plus_G1','G1','G2M','S'),], aes(x = method, y = spear_all)) + geom_point(position = position_jitterdodge(jitter.width=0.1, dodge.width = 0.5), aes(group = genelist, color = genelist)) + geom_boxplot(aes(group = method), alpha = 0) + theme_clean() + coord_flip() + labs(color = 'Gene list') + xlab("Method") + ylab("Single-Cell Spearman Correlation") + ggtitle("Correlation of gene set score vs. expression stochasticity with ground truth potency") + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14)) + scale_color_discrete(labels = c('S + G2M','S + G2M + G1','G1 only','G2M only','S only'))
stat.test <- df_all[df_all$method %in% c('stemFinder','GeneSetScore') & df_all$genelist %in% c('Full_Regev','Full_Regev_Plus_G1','G1','G2M','S'),] %>% pairwise_t_test(spear_all ~ method, paired = T, p.adjust.method = 'bonferroni')

  #S1B: box plot, performance of stemFinder vs. gene set score 
s1b_ss <- ggplot(df_all[df_all$method %in% c('stemFinder','GeneSetScore') & df_all$genelist == 'Full_Regev',], aes(x = spear_all, y = method)) + geom_boxplot(aes(group = method, fill = method)) + geom_point(shape=1) + ggtitle("Single-cell Spearman correlation with ground truth potency") + xlab("Single-cell Spearman Correlation") + ylab("Method") + scale_fill_discrete(name = 'Method', label = c('Gene Set Score','stemFinder')) + theme_clean() + theme(plot.title=element_text(size=12))
s1b_sp <- ggplot(df_all[df_all$method %in% c('stemFinder','GeneSetScore') & df_all$genelist == 'Full_Regev',], aes(x = spear_pheno, y = method)) + geom_boxplot(aes(group = method, fill = method)) + geom_point(shape=1) + ggtitle("Phenotypic Spearman correlation with ground truth potency") + xlab("Phenotypic Spearman Correlation") + ylab("Method") + scale_fill_discrete(name = 'Method', label = c('Gene Set Score','stemFinder')) + theme_clean() + theme(plot.title=element_text(size=12))
s1b_auc <- ggplot(df_all[df_all$method %in% c('stemFinder','GeneSetScore') & df_all$genelist == 'Full_Regev',], aes(x = auc, y = method)) + geom_boxplot(aes(group = method, fill = method)) + geom_point(shape=1) + ggtitle("Discrimination Accuracy") + xlab("AUC") + ylab("Method") + scale_fill_discrete(name = 'Method', label = c('Gene Set Score','stemFinder')) + theme_clean() + theme(plot.title=element_text(size=12))
s1b_pctrec <- ggplot(df_all[df_all$method %in% c('stemFinder','GeneSetScore') & df_all$genelist == 'Full_Regev',], aes(x = pct.recov, y = method)) + geom_boxplot(aes(group = method, fill = method)) + geom_point(shape=1) + ggtitle("Percent recovery of highly potent cells") + xlab("Percentage Recovery") + ylab("Method") + scale_fill_discrete(name = 'Method', label = c('Gene Set Score','stemFinder')) + theme_clean() + theme(plot.title=element_text(size=12))

#####################################################

  #Supplemental Fig 2: robustness heat maps
      ## order of all figures: single-cell Spearman, phenotypic Spearman, AUC, pct.recov
      ## .rda files generated using automated_performance.R and compute_performance.R files, reformatted using robustness_format.R

#S2A-D: K robustness
load("df_robustness_k.rda")
s2a <- ggplot(df_robust_k, aes(x = kratio, y = dataset_forplot, fill = neg_deviation_spearss)) + geom_tile() + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-2, 2, 0.5), limits=c(-2, 2), name = 'Deviation in Performance') + ggtitle("Robustness to changes in K - Single-cell Spearman Correlation") + ylab("Validation dataset") + xlab("K / Kideal")
s2b <- ggplot(df_robust_k, aes(x = kratio, y = dataset_forplot, fill = neg_deviation_spearp)) + geom_tile() + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-2, 2, 0.5), limits=c(-2, 2), name = 'Deviation in Performance') + ggtitle("Robustness to changes in K - Phenotypic Spearman Correlation") + ylab("Validation dataset") + xlab("K / Kideal")
s2c <- ggplot(df_robust_k, aes(x = kratio, y = dataset_forplot, fill = neg_deviation_auc)) + geom_tile() + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-1, 1, 0.5), limits=c(-1, 1), name = 'Deviation in Performance') + ggtitle("Robustness to changes in K - Discrimination Accuracy") + ylab("Validation dataset") + xlab("K / Kideal")
s2d <- ggplot(df_robust_k, aes(x = kratio, y = dataset_forplot, fill = neg_deviation_pctrecov)) + geom_tile() + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-100, 100, 50), limits=c(-100, 100), name = 'Deviation in Performance') + ggtitle("Robustness to changes in K - Percent Recovery") + ylab("Validation dataset") + xlab("K / Kideal")

#S2E-H: Cell cycle gene list (GO, Kegg, Regev)
load("df_robustness_genelist.rda")
s2e <- ggplot(df_robust_genelist[df_robust_genelist$genelist %in% c('Full_Regev','KEGG','GeneOntology'),], aes(x = genelist, y = dataset_forplot, fill = neg_deviation_spearss)) + geom_tile() + ggtitle("Robustness to changes in cell cycle gene list - Single-cell Spearman Correlation") + ylab("Validation dataset") + xlab("Gene List") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-2, 2, 0.5), limits=c(-2, 2), name = 'Deviation in Performance') + theme_bw()
s2f <- ggplot(df_robust_genelist[df_robust_genelist$genelist %in% c('Full_Regev','KEGG','GeneOntology'),], aes(x = genelist, y = dataset_forplot, fill = neg_deviation_spearp)) + geom_tile() + ggtitle("Robustness to changes in cell cycle gene list - Phenotypic Spearman Correlation") + ylab("Validation dataset") + xlab("Gene List") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-2, 2, 0.5), limits=c(-2, 2), name = 'Deviation in Performance') + theme_bw()
s2g <- ggplot(df_robust_genelist[df_robust_genelist$genelist %in% c('Full_Regev','KEGG','GeneOntology'),], aes(x = genelist, y = dataset_forplot, fill = neg_deviation_auc)) + geom_tile() + ggtitle("Robustness to changes in cell cycle gene list - Discrimination Accuracy") + ylab("Validation dataset") + xlab("Gene List") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-1, 1, 0.5), limits=c(-1, 1), name = 'Deviation in Performance') + theme_bw()
s2h <- ggplot(df_robust_genelist[df_robust_genelist$genelist %in% c('Full_Regev','KEGG','GeneOntology'),], aes(x = genelist, y = dataset_forplot, fill = neg_deviation_pctrecov)) + geom_tile() + ggtitle("Robustness to changes in cell cycle gene list - Percent Recovery") + ylab("Validation dataset") + xlab("Gene List") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-100, 100, 50), limits=c(-100, 100), name = 'Deviation in Performance') + theme_bw()

#S2I-L: Proportional downsampling of all cells
load("df_robustness_downsample_allandpotent.rda")
df_robust_ds_all <- df_robust_ds[df_robust_ds$ds_type == 'All cell types',]
s2i <- ggplot(df_robust_ds_all, aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_spearss)) + geom_tile() + ggtitle("Robustness to proportional downsampling - Single-cell Spearman Correlation") + ylab("Validation dataset") + xlab("Downsampling ratio") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-2, 2, .5), limits=c(-2, 2), name = 'Deviation in Performance') + theme_bw()
s2j <- ggplot(df_robust_ds_all, aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_spearp)) + geom_tile() + ggtitle("Robustness to proportional downsampling - Phenotypic Spearman Correlation") + ylab("Validation dataset") + xlab("Downsampling ratio") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-2, 2, .5), limits=c(-2, 2), name = 'Deviation in Performance') + theme_bw()
s2k <- ggplot(df_robust_ds_all, aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_auc)) + geom_tile() + ggtitle("Robustness to proportional downsampling - Discrimination Accuracy") + ylab("Validation dataset") + xlab("Downsampling ratio") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-1, 1, .5), limits=c(-1, 1), name = 'Deviation in Performance') + theme_bw()
s2l <- ggplot(df_robust_ds_all, aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_pctrecov)) + geom_tile() + ggtitle("Robustness to proportional downsampling - Percent Recovery") + ylab("Validation dataset") + xlab("Downsampling ratio") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-100, 100, 50), limits=c(-100, 100), name = 'Deviation in Performance') + theme_bw()

#S2M-P: Downsampling of most potent population
df_robust_ds_pot <- df_robust_ds[df_robust_ds$ds_type == 'Highly potent only',]
s2m <- ggplot(df_robust_ds_pot, aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_spearss)) + geom_tile() + ggtitle("Robustness to downsampling of highly potent cells - Single-cell Spearman Correlation") + ylab("Validation dataset") + xlab("Downsampling ratio") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-2, 2, .5), limits=c(-2, 2), name = 'Deviation in Performance') + theme_bw()
s2n <- ggplot(df_robust_ds_pot, aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_spearp)) + geom_tile() + ggtitle("Robustness to downsampling of highly potent cells - Phenotypic Spearman Correlation") + ylab("Validation dataset") + xlab("Downsampling ratio") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-2, 2, .5), limits=c(-2, 2), name = 'Deviation in Performance') + theme_bw()
s2o <- ggplot(df_robust_ds_pot, aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_auc)) + geom_tile() + ggtitle("Robustness to downsampling of highly potent cells - Discrimination Accuracy") + ylab("Validation dataset") + xlab("Downsampling ratio") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-1, 1, .5), limits=c(-1, 1), name = 'Deviation in Performance') + theme_bw()
s2p <- ggplot(df_robust_ds_pot, aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_pctrecov)) + geom_tile() + ggtitle("Robustness to downsampling of highly potent cells - Percent Recovery") + ylab("Validation dataset") + xlab("Downsampling ratio") + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0, breaks = seq(-100, 100, 50), limits=c(-100, 100), name = 'Deviation in Performance') + theme_bw()

#S2Q-T: Regev vs. random gene list (box plot)
s2q <- ggplot(df_robust_genelist[df_robust_genelist$genelist %in% c('Full_Regev','Random'),], aes(x = genelist, y = spear_all_sF)) + geom_point(aes(color = dataset_forplot)) + geom_boxplot(aes(group = genelist), alpha = 0) + ggtitle("Random vs. cell cycle genes") + xlab("Validation dataset") + ylab("Single-cell Spearman Correlation") + theme_bw() + scale_color_discrete(name = 'Validation Dataset') + scale_x_discrete(labels = c('S + G2M genes','Random'))
s2r <- ggplot(df_robust_genelist[df_robust_genelist$genelist %in% c('Full_Regev','Random'),], aes(x = genelist, y = spear_pheno_sF)) + geom_point(aes(color = dataset_forplot)) + geom_boxplot(aes(group = genelist), alpha = 0) + ggtitle("Random vs. cell cycle genes") + xlab("Validation dataset") + ylab("Phenotypic Spearman Correlation") + theme_bw() + scale_color_discrete(name = 'Validation Dataset') + scale_x_discrete(labels = c('S + G2M genes','Random'))
s2s <- ggplot(df_robust_genelist[df_robust_genelist$genelist %in% c('Full_Regev','Random'),], aes(x = genelist, y = auc_sF)) + geom_point(aes(color = dataset_forplot)) + geom_boxplot(aes(group = genelist), alpha = 0) + ggtitle("Random vs. cell cycle genes") + xlab("Validation dataset") + ylab("Discrimination Accuracy") + theme_bw() + scale_color_discrete(name = 'Validation Dataset') + scale_x_discrete(labels = c('S + G2M genes','Random'))
s2t <- ggplot(df_robust_genelist[df_robust_genelist$genelist %in% c('Full_Regev','Random'),], aes(x = genelist, y = pct.recov_sF)) + geom_point(aes(color = dataset_forplot)) + geom_boxplot(aes(group = genelist), alpha = 0) + ggtitle("Random vs. cell cycle genes") + xlab("Validation dataset") + ylab("Percentage Recovery") + theme_bw() + scale_color_discrete(name = 'Validation Dataset') + scale_x_discrete(labels = c('S + G2M genes','Random'))

# example of significance tests: one way anova + tukey HSD
res.aov <- aov(pct.recov_sF ~ ds_ratio, data = df_robust_ds_pot)
summary(res.aov) #Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
TukeyHSD(res.aov) #do this if p adj < 0.05 for anova 

#####################################

#Supplementary figure 3: oligodendrocyte validation dataset

adata = readRDS("Oligodendrocytephenotypes_C1mouse.rds")

#S3A: UMAP cell phenotype
s3a <- DimPlot(adata, group.by = 'Phenotype', cols = c('blue','brown','orange','navy','green','black','darkgreen','purple','magenta','tan','pink','red')) + ggtitle("Cell phenotype")

#S3B-C: violin plot of sF, CytoTRACE, and CCAT scores
s3b <- VlnPlot(adata, features = 'stemFinder', group.by = 'Phenotype', cols = c('blue','brown','orange','navy','green','black','darkgreen','purple','magenta','tan','pink','red')) + ggtitle("Inverted stemFinder score") + xlab("Cell phenotype")
s3c_cyto <- VlnPlot(adata, features = 'CytoTRACE_invert', group.by = 'Phenotype', cols = c('blue','brown','orange','navy','green','black','darkgreen','purple','magenta','tan','pink','red')) + ggtitle("Inverted CytoTRACE score")+ xlab("Cell phenotype")
s3c_ccat <- VlnPlot(adata, features = 'ccat_invert', group.by = 'Phenotype', cols = c('blue','brown','orange','navy','green','black','darkgreen','purple','magenta','tan','pink','red')) + ggtitle("Inverted CCAT score")+ xlab("Cell phenotype")

#S3D: feature plot of potency scores
s3d <- FeaturePlot(adata, features = 'stemFinder', cols = c('blue','red'))
