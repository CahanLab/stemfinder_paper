### Code for stemFinder manuscript figures

# updated 4/24/23 K.Noller

#Setup
library(ggplot2)
library(ggthemes)
library(Seurat)
library(RColorBrewer)
setwd("~/Dropbox (CahanLab)/Kathleen.Noller/Stochasticity/Paper&Presentations/ManuscriptFigures_300dpiTIFF/")

#Example code for saving as TIFF 300dpi figures
tiff("test.tiff", units="in", width=6, height=5, res=300)
#insert ggplot code
dev.off()

########################################

# Fig 1 - schematic: majority generated with BioRender

    ##1C: Feature plots of BMMC 10X Tabula Muris, scatter plot for Spearman and pct.recov
adata = readRDS("../../ValidationandGeneLists/bmmc_mouse_fromcytotrace_10X_ginicytoandccat.rds") + ggtitle("Inverted stemFinder score") + xlab("UMAP_1") + ylab("UMAP_2")
figc.1 = FeaturePlot(adata, features = 'gini_cellcycle_normandinvert', cols = c('blue','red')) + ggtitle("Ground truth potency") + xlab("UMAP_1") + ylab("UMAP_2")
figc.2 = FeaturePlot(adata, features = 'Ground_truth', cols = c('blue','red'))
      
    ###single-cell spearman correlation representation
figc.3 = ggplot(adata@meta.data, aes(x = Ground_truth, y = gini_cellcycle_normandinvert)) + geom_jitter(aes(color = Ground_truth)) + stat_smooth(method = 'lm', formula = y ~ x, geom = 'smooth') + xlab("Ground truth potency") + ylab("Inverted stemFinder score") + ggtitle("Single-cell Spearman Correlation") + stat_cor(method = 'spearman', label.x = 0.6, label.y = 1)
    
    ###phenotypic Spearman correlation representation 
clusters = as.character(unique(adata$Phenotype))
meanpotency = data.frame("cluster" = clusters, "sf_potency" = rep(NA, length(clusters)), "ground_truth" = rep(NA, length(clusters)))
for (i in clusters){
  meanpotency$sf_potency[which(clusters == i)] = mean(adata$gini_cellcycle_normandinvert[adata$Phenotype == i])
  meanpotency$ground_truth[which(clusters == i)] = mean(adata$Ground_truth[adata$Phenotype == i])
}
figc.4 = ggplot(meanpotency, aes(x = ground_truth, y = sf_potency)) + geom_point(aes(color = cluster)) + theme_classic() + geom_label_repel(aes(label = cluster, color = cluster), box.padding = 0.35, point.padding = 0.3, label.padding = 0.1) + ggtitle("Phenotypic Spearman correlation") + ylab("Mean inverted stemFinder score") + xlab("Ground truth potency") + stat_smooth(method = 'lm', formula = y ~ x, geom = 'smooth', alpha = 0) + stat_cor(method = 'spearman', label.x = 0.6, label.y=1)

    ###pct.recov representation
thresh_cells = rownames(adata@meta.data[adata$gini_cellcycle_normandinvert < quantile(adata$gini_cellcycle_normandinvert, (1/length(unique(adata$Ground_truth)))),])
figc.5 = DimPlot(adata, cells.highlight = thresh_cells) + ggtitle("Recovered potent cells") + xlab("UMAP_1") + ylab("UMAP_2")
  
########################################

#Fig 2 - C. elegans hypodermis 

adata = readRDS("../../ValidationandGeneLists/Celegans_Hypodermis_10X_giniandcytoandccat.rds")
fig2a.1 = ggplot(adata@meta.data, aes(x = Ground_truth, y = gini_cellcycle_normandinvert)) + geom_point(color = 'darkblue', size = 0.5) + ggtitle("stemFinder") + xlab("Ground truth potency") + ylab("Inverted potency score") + theme_classic()
fig2a.2 = ggplot(adata@meta.data, aes(x = Ground_truth, y = CytoTRACE_normandinvert)) + geom_point(color = 'darkgreen', size = 0.5) + ggtitle("CytoTRACE") + xlab("Ground truth potency") + ylab("Inverted potency score") + theme_classic()
fig2a.3 = ggplot(adata@meta.data, aes(x = Ground_truth, y = ccat_normandinvert)) + geom_point(color = 'darkred', size = 0.5) + ggtitle("CCAT") + xlab("Ground truth potency") + ylab("Inverted potency score") + theme_classic()
fig2a.4 = FeaturePlot(adata, features = 'Ground_truth', cols = c('blue','red')) + ggtitle("Ground truth potency")
  #etc. for feature plots of sF, Cyto, and CCAT inverted scores

########################################

#Fig 3 - Performance comparison for UMI data 

  ## four performance metrics box plots
load("../../PerformanceResults/df_performanceNov1822_noregress_genelisttestandcomp_newgini_withgeneset.rda")
fig3a = ggplot(df_all[df_all$iter == 1 & df_all$genelist == 'Full_Regev' & df_all$UMI == 1 & df_all$method %in% c('stemFinder','CytoTRACE','CCAT'),], aes(x = method, y = spear_all)) + geom_boxplot(aes(group = method)) + geom_point(aes(fill = dataset_forplotting), size=2, shape=21,position = position_dodge(0.1)) + theme_clean() + ggtitle("Single-Cell Spearman Correlation") + ylab("Single-cell Spearman Correlation") + xlab("Potency quantification method") + scale_fill_discrete(name = 'Validation Datasets') + scale_colour_discrete(guide = 'none') + theme(axis.text.x=element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), plot.title = element_text(size=20)) 
fig3b = ggplot(df_all[df_all$iter == 1 & df_all$genelist == 'Full_Regev' & df_all$UMI == 1 & df_all$method %in% c('stemFinder','CytoTRACE','CCAT'),], aes(x = method, y = spear_pheno)) + geom_boxplot(aes(group = method)) + geom_point(aes(fill = dataset_forplotting), size=2, shape=21,position = position_dodge(0.1)) + theme_clean() + ggtitle("Phenotypic Spearman Correlation") + ylab("Phenotypic Spearman Correlation") + xlab("Potency quantification method") + scale_fill_discrete(name = 'Validation Datasets') + scale_colour_discrete(guide = 'none') + theme(axis.text.x=element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), plot.title = element_text(size=20)) 
fig3c = ggplot(df_all[df_all$iter == 1 & df_all$genelist == 'Full_Regev' & df_all$UMI == 1 & df_all$method %in% c('stemFinder','CytoTRACE','CCAT'),], aes(x = method, y = auc)) + geom_boxplot(aes(group = method)) + geom_point(aes(fill = dataset_forplotting), size=2, shape=21,position = position_dodge(0.1)) + theme_clean() + ggtitle("Discrimination Accuracy") + ylab("AUC") + xlab("Potency quantification method") + scale_fill_discrete(name = 'Validation Datasets') + scale_colour_discrete(guide = 'none') + theme(axis.text.x=element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), plot.title = element_text(size=20)) 
fig3d = ggplot(df_all[df_all$iter == 1 & df_all$genelist == 'Full_Regev' & df_all$UMI == 1 & df_all$method %in% c('stemFinder','CytoTRACE','CCAT'),], aes(x = method, y = pct.recov)) + geom_boxplot(aes(group = method)) + geom_point(aes(fill = dataset_forplotting), size=2, shape=21,position = position_dodge(0.1)) + theme_clean() + ggtitle("Percent recovery of highly potent cells") + ylab("Percent recovery") + xlab("Potency quantification method") + scale_fill_discrete(name = 'Validation Datasets') + scale_colour_discrete(guide = 'none') + theme(axis.text.x=element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), plot.title = element_text(size=20)) 

  ## memory & run time
load("../../PerformanceResults/runandmem_ubuntus3_periphblood.rda")
fig3e = ggplot(runandmem_copy, aes(x = NumCells, y = MemoryKB)) + geom_line(aes(color = Method), size=2, alpha = 0.5) + geom_point(aes(color = Method), size=3, alpha = 0.5) + theme_clean() + ggtitle("Memory usage") + ylab("Memory (kB)") + xlab("Number of cells") + theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
fig3f = ggplot(runandmem_copy, aes(x = NumCells, y = RunTime.sec.)) + geom_line(aes(color = Method), size=2, alpha = 0.5) + geom_point(aes(color = Method), size=3, alpha = 0.5) + theme_clean() + ggtitle("Run time") + ylab("Run time (sec)") + xlab("Number of cells") + theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))

########################################

#Fig 4 - comparative sF score

load("../../PerformanceResults/df_mostandleastpotentpops_comparativesF_120622.rda")
fig4a = ggplot(df_stem[df_stem$umi == 1,], aes(x = tissue, y = stemFinder_comp)) + geom_boxplot(aes(fill = tissue, color = tissue)) + xlab("Developmental stage") + ylab("Comparative stemFinder score") + scale_color_discrete(name = 'Developmental stage') + scale_fill_discrete(guide = 'none')

load("../../PerformanceResults/df_escs_compscore.rda")
fig4b_mouse = ggplot(df_esc[df_esc$dataset == 'mESC STARR-seq',], aes(x = stemFinder_comp, y = phenotype)) + geom_jitter(aes(color = phenotype), height=0.1, size = 0.1) + geom_boxplot(aes(group = phenotype, color = phenotype), alpha = 0.2, outlier.colour = NA) + theme_classic() + ylab("Sample group") + xlab("Comparative stemFinder score") + ggtitle("Mouse") + scale_color_discrete(name = 'Sample group')
fig4b_human = ggplot(df_esc[df_esc$dataset != 'mESC STARR-seq',], aes(x = stemFinder_comp, y = phenotype)) + geom_jitter(aes(color = phenotype), height=0.1, size = 0.1) + geom_boxplot(aes(group = phenotype, color = phenotype), alpha = 0.2, outlier.colour = NA) + theme_classic() + ylab("Sample group") + xlab("Comparative stemFinder score") + ggtitle("Human") + scale_color_discrete(name = 'Sample group')

########################################

#Fig 5 - CC gene expression in TM BMMC 10X dataset

adata = readRDS("../../ValidationandGeneLists/bmmc_mouse_fromcytotrace_10X_ginicytoandccat.rds")
Idents(adata) = 'Phenotype'

  #5A: GSEA for stem cells vs. all other cells
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
fig5a <- ggplot(dot_df[dot_df$type == 'upregulated',][c(2,4,6,7,10,21,43,46,59),], aes(x = NES, y = fct_reorder(Description, NES))) + 
  geom_point(aes(size=2, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) +
  ggtitle("MSigDB GSEA results") +
  scale_y_discrete(labels = function(y) str_wrap(y, width=45)) +
  theme(axis.text.y = element_text(size=11)) +
  guides(size="none") 

  #5B-D heat maps of binarized scaled CC gene expression
load("../../ValidationandGeneLists/s_genes_mouse.rda")
load("../../ValidationandGeneLists/g2m_genes_mouse.rda")
load("../../ValidationandGeneLists/G1_GSEAMsigDB_mouse.rda")
markers_cc = g1_mouse #set as needed for figB-D

expDat_sub = as.matrix(adata@assays$RNA@scale.data)[markers_cc[markers_cc %in% rownames(adata)],]
sf = adata@meta.data$gini_cellcycle_normandinvert
names(sf) = rownames(adata@meta.data)
sf = sf[order(sf, decreasing = F)]
toplot = expDat_sub[,order(match(colnames(expDat_sub), names(sf)))]
toplot = (as.matrix(toplot > 0))*1 #binarized scaled expression matrix with cells in order of ascending inverted stemFinder score
rows.cor <- cor(t(expDat_sub), use = "pairwise.complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor)) #dendrogram

bar = as.numeric(as.factor(names(sf))) #colors
newcol <- colorRampPalette(c('blue','orange'))
ncols <- length(sf)
colSide <- newcol(ncols)
dev.off()
heatmap(toplot, Colv=NA, Rowv = as.dendrogram(hclust.row), scale = 'none', col = c('black','red'), ColSideColors = colSide, cexRow = 0.5) # distfun=function (y) dist(y,method = "euclidean"))
      #note: to save as TIFF, just run heatmap code after tiff() line, not object 
      #cexRow is rowname font size

  #5E: gene set score vs. inverted sF score 
source("../../stemFinder_code/stemFinder/R/gene_set_score.R")
adata = gene_set_score(gene.set = s_genes_mouse, adata) #returns in meta column called gene.set.score
fig5e_s <- ggplot(adata@meta.data, aes(x = gini_cellcycle_normandinvert, y = gene.set.score)) + geom_point(aes(color = gini_cellcycle_normandinvert)) + scale_color_gradient(low = 'blue', high = 'orange', breaks = c(0, max(adata$gini_cellcycle_normandinvert)), labels = c('Min','Max'), name = 'Inverted stemFinder score') + geom_smooth(method = 'loess', colour = 'blue') + ylab("S phase gene set score") + xlab("Inverted stemFinder potency score") 
    ## note: "gini_cellcycle_normandinvert" was the old column notation for "stemFinder_invert"
    ## etc. for G2M and G1 gene sets

  #5F-G: histograms of CC phase
toplot = adata@meta.data[,c('Phase','gini_cellcycle_normandinvert')]
fig5f <- ggplot(toplot, aes(gini_cellcycle_normandinvert)) + geom_histogram(binwidth = 0.05, aes(fill = Phase)) + ggtitle("Distribution of CC Phase over stemFinder Trajectory") + theme_classic() + ylab("Cell count") + xlab("Inverted stemFinder potency score")
toplot = adata@meta.data[,c('Phase','Ground_truth')]
fig5g <- ggplot(toplot, aes(Ground_truth)) + geom_histogram(binwidth = .5, aes(fill = Phase)) + ggtitle("Distribution of CC Phase over Ground truth potency") + theme_classic() + ylab("Cell count") + xlab("Ground truth potency")

#######################################################

# Fig 6 - murine dentate gyrus example

adata = readRDS("../../ValidationandGeneLists/DentateGyrus_mouse_phenotypes_10X_ginicytoandccat.rds")

  #6A - B: UMAP of cell phenotype (A) and CC phase (B)
fig6a <- DimPlot(adata, group.by = 'Phenotype', cols = c('pink','red','yellow','orange','green','darkgreen','black','grey','blue','darkblue','brown','darkred','purple','magenta','tan')) + ggtitle("Original cell phenotypes")
fig6b <- DimPlot(adata, group.by = 'Phase') + ggtitle("Cell cycle phase")

  #6C - D: UMAP of potency values (GT and sF)
fig6c <- FeaturePlot(adata, features = 'Order', cols = c('blue','red')) + ggtitle("Ground truth potency")
fig6d <- FeaturePlot(adata, features = 'gini_cellcycle_normandinvert', cols = c('blue','red')) + ggtitle("Inverted stemFinder potency")

  #6E: UMAP with highlighted non-selected RGL cells
highlypot.sF <-rownames(adata@meta.data[adata$gini_cellcycle_normandinvert < quantile(adata$gini_cellcycle_normandinvert, (1/length(unique(adata$Order)))),]) #use same threshold in pct.recov() to get highly potent cells per potency calculator
highlypot.cyto <-rownames(adata@meta.data[adata$CytoTRACE_normandinvert < quantile(adata$CytoTRACE_normandinvert, (1/length(unique(adata$Order)))),]) 
highlypot.ccat <-rownames(adata@meta.data[adata$ccat_normandinvert < quantile(adata$ccat_normandinvert, (1/length(unique(adata$Order)))),]) 
truepot <- rownames(adata@meta.data[adata$Order == min(adata$Order),]) #cells with min ground truth potency
computedpot <- union(union(highlypot.sF, highlypot.cyto), highlypot.ccat) #cells marked as highly potent by any calculator
length(truepot[!truepot %in% computedpot]) #num highly potent cells missed by all calculators
fig6e <- DimPlot(adata, cells.highlight = truepot[!truepot %in% computedpot], pt.size = 0.5, sizes.highlight = 0.5) + ggtitle("Non-selected RGL cells") + NoLegend()

  #6F: Violin plots of potency scores by phenotype
fig6f_sF <- VlnPlot(adata, features = 'gini_cellcycle_normandinvert', group.by = 'Phenotype') + ggtitle("Inverted stemFinder potency")
fig6f_cyto <- VlnPlot(adata, features = 'CytoTRACE_normandinvert', group.by = 'Phenotype') + ggtitle("Inverted CytoTRACE potency")
fig6f_ccat <- VlnPlot(adata, features = 'ccat_normandinvert', group.by = 'Phenotype') + ggtitle("Inverted CCAT potency")

  #6G - I: RGL cells: Leiden clusters, potency bin, and CC phase
rgl <- readRDS("../../ValidationandGeneLists/dentategyrus_rglonly.rds") #subsetted RGL cells, re-analyzed using BasicSeuratPipeline.R, clustered with Leiden resolution 0.25
fig6g <- DimPlot(rgl, group.by = 'seurat_clusters') + ggtitle("RGL cells, Leiden clusters")  
rgl$highly.pot <- 'Not marked as highly potent'
rgl@meta.data[colnames(rgl)[colnames(rgl) %in% computedpot],]$highly.pot = 'Marked as highly potent'
fig6h <- DimPlot(rgl, group.by = 'highly.pot') + ggtitle("RGL cells, potency designation")
fig6i <- DimPlot(rgl, group.by = 'Phase') + ggtitle("RGL cells, cell cycle phase")
  
  #6J: RGL cells: number of counts and features per cluster
fig6j_counts <- VlnPlot(rgl, features = 'nCount_RNA', group.by = 'seurat_clusters') + ggtitle("Counts") + xlab("RGL subcluster")
fig6j_feat <- VlnPlot(rgl, features = 'nFeature_RNA', group.by = 'seurat_clusters') + ggtitle("Features") + xlab("RGL subcluster")

  #6K: heat map of DE genes
Idents(rgl) <- 'seurat_clusters'
de.genes <- FindMarkers(rgl, ident.1 = '0', only.pos = F, min.pct = 0, logfc.threshold = 0)
de.genes$p_val_adj_fdr <- p.adjust(p = de.genes$p_val, method = 'fdr') #FDR adjustment (default by Seurat is Bonferroni)
de.genes[c('S100a16','Gfap','Sparcl1','Pla2g7','Mycn','Hes6'),] #genes in heat map
fig6k <- DoHeatmap(adata, features = c('S100a16','Gfap','Sparcl1','Pla2g7','Mycn','Hes6'), group.by = 'Phenotype', group.colors = c('pink','red','yellow','orange','green','darkgreen','black','grey','blue','darkblue','brown','darkred','purple','magenta','tan'), angle = 90, size = 3)

  #6L: revised phenotypes
rev.pheno <- read.csv("../../ValidationandGeneLists/DentateGyrusPhenotype_RevisedAnnotandOrder.csv", header = T, row.names = 1)
adata$Phenotype_revised = rev.pheno$Revised_Phenotype
fig6l <- DimPlot(adata, group.by = 'Phenotype_revised', pt.size = 0.5, cols = c('pink','red','black','yellow','orange','green','darkgreen','black','grey','blue','darkblue','brown','darkred','purple','magenta','tan')) + ggtitle("Revised cell phenotypes")
  
  #6M: GSEA, RGL clust 0 vs. 1
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

#################################################

#Fig 7 - Lineage tracing
adata = readRDS("../../data_from_pc/adata_cospar_seurat_030823.rds")

    #7A-B: UMAP cell type, time point, sF score
fig7a <- DimPlot(adata, group.by = 'state_info') + ggtitle("Cell type")
fig7b <- DimPlot(adata, group.by = 'time_info') + ggtitle("Time point")
fig7c <- FeaturePlot(adata, features = 'stemFinder_invert') + ggtitle("Inverted stemFinder potency")

    #7D: box plot, sF vs number lineages
      ## Run code in LineageTracing_sF_CoSpar.R to get lin_counts
fig7d <- ggplot(lin_counts[,c(2,4)], aes(y = as.numeric(Min_sFinvert_undiff), x = as.numeric(Num_lineages))) + geom_boxplot(aes(group = Num_lineages)) + guides(fill = 'none', color = 'none') + geom_point(aes(color = as.numeric(Min_sFinvert_undiff))) + xlab("Number of downstream lineages") + ylab("Minimum inverted stemFinder score") + ggtitle("Inverted stemFinder vs. number of lineages")

    #7E-F: UMAP marker gene expression
fig7e <- FeaturePlot(adata, features = 'Ngp') 
fig7f <- FeaturePlot(adata, features = 'Mmp8') 

##################################################

# SUPPLEMENTAL FIGURES


# Fig S1: CC gene expression correlates with GT potency

  #S1A-B: GSEA of TFs ranked by stochasticity of expression -- NOT READY YET **********

    ## A: murine BMMC 10X
load("../../ValidationandGeneLists/mmTFs.rda")
adata = readRDS("../../ValidationandGeneLists/bmmc_mouse_fromcytotrace_10X_ginicytoandccat.rds")

        #compute gene expression stochasticity of each TF for each ground truth potency-defined cluster
markers = mmTFs[mmTFs %in% rownames(adata)]
idents = unique(adata$Ground_truth) 
expDat = as.matrix(adata@assays$RNA@scale.data)
sampTab = adata@meta.data

gini_agg = data.frame("cluster" = rep(idents, length(markers)), "gini_sum" = rep(NA, length(idents)*length(markers)), "TF" = rep("", length(idents)*length(markers)))
for (i in idents){
  gini_agg[gini_agg$cluster ==i,]$TF = markers
  neigh = rownames(sampTab[sampTab$Ground_truth == i,]) 
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
        #rank each TF by Pearson correlation of GT potency and gene expression stochasticity
pearsons_agg <- data.frame("Pearson" = rep(NA, length(markers)), "TF" = markers)
for (t in unique(gini_agg$TF)){
  gini_sub = gini_agg[gini_agg$TF == t,]
  pearsons_agg[pearsons_agg$TF == t,]$Pearson = cor.test(gini_sub$cluster, gini_sub$gini_sum, method = 'pearson')$estimate
}
    ## recall: here, gini_sum is high if high stochasticity of expression, aka low GT potency score
    ## multiplied all correlations by -1 because neg correlation is desired for key TFs
pearsons_agg$Pearson = -(pearsons_agg$Pearson)
save(pearsons_agg, file = '../../ValidationandGeneLists/pearsons_agg_tfs_bmmc10x.rda')

        #run GSEA on TFs ranked by Pearson correlation
m_df_c5 = msigdbr(species = "Mus musculus", category = c('C5')) %>%   #GO, KEGG
  dplyr::select(gs_name, gene_symbol)
m_df_h = msigdbr(species = "Mus musculus", category = c('H')) %>%   #hallmark
  dplyr::select(gs_name, gene_symbol)
m_df_c8 = msigdbr(species = "Mus musculus", category = c('C8')) %>%   #hallmark
  dplyr::select(gs_name, gene_symbol)
m_df = rbind(m_df_c5, m_df_h, m_df_c8)
genelist_rf = data.frame("Pearson" = pearsons_agg$Pearson, row.names = pearsons_agg$TF)
genelist_rf = na.omit(genelist_rf)
genelist_rf2 <- genelist_rf %>% mutate(rank = rank(Pearson,  ties.method = "random")) %>%
  arrange(desc(rank)) # handle ties in ranked list
genelist_vect <- as.vector(genelist_rf2$rank) #must be vector for GSEA() function
names(genelist_vect) = rownames(genelist_rf2)
gsea_res = GSEA(genelist_vect, TERM2GENE = m_df, pvalueCutoff=0.05, scoreType = 'pos', verbose = TRUE, nPermSimple=10000) #increased nPermSimple based on error msg suggestion for p value calculation
        #dot plot of top 20 terms ranked by NES
dot_df = as.data.frame(gsea_res[gsea_res$p.adjust < 0.05,])
dot_df$absNES = abs(dot_df$NES)
dot_df = arrange(dot_df, desc(absNES),p.adjust) 
dot_df$type = "upregulated"
dot_df$type[dot_df$NES < 0] = "downregulated"
dot_df$Description = stri_replace_all_fixed(dot_df$Description, "_", " ") #adds space

s1a <- ggplot(dot_df[dot_df$type == 'upregulated',][1:20,], aes(x = NES, y = fct_reorder(Description, NES))) + 
  geom_point(aes(size=2, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) +
  ggtitle("GSEA, Murine BMMC") +
  scale_y_discrete(labels = function(y) str_wrap(y, width=45)) +
  theme(axis.text.y = element_text(size=11)) +
  guides(size="none") 

  #S1C: box plot, single-cell Spearman for stochastic vs. mean CC expression (UMI only)
load("../../PerformanceResults/df_performanceNov1822_noregress_genelisttestandcomp_newgini_withgeneset.rda")
s1c <- ggplot(df_all[df_all$UMI == 1 & df_all$method %in% c('stemFinder','GeneSetScore') & df_all$genelist %in% c('Full_Regev','Full_Regev_Plus_G1','G1','G2M','S'),], aes(x = method, y = spear_all)) + geom_point(position = position_jitterdodge(jitter.width=0.1, dodge.width = 0.5), aes(group = genelist, color = genelist)) + geom_boxplot(aes(group = method), alpha = 0) + theme_clean() + coord_flip() + labs(color = 'Gene list') + xlab("Method") + ylab("Single-Cell Spearman Correlation") + ggtitle("Correlation of gene set score vs. expression stochasticity with ground truth potency") + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14)) + scale_color_discrete(labels = c('S + G2M','S + G2M + G1','G1 only','G2M only','S only'))
stat.test <- df_all[df_all$UMI == 1 & df_all$method %in% c('stemFinder','GeneSetScore') & df_all$genelist %in% c('Full_Regev','Full_Regev_Plus_G1','G1','G2M','S'),] %>% pairwise_t_test(spear_all ~ method, paired = T, p.adjust.method = 'bonferroni')

  #S1D: box plot, performance of stemFinder vs. gene set score 

##
##

#####################################################

  #Supplemental Fig 2: performance in non-UMI data

#####################################################

  #Supplemental Fig 3: robustness heat maps
      ## order of all figures: single-cell Spearman, phenotypic Spearman, AUC, pct.recov

#S3A-D: K robustness
load("../../PerformanceResults/df_robustness_k_041123.rda")
s3a <- ggplot(df_robust_k[df_robust_k$UMI == 1,], aes(x = kratio, y = dataset_forplot, fill = neg_deviation_spearss)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to changes in K - Single-cell Spearman Correlation") + ylab("Validation dataset") + xlab("Deviation from Single-cell Spearman Correlation at Kideal") 
s3b <- ggplot(df_robust_k[df_robust_k$UMI == 1,], aes(x = kratio, y = dataset_forplot, fill = neg_deviation_spearp)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to changes in K - Phenotypic Spearman Correlation") + ylab("Validation dataset") + xlab("Deviation from Phenotypic Spearman Correlation at Kideal") 
s3c <- ggplot(df_robust_k[df_robust_k$UMI == 1,], aes(x = kratio, y = dataset_forplot, fill = neg_deviation_auc)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to changes in K - AUC") + ylab("Validation dataset") + xlab("Deviation from AUC at Kideal") 
s3d <- ggplot(df_robust_k[df_robust_k$UMI == 1,], aes(x = kratio, y = dataset_forplot, fill = neg_deviation_pctrecov)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to changes in K - Percent recovery") + ylab("Validation dataset") + xlab("Deviation from Percent Recovery at Kideal") 

#S3E-H: Cell cycle gene list (GO, Kegg, Regev)
load("../../PerformanceResults/df_robustness_genelist_04202023.rda")
s3e <- ggplot(df_robust_genelist[df_robust_genelist$UMI == 1 & df_robust_genelist$genelist %in% c('Full_Regev','KEGG','GeneOntology'),], aes(x = genelist, y = dataset_forplot, fill = neg_deviation_spearss)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to changes in cell cycle gene list - Single-cell Spearman Correlation") + ylab("Validation dataset") + xlab("Deviation from Single-cell Spearman Correlation") 
s3f <- ggplot(df_robust_genelist[df_robust_genelist$UMI == 1 & df_robust_genelist$genelist %in% c('Full_Regev','KEGG','GeneOntology'),], aes(x = genelist, y = dataset_forplot, fill = neg_deviation_spearp)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to changes in cell cycle gene list - Phenotypic Spearman Correlation") + ylab("Validation dataset") + xlab("Deviation from Phenotypic Spearman Correlation") 
s3g <- ggplot(df_robust_genelist[df_robust_genelist$UMI == 1 & df_robust_genelist$genelist %in% c('Full_Regev','KEGG','GeneOntology'),], aes(x = genelist, y = dataset_forplot, fill = neg_deviation_auc)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to changes in cell cycle gene list - AUC") + ylab("Validation dataset") + xlab("Deviation from AUC") 
s3h <- ggplot(df_robust_genelist[df_robust_genelist$UMI == 1 & df_robust_genelist$genelist %in% c('Full_Regev','KEGG','GeneOntology'),], aes(x = genelist, y = dataset_forplot, fill = neg_deviation_pctrecov)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to changes in cell cycle gene list - Percent recovery") + ylab("Validation dataset") + xlab("Deviation from Percent Recovery") 

#S3I-L: Proportional downsampling of all cells
load("../../PerformanceResults/df_robustness_downsample_allandpotent_042423.rda")
df_robust_ds_all <- df_robust_ds[df_robust_ds$ds_type == 'All cell types',]
s3i <- ggplot(df_robust_ds_all[df_robust_ds_all$UMI == 1,], aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_spearss)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to proportional downsampling - Single-cell Spearman Correlation") + ylab("Validation dataset") + xlab("Deviation from Single-cell Spearman Correlation") 
s3j <- ggplot(df_robust_ds_all[df_robust_ds_all$UMI == 1,], aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_spearp)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to proportional downsampling - Phenotypic Spearman Correlation") + ylab("Validation dataset") + xlab("Deviation from Phenotypic Spearman Correlation") 
s3k <- ggplot(df_robust_ds_all[df_robust_ds_all$UMI == 1,], aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_auc)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to proportional downsampling - AUC") + ylab("Validation dataset") + xlab("Deviation from AUC") 
s3l <- ggplot(df_robust_ds_all[df_robust_ds_all$UMI == 1,], aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_pctrecov)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to proportional downsampling - Percent recovery") + ylab("Validation dataset") + xlab("Deviation from Percent Recovery") 

#S3M-P: Downsampling of most potent population
df_robust_ds_pot <- df_robust_ds[df_robust_ds$ds_type == 'Highly potent only',]
s3m <- ggplot(df_robust_ds_pot[df_robust_ds_pot$UMI == 1,], aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_spearss)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to downsampling of highly potent cells - Single-cell Spearman Correlation") + ylab("Validation dataset") + xlab("Deviation from Single-cell Spearman Correlation") 
s3n <- ggplot(df_robust_ds_pot[df_robust_ds_pot$UMI == 1,], aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_spearp)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to downsampling of highly potent cells - Phenotypic Spearman Correlation") + ylab("Validation dataset") + xlab("Deviation from Phenotypic Spearman Correlation") 
s3o <- ggplot(df_robust_ds_pot[df_robust_ds_pot$UMI == 1,], aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_auc)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to downsampling of highly potent cells - AUC") + ylab("Validation dataset") + xlab("Deviation from AUC") 
s3p <- ggplot(df_robust_ds_pot[df_robust_ds_pot$UMI == 1,], aes(x = ds_ratio, y = dataset_forplot, fill = neg_deviation_pctrecov)) + geom_tile() + scale_fill_viridis(discrete = F, name = 'Deviation in Performance') + ggtitle("Robustness to downsampling of highly potent cells - Percent recovery") + ylab("Validation dataset") + xlab("Deviation from Percent Recovery") 

#S3Q-T: Regev vs. random gene list (box plot)
s3q <- ggplot(df_robust_genelist[df_robust_genelist$UMI == 1 & df_robust_genelist$genelist %in% c('Full_Regev','Random'),], aes(x = genelist, y = spear_all_sF)) + geom_point(aes(color = dataset_forplot)) + geom_boxplot(aes(group = genelist), alpha = 0) + ggtitle("Random vs. cell cycle genes") + xlab("Validation dataset") + ylab("Single-cell Spearman Correlation") + theme_bw() + scale_color_discrete(name = 'Validation Dataset') + scale_x_discrete(labels = c('S + G2M genes','Random'))
s3r <- ggplot(df_robust_genelist[df_robust_genelist$UMI == 1 & df_robust_genelist$genelist %in% c('Full_Regev','Random'),], aes(x = genelist, y = spear_pheno_sF)) + geom_point(aes(color = dataset_forplot)) + geom_boxplot(aes(group = genelist), alpha = 0) + ggtitle("Random vs. cell cycle genes") + xlab("Validation dataset") + ylab("Phenotypic Spearman Correlation") + theme_bw() + scale_color_discrete(name = 'Validation Dataset') + scale_x_discrete(labels = c('S + G2M genes','Random'))
s3s <- ggplot(df_robust_genelist[df_robust_genelist$UMI == 1 & df_robust_genelist$genelist %in% c('Full_Regev','Random'),], aes(x = genelist, y = auc_sF)) + geom_point(aes(color = dataset_forplot)) + geom_boxplot(aes(group = genelist), alpha = 0) + ggtitle("Random vs. cell cycle genes") + xlab("Validation dataset") + ylab("Discrimination Accuracy") + theme_bw() + scale_color_discrete(name = 'Validation Dataset') + scale_x_discrete(labels = c('S + G2M genes','Random'))
s3t <- ggplot(df_robust_genelist[df_robust_genelist$UMI == 1 & df_robust_genelist$genelist %in% c('Full_Regev','Random'),], aes(x = genelist, y = pct.recov_sF)) + geom_point(aes(color = dataset_forplot)) + geom_boxplot(aes(group = genelist), alpha = 0) + ggtitle("Random vs. cell cycle genes") + xlab("Validation dataset") + ylab("Percentage Recovery") + theme_bw() + scale_color_discrete(name = 'Validation Dataset') + scale_x_discrete(labels = c('S + G2M genes','Random'))

#####################################

#Supplementary figure 4: SCN to rescue non-UMI datasets --> removed from paper

#####################################

#Supplementary figure 5: identifying TFs in oligodendrocytes

adata = readRDS("../../ValidationandGeneLists/Oligodendrocytephenotypes_C1mouse_giniandcyto.rds")
load("../../ValidationandGeneLists/sporadtfs_oligodendro.rda") 
sporad_tfs <- finals_withcc[[6]] #TFs with high stochasticity of expression in one phenotype-defined cluster & low in another AND overlap with top DE genes

#S5A: UMAP cell phenotype
s5a <- DimPlot(adata, group.by = 'Phenotype', cols = c('blue','brown','orange','navy','green','black','darkgreen','purple','magenta','tan','pink','red')) + ggtitle("Cell phenotype")

#S5B-C: violin plot of inverted potency score
s5b <- VlnPlot(adata, features = 'gini_cellcycle_normandinvert', group.by = 'Phenotype', cols = c('blue','brown','orange','navy','green','black','darkgreen','purple','magenta','tan','pink','red')) + ggtitle("Inverted stemFinder score") + xlab("Cell phenotype")
s5c_cyto <- VlnPlot(adata, features = 'CytoTRACE_normandinvert', group.by = 'Phenotype', cols = c('blue','brown','orange','navy','green','black','darkgreen','purple','magenta','tan','pink','red')) + ggtitle("Inverted CytoTRACE score")+ xlab("Cell phenotype")
s5c_ccat <- VlnPlot(adata, features = 'ccat_normandinvert', group.by = 'Phenotype', cols = c('blue','brown','orange','navy','green','black','darkgreen','purple','magenta','tan','pink','red')) + ggtitle("Inverted CCAT score")+ xlab("Cell phenotype")

#S5D: UMAP of sF score with sporadically expressed TFs as input
adata = run_stemFinder(adata, nn = adata@graphs$RNA_nn, k = round(sqrt(ncol(adata))), thresh = 0, markers = sporad_tfs)
s5d <- FeaturePlot(adata, features = 'stemFinder_invert', cols = c('blue','red')) + ggtitle("Inverted stemFinder score, Sporadically expressed TFs")

#S5E: Performance of CC genes vs. sporadically expressed TFs for all UMI datasets 
load("../../PerformanceResults/df_performanceNov1822_noregress_genelisttestandcomp_newgini_withgeneset.rda")
s5e <- ggplot(df_all[df_all$iter == 1 & df_all$method == 'stemFinder' & df_all$genelist %in% c('Full_Regev','SporadTFonly') & df_all$UMI == 1,], aes(x = spear_all, y = genelist)) + geom_boxplot(aes(group = genelist)) + geom_point(aes(fill = dataset_forplotting), size=2, shape=21,position = position_dodge(0.1)) + theme_bw() + ggtitle("Performance of Sporadically Expressed TFs vs. Cell Cycle Genes") + geom_line(aes(group = dataset_forplotting, color = dataset_forplotting, alpha = 0.8), position = position_dodge(0.1)) + xlab("Single-cell Spearman Correlation") + ylab("Input Marker Gene List") + scale_y_discrete(labels = c('SporadTFonly' = 'Sporadically Expressed TFs','Full_Regev' = 'Regev Cell Cycle Genes')) + scale_fill_discrete(name = 'Validation Datasets') + scale_colour_discrete(guide = 'none') + theme(axis.text.x=element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), plot.title = element_text(size=20)) + coord_flip()

