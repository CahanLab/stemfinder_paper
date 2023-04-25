### GSEA analysis with clusterProfiler

# Setup
library(plyr)
library(enrichplot)
library(forcats)
library(stringr)
library(stringi)
library(clusterProfiler)
library(msigdbr)


    # Get input gene list: need ranked list of ALL genes in dataset
          ## Set cutoffs of min.pct and logfc.threshold to zero to get all genes
markers = FindMarkers(adata, ident.1 = "CIAtreat", min.pct = 0, logfc.threshold=0, only.pos=FALSE)
dim(markers[markers$p_val_adj < 0.05,])

    # Get reference database
cat = "C5"
m_df = msigdbr(species = "Homo sapiens", category = cat, subcategory = "BP") %>% 
  dplyr::select(gs_name, gene_symbol)

    # Reformat input gene list - order by log fold change
genelist = markers
genelist_rf = genelist$avg_log2FC
names(genelist_rf) = rownames(genelist)
genelist_rf = na.omit(genelist_rf)
genelist_rf = sort(genelist_rf, decreasing=TRUE)
length(genelist_rf)

    # Perform GSEA
gsea_res = GSEA(genelist_rf, TERM2GENE = m_df, pvalueCutoff=0.05, verbose = TRUE)

    # Dot plot 
n = 20 #number of terms to plot
dot_df = as.data.frame(gsea_res[gsea_res$p.adjust < 0.05,])
dot_df$absNES = abs(dot_df$NES)
dot_df = arrange(dot_df, desc(absNES),p.adjust)  #reorder by NES then by p value
dot_df$type = "upregulated"
dot_df$type[dot_df$NES < 0] = "downregulated"
dot_df$Description = stri_replace_all_fixed(dot_df$Description, "_", " ") #adds space

ggplot(dot_df[1:n,], aes(x = NES, y = fct_reorder(Description, NES))) + 
  geom_point(aes(size=2, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) +
  ggtitle("MSigDB GSEA results") +
  scale_y_discrete(labels = function(y) str_wrap(y, width=45)) +
  theme(axis.text.y = element_text(size=11)) +
  guides(size="none") 
