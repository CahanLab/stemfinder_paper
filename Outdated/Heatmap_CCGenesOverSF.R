

## Base stats package

  #Prepare scaled expression matrix with cells in order of ascending inverted stemFinder score
expDat_sub = as.matrix(adata@assays$RNA@scale.data)[markers[markers %in% rownames(adata)],]
sf = adata@meta.data$gini_cellcycle_normandinvert
names(sf) = rownames(adata@meta.data)
sf = sf[order(sf, decreasing = F)]
toplot = expDat_sub[,order(match(colnames(expDat_sub), names(sf)))]
toplot = (as.matrix(toplot > 0))*1

  #Dendrogram
rows.cor <- cor(t(expDat_sub), use = "pairwise.complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))

  #Color bar
library(RColorBrewer)
#bar = as.numeric(as.factor(substr(names(sf), 1, 1)))
#colSide = brewer.pal(length(unique(bar)), "Set1")[bar]

bar = as.numeric(as.factor(names(sf)))
newcol <- colorRampPalette(c('blue','orange'))
ncols <- length(sf)
colSide <- newcol(ncols)


  #Heat map 
heatmap(toplot, Colv=NA, Rowv = as.dendrogram(hclust.row), scale = 'none', col = c('black','red'), ColSideColors = colSide, cexRow = 0.5) # distfun=function (y) dist(y,method = "euclidean"))
    ##cexRow is rowname font size
#legend(x = "bottomright", title = "Scaled Gene Expression", legend = c("<= 0", "> 0"), cex = 0.8, fill = c('black','red'))


## Alternative: ggplot2

  #Prepare data frame
library(reshape2)
expMelt = melt(expDat_sub, na.rm = T)
colnames(expMelt) = c('Gene','Cell','Scaled_Expression')
expMelt = expMelt[order(match(expMelt$Cell, names(sf))),]

ggplot(data = expMelt, mapping = aes(x = Cell, y = Gene, fill = Scaled_Expression)) + geom_tile() + theme(axis.title.x = element_blank()) + scale_fill_gradient2(low = 'black', mid = 'black', high = 'red', midpoint = 0, breaks = seq(-20, 0, 20))

