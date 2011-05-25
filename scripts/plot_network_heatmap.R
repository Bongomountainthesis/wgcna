#!/usr/local/bin/Rscript

require("WGCNA")

options(stringsAsFactors = FALSE)

load("results/E.RData")
load("results/WCGNA_part1.RData")
load("results/results.RData")
load("data/TOM.RData")

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function

#postscript(file="TOMplot_again.ps",horizontal=FALSE)
#TOMplot(plotTOM, geneTree, moduleColors,setLayout=TRUE, main = "Network heatmap plot")
#dev.off()

labeltree <- as.character(moduleColors)
labelrow <- as.character(moduleColors)
labelrow[geneTree$order[length(labeltree):1]] <- labelrow[geneTree$order]
options(expressions = 10000)
geneTree$height <- (geneTree$height - min(geneTree$height))/(1.15 *
            (max(geneTree$height) - min(geneTree$height)))

postscript(file="network_heatmap_heat2.ps",horizontal=FALSE)
heatmap(plotTOM, 
	Rowv = NULL, 
	Colv = NULL,
        scale = "none", 
	revC = T, 
	ColSideColors = as.character(labeltree),
        RowSideColors = as.character(labelrow), 
	labRow = F,
        labCol = F, 
	col = heat.colors(100)
)
dev.off()

#dendro <- as.dendrogram(geneTree, hang=0.1)

postscript(file="network_heatmapFULL.ps",horizontal=FALSE)
heatmap(as.matrix(plotTOM), Rowv = as.dendrogram(geneTree, hang=0.1), Colv = as.dendrogram(geneTree, hang=0.1),
                scale = "none", revC = T, ColSideColors = as.character(labeltree),
                RowSideColors = as.character(labelrow), labRow = F,
                labCol = F, col = heat.colors(100),setLayout=T
)
dev.off()

function (dissim, dendro, colors = NULL, colorsLeft = colors,
    terrainColors = FALSE, setLayout = TRUE, ...)
{
    if (is.null(colors))
        colors = rep("white", dim(as.matrix(dissim))[[1]])
    if (is.null(colorsLeft))
        colorsLeft = colors
    nNodes = length(colors)
    if (nNodes < 2) {
        warning("You have only 1 or 2 genes in TOMplot. No plot will be produced")
    }
    else {
        if (nNodes != length(colorsLeft))
            stop("ERROR: number of (top) color labels does not equal number of left color labels")
        if (nNodes != dim(dissim)[[1]])
            stop(paste("ERROR: number of color labels does not equal number of nodes in dissim.\n",
                "     nNodes != dim(dissim)[[1]] "))
        labeltree = as.character(colors)
        labelrow = as.character(colorsLeft)
        labelrow[dendro$order[length(labeltree):1]] = labelrow[dendro$order]
        options(expressions = 10000)
        dendro$height = (dendro$height - min(dendro$height))/(1.15 *
            (max(dendro$height) - min(dendro$height)))
        if (terrainColors) {
            .heatmap(as.matrix(dissim), Rowv = as.dendrogram(dendro,
                hang = 0.1), Colv = as.dendrogram(dendro, hang = 0.1),
                scale = "none", revC = T, ColSideColors = as.character(labeltree),
                RowSideColors = as.character(labelrow), labRow = F,
                labCol = F, col = terrain.colors(100), setLayout = setLayout,
                ...)
        }
        else {
            .heatmap(as.matrix(dissim), Rowv = as.dendrogram(dendro,
                hang = 0.1), Colv = as.dendrogram(dendro, hang = 0.1),
                scale = "none", revC = T, ColSideColors = as.character(labeltree),
                RowSideColors = as.character(labelrow), labRow = F,
                labCol = FALSE, setLayout = setLayout, ...)
        }
    }
}
