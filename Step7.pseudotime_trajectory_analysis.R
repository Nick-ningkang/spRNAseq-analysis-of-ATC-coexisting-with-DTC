######################################################################
###############Pseudotime analysis (monocle3)#########################
######################################################################
# =============================================
# Pseudotime trajectory inference using Monocle3
# =============================================

library(monocle3)
library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)

setwd('/datadisk1/person/NICK/spRNA/pseudotime_analysis')  # Set working directory

# ---------------------------------------------------
# Load Seurat object (after spatial niche annotation)
# ---------------------------------------------------
spRNA <- readRDS("/datadisk1/person/NICK/spRNA/data/after_nicheidentified_spRNA.rds")

# ---------------------------------------------------
# Build a Monocle3 cell_data_set (CDS) object
# ---------------------------------------------------
data <- spRNA@assays$SCT@counts                       # Raw count matrix from SCT assay
cell_metadata <- spRNA@meta.data                      # Cell-level metadata
gene_annotation <- data.frame(gene_short_name = rownames(data))  # Gene annotations
rownames(gene_annotation) <- rownames(data)

# Create Monocle3 CDS
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# ---------------------------------------------------
# Preprocessing and dimensionality reduction
# ---------------------------------------------------
cds <- preprocess_cds(cds, num_dim = 50)               # PCA and scaling
cds <- align_cds(cds, alignment_group = "sample")      # Batch correction across samples
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# ---------------------------------------------------
# Optional: Replace Monocle UMAP with integrated UMAP from Seurat
# ---------------------------------------------------
cds_embed <- cds@int_colData$reducedDims$UMAP
int_embed <- Embeddings(spRNA, reduction = "umap")
int_embed <- int_embed[rownames(cds_embed), ]
cds@int_colData$reducedDims$UMAP <- int_embed

# Plot UMAP with niche annotations
p <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "niche0.3") + ggtitle("Integrated UMAP")
print(p)

# ---------------------------------------------------
# Learn trajectory graph
# ---------------------------------------------------
cds <- cluster_cells(cds, cluster_method = "louvain")  # Cluster before graph learning
plot_cells(cds, color_cells_by = "niche0.3")           # Visualize clusters
cds <- learn_graph(cds)                                # Infer trajectory graph

# Plot trajectory without labels
plot_cells(cds, color_cells_by = "niche0.3", label_groups_by_cluster = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE)

# ---------------------------------------------------
# Order cells in pseudotime (selecting root node)
# ---------------------------------------------------
cds <- order_cells(cds)  # Manual root selection

# Plot pseudotime coloring
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE)

# ---------------------------------------------------
# Set "CITED1+DTC" cells as trajectory root
# ---------------------------------------------------
myselect <- function(cds, select.classify, my_select) {
  cell_ids <- which(colData(cds)[, select.classify] == my_select)
  closest_vertex <- as.matrix(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
  return(root_pr_nodes)
}

cds <- order_cells(cds, root_pr_nodes = myselect(cds, select.classify = 'niche0.3', my_select = "CITED1+DTC"))

# Plot updated pseudotime
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE)

# ---------------------------------------------------
# Plot with manual color scale for niche identity
# ---------------------------------------------------
plot_cells(cds, color_cells_by = "niche0.3", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE,
           graph_label_size = 2) +
  scale_color_manual(values = c(
    'CITED1+DTC' = "#4b6aa8", 'S100A1+DTC' = "#ece399", 'DNAJA4+DTC' = "#2d3462", 'NPW+DTC' = "#c376a7",
    'BMP8A+DHGTC' = "#408444", 'SIGLEC6+DHGTC' = "#61bada", 'SLC34A21+DHGTC' = '#b05545',
    'SERPINE1+ATC' = "#df5734", 'KRT5+ATC' = "#9d3b62"))

# ---------------------------------------------------
# Export plots to PDF
# ---------------------------------------------------
pdf("monocle3_umap.pdf", width = 4.5, height = 3)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 2)
plot_cells(cds, color_cells_by = "niche0.3", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 2) +
  scale_color_manual(values = c(
    'CITED1+DTC' = "#4b6aa8", 'S100A1+DTC' = "#ece399", 'DNAJA4+DTC' = "#2d3462", 'NPW+DTC' = "#c376a7",
    'BMP8A+DHGTC' = "#408444", 'SIGLEC6+DHGTC' = "#61bada", 'SLC34A21+DHGTC' = '#b05545',
    'SERPINE1+ATC' = "#df5734", 'KRT5+ATC' = "#9d3b62"))
plot_cells(cds, color_cells_by = "region_intergration", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 2) +
  scale_color_manual(values = c("#E64B35FF", "#F39B7FFF", "#00A087FF"))
dev.off()

# Export pseudotime-only plot
pdf("monocle3_pseudotime.pdf", width = 4.5, height = 3)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 2)
dev.off()

# Export niche-level plot
pdf("monocle3_niche0.3.pdf", width = 5, height = 3)
plot_cells(cds, color_cells_by = "niche0.3", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 2) +
  scale_color_manual(values = c(
    'CITED1+DTC' = "#4b6aa8", 'S100A1+DTC' = "#ece399", 'DNAJA4+DTC' = "#2d3462", 'NPW+DTC' = "#c376a7",
    'BMP8A+DHGTC' = "#408444", 'SIGLEC6+DHGTC' = "#61bada", 'SLC34A21+DHGTC' = '#b05545',
    'SERPINE1+ATC' = "#df5734", 'KRT5+ATC' = "#9d3b62"))
dev.off()

# Export region integration plot
pdf("monocle3_region_intergration.pdf", width = 5, height = 3)
plot_cells(cds, color_cells_by = "region_intergration", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 2) +
  scale_color_manual(values = c("#E64B35FF", "#F39B7FFF", "#00A087FF"))
dev.off()

# ---------------------------------------------------
# Save the Monocle3 object
# ---------------------------------------------------
saveRDS(cds, "cds_monocle3.rds")  # Save CDS object for future reuse



# =============================================
# Monocle3 Data Visualization
# =============================================

# Load necessary libraries
library(monocle3)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# Step 1: Identify pseudotime-related genes from Monocle3 trajectory
Track_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)
genes <- rownames(subset(Track_genes, q_value < 0.01 & morans_I > 0.2))

# Step 2: Extract expression matrix for significant genes ordered by pseudotime
plot_matrix <- monocle3::exprs(cds)[match(genes, rownames(rowData(cds))), order(pseudotime(cds))]

# Step 3: Smooth and Z-score normalize expression profiles
plot_matrix <- t(apply(plot_matrix, 1, function(x) smooth.spline(x, df = 3)$y))
plot_matrix <- t(apply(plot_matrix, 1, function(x) (x - mean(x)) / sd(x)))
rownames(plot_matrix) <- genes

# Step 4: Bin expression matrix by pseudotime for heatmap (average every N columns)
cutColumn_Means <- function(data_exp, cut) {
  plot_matrix_combin <- list()
  nums <- ncol(data_exp) / cut
  
  if (nums - round(nums, 0) == 0) {
    for (i in seq(1, ncol(data_exp), cut)) {
      A <- rowMeans(data_exp[, i:(i + cut - 1)])
      plot_matrix_combin[[length(plot_matrix_combin) + 1]] <- A
    }
  } else {
    for (i in seq(1, ncol(data_exp) - cut, cut)) {
      A <- rowMeans(data_exp[, i:(i + cut - 1)])
      plot_matrix_combin[[length(plot_matrix_combin) + 1]] <- A
    }
    last <- rowMeans(data_exp[, (max(seq(1, ncol(data_exp) - cut, cut)) + cut):ncol(data_exp)])
    plot_matrix_combin[[length(plot_matrix_combin) + 1]] <- last
  }
  
  plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
  rownames(plot_matrix_combin) <- rownames(data_exp)
  colnames(plot_matrix_combin) <- seq_len(ncol(plot_matrix_combin))
  return(plot_matrix_combin)
}

plot_test <- cutColumn_Means(plot_matrix, cut = 25)

# Step 5: Heatmap with hierarchical clustering
callback <- function(hc, mat) {
  sv <- svd(t(mat))$v[, 1]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

# Basic heatmap clustering
p1 <- pheatmap(plot_test, cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE,
               show_colnames = FALSE, cutree_rows = 4, clustering_method = "ward.D2",
               color = colorRampPalette(c('#1A5592', 'white', '#B83D3D'))(100),
               clustering_callback = callback)

# Annotate heatmap rows with cluster assignments
annotation_row <- data.frame(Cluster = factor(cutree(p1$tree_row, 4)))
rownames(annotation_row) <- rownames(plot_test)
rowcolor <- pal_jama()(4)
names(rowcolor) <- as.character(1:4)
ann_colors <- list(Cluster = rowcolor)

# Final heatmap with annotation
p3 <- pheatmap(plot_test, cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE,
               annotation_row = annotation_row, annotation_colors = ann_colors,
               cutree_rows = 4, clustering_method = "ward.D2",
               color = colorRampPalette(c('#1A5592', 'white', '#B83D3D'))(100),
               clustering_callback = callback,
               main = "Pseudotime")

# Step 6: Manually label genes of interest
gene <- c("TG", "DIO1", "DIO2", "NKX2-1", "PAX8", "FOXE1", "DUOX2", "DUOXA1", "SLC5A3", "SLC5A5", "SLC16A2",
          "COL4A1", "COL4A2", "C1QA", "C1QB", "C1QC", "C3", "C3AR1", "CD14", "CD68", "TYROBP", "CSF1", "CSF1R",
          "MARCO", "MSR1", "MRC1", "MRC2", "CLEC7A", "VSIG4", "TLR2", "TYMP", "CXCL8", "CXCL10", "CXCL11",
          "CXCL13", "CCL2", "CCL5", "CCL13", "TNFSF10", "IL7R", "TRAC", "TRDC", "IGHA1", "IGHG1", "IGHG3",
          "FCER1G", "FCGR2A", "FCGR3A", "TNFRSF11B")

source("add.flag.R")
p5 <- add.flag(p3, kept.labels = gene, repel.degree = 0.2)
ggsave("pheatmap.pdf", p5, width = 8, height = 8)

# Step 7: GO enrichment analysis for each gene module
module_gene <- as.data.frame(cutree(p3$tree_row, k = 4))
colnames(module_gene) <- "Module"
module_gene$gene <- rownames(module_gene)

Module_GO <- data.frame()

for (i in unique(module_gene$Module)) {
  data <- subset(module_gene, Module == i)
  df <- bitr(data$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  go <- enrichGO(gene = unique(df$ENTREZID),
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = TRUE)
  
  if (nrow(go@result) > 0) {
    go@result$cluster <- i
    Module_GO <- rbind(Module_GO, go@result)
  }
}

# Save significant GO terms
Module_GO <- Module_GO[Module_GO$qvalue <= 0.05, ]
write.csv(Module_GO, file = "Module_GO.csv")





######################################################################
###############Pseudotime analysis (Slingshot)#########################
######################################################################
# Load required libraries
suppressPackageStartupMessages({
  library(slingshot)
  library(SingleCellExperiment)
  library(qs)
  library(tidyverse)
  library(RColorBrewer)
})

# Load filtered Seurat object
spRNA <- readRDS("/datadisk1/person/NICK/spRNA/data/after_nicheidentified_spRNA.rds")

# Subset specific niches for trajectory analysis
niches_to_use <- c("CITED1+DTC", "S100A1+DTC", "DNAJA4+DTC", "NPW+DTC",
                   "BMP8A+DHGTC", "SIGLEC6+DHGTC", "SLC34A21+DHGTC",
                   "SERPINE1+ATC", "KRT5+ATC")
spRNAsub <- subset(spRNA, niche0.3 %in% niches_to_use)

# Extract highly variable genes for dimensionality reduction and trajectory
scale.data <- spRNAsub@assays$SCT@scale.data
scale.gene <- rownames(scale.data)
counts <- spRNAsub@assays$SCT@counts[scale.gene, ]

# Create SingleCellExperiment object
sim <- SingleCellExperiment(assays = List(counts = counts))

# Add UMAP embedding
umap <- spRNAsub@reductions$umap@cell.embeddings
colnames(umap) <- c('UMAP-1', 'UMAP-2')
reducedDims(sim) <- SimpleList(UMAP = umap)

# Add metadata
meta <- spRNAsub@meta.data
colData(sim)$region <- meta$region_intergration
colData(sim)$celltype <- meta$niche0.3

# Run Slingshot lineage inference
sim <- slingshot(sim,
                 clusterLabels = sim$celltype,
                 reducedDim = 'UMAP',
                 start.clus = "CITED1+DTC",
                 end.clus = "KRT5+ATC")

# Save result
saveRDS(sim, "LZY/spRNA/trajectory/slingshot/spRNAsub_slingshot_sim.RDS")

# Plot UMAP + lineages
pdf("LZY/spRNA/trajectory/slingshot/spRNAsub_slingshot_All.pdf", width = 6, height = 6)
cell_colors <- setNames(brewer.pal(length(unique(sim$celltype)), "Set1"), unique(sim$celltype))
plot(reducedDims(sim)$UMAP, pch = 16, col = cell_colors[sim$celltype], asp = 1)
lines(SlingshotDataSet(sim), lwd = 2, col = brewer.pal(9, "Set1"))
dev.off()

# Visualize pseudotime gradient
colors <- colorRampPalette(brewer.pal(11, 'Spectral')[-6])(100)
plotcol <- colors[cut(sim$slingPseudotime_1, breaks = 100)]
plotcol[is.na(plotcol)] <- "lightgrey"
pdf("LZY/spRNA/trajectory/slingshot/spRNAsub_slingshot_lineage1.pdf", width = 6, height = 6)
plot(reducedDims(sim)$UMAP, col = plotcol, pch = 16, asp = 1)
lines(SlingshotDataSet(sim), lwd = 1, col = brewer.pal(9, "Set1"))
dev.off()

# Extract cells from lineage 1
coldata <- data.frame(celltype = sim$celltype,
                      region = sim$region,
                      pseudotime = sim$slingPseudotime_1)
filter_cell <- rownames(coldata)[!is.na(coldata$pseudotime)]
counts <- assays(sim)$counts[, filter_cell]

# Randomly sample 2000 cells for downstream fitting
set.seed(123)
scell <- sample(colnames(counts), size = 2000)
counts <- counts[, scell]

# Construct SingleCellExperiment for sampled cells
filter_sim <- SingleCellExperiment(assays = List(counts = counts))
filter_sim@colData <- colData(sim)[scell, 1:3]
reducedDims(filter_sim) <- SimpleList(UMAP = reducedDims(sim)[scell, ])

# K-means clustering for lineage inference
set.seed(4)
kmeans_res <- kmeans(reducedDims(filter_sim)$UMAP, centers = 4)$cluster
filter_sim$kmeans <- kmeans_res

# Run Slingshot again on clustered data
filter_sim <- slingshot(filter_sim,
                        clusterLabels = 'kmeans',
                        reducedDim = 'UMAP',
                        start.clus = "2",
                        end.clus = "1")

# Save processed object
saveRDS(filter_sim, "LZY/spRNA/trajectory/slingshot/filter_sim_K_Means_slingPseudotime.RDS")

# Run tradeSeq to model gene expression along pseudotime
library(tradeSeq)
crv <- SlingshotDataSet(filter_sim)
counts <- assays(filter_sim)$counts

# Choose number of knots (basis functions)
icMat <- evaluateK(counts = counts, sds = crv, k = 3:10, nGenes = 500)
saveRDS(icMat, "LZY/spRNA/trajectory/slingshot/spRNAsub_icMat.RDS")

# Fit GAM with chosen knots (e.g., 6)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
sce <- fitGAM(counts = counts,
              pseudotime = pseudotime,
              cellWeights = cellWeights,
              nknots = 6)
saveRDS(sce, "LZY/spRNA/trajectory/slingshot/spRNAsub_sce.RDS")

# Test for differential expression
assoRes <- associationTest(sce)
startRes <- startVsEndTest(sce)
write.csv(startRes, "LZY/spRNA/trajectory/slingshot/spRNAsub_startRes.csv")

# Plot selected top DE genes along pseudotime
sigGenes <- names(sce)[order(startRes$waldStat, decreasing = TRUE)[1:6]]
sigGenes <- c(sigGenes, "PDCD4", "TYMP", "CD163")
exp_mat <- log2(assays(sce)$counts[sigGenes, ] + 1)
plot_data <- data.frame(colData(sce), t(exp_mat))

# Plot smoothers
library(ggpubr)
library(patchwork)
colors_use <- c("#4b6aa8", "#ece399", "#2d3462", "#c376a7", "#408444", "#61bada",
                "#b05545", "#df5734", "#9d3b62")
plot_list <- lapply(sigGenes, function(gene) {
  ggscatter(plot_data, x = 'pseudotime.Lineage1', y = gene,
            color = 'celltype', size = 0.6) +
    geom_smooth(se = FALSE, color = 'orange') +
    scale_color_manual(values = colors_use) +
    theme_bw() +
    theme(legend.position = 'none')
})
wrap_plots(plot_list)
ggsave("LZY/spRNA/trajectory/slingshot/spRNAsub_select_genes.pdf", width = 9, height = 8)




