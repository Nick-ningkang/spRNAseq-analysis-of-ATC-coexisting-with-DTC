################################################################################
############Pathway analysis and tumor heterogeneity analysis###################
################################################################################
# ==== Load Required Libraries ====
library(Seurat)
library(COSG)
library(data.table)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)

# ==== Set Paths ====
input_rds_path <- "/datadisk1/person/NICK/spRNA/data/after_nicheidentified_spRNA.rds"
output_dir <- "LZY/spRNA/上皮异质性"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==== Step 1: Load and Subset Data ====
spRNA <- readRDS(input_rds_path)
niches_of_interest <- c("CITED1+DTC", "S100A1+DTC", "DNAJA4+DTC", "NPW+DTC",
                        "BMP8A+DHGTC", "SIGLEC6+DHGTC", "SLC34A21+DHGTC",
                        "SERPINE1+ATC", "KRT5+ATC")
spRNAsub <- subset(spRNA, niche0.3 %in% niches_of_interest)
Idents(spRNAsub) <- "niche0.3"

# ==== Step 2: Identify Marker Genes using COSG ====
marker_cosg <- cosg(spRNAsub,
                    groups = 'all',
                    assay = 'SCT',
                    slot = 'data',
                    mu = 1,
                    n_genes_user = 100)
write.csv(marker_cosg$names, file.path(output_dir, "diff_cellype_genes.csv"), row.names = FALSE)

# ==== Step 3: Define Custom Heatmap Function ====
plot_marker_heatmap <- function(top_n = 5, out_prefix = "Heatmap_celltype_top") {
  # Extract top-N marker genes
  cgs_df <- marker_cosg$names[1:top_n, ]
  top_genes <- unlist(cgs_df)
  
  # Average expression matrix
  avg_expr <- AverageExpression(spRNAsub,
                                features = top_genes,
                                assay = 'SCT',
                                group.by = 'niche0.3',
                                slot = 'data') %>%
    data.frame() %>%
    as.matrix()
  
  # Clean column names
  colnames(avg_expr) <- str_remove(colnames(avg_expr), "SCT.")
  celltype_colors <- c("#4b6aa8", "#ece399", "#2d3462", "#c376a7",
                       "#408444", "#61bada", "#b05545", "#df5734", "#9d3b62")
  names(celltype_colors) <- colnames(avg_expr)
  
  # Z-score normalization across genes
  zscore_mat <- t(scale(t(avg_expr), center = TRUE, scale = TRUE))
  
  # Color scale
  col_fun <- colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))
  
  # Top annotation
  column_ha <- HeatmapAnnotation(celltype = colnames(zscore_mat),
                                 col = list(celltype = celltype_colors))
  
  # Plot and save heatmap
  pdf(file.path(output_dir, sprintf("%s%d_marker.pdf", out_prefix, top_n)),
      width = 8, height = ifelse(top_n == 5, 8, 6))
  Heatmap(zscore_mat,
          name = "Z-score",
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          row_title = "Marker genes",
          column_title = "Celltype",
          row_names_gp = gpar(fontface = 'italic', fontsize = 10),
          row_names_side = 'left',
          border = TRUE,
          rect_gp = gpar(col = "white", lwd = 1),
          column_names_side = 'top',
          column_names_rot = 45,
          top_annotation = column_ha,
          col = col_fun)
  dev.off()
}

# ==== Step 4: Generate Heatmaps ====
plot_marker_heatmap(top_n = 5)  # Top 5 marker genes
plot_marker_heatmap(top_n = 3)  # Top 3 marker genes


# ================================================
# GSVA Analysis using Hallmark Gene Sets (Simplified and Annotated)
# ================================================

# ---- Load Required Libraries ----
library(Seurat)
library(GSVA)
library(msigdbr)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)

# ---- Set Global Options ----
options(stringsAsFactors = FALSE)

# ---- Define Output Directory ----
output_dir <- "LZY/spRNA/上皮异质性"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Prepare Hallmark Gene Sets ----
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol) %>%
  split(~ gs_name)

# ---- Calculate Average Gene Expression ----
Idents(spRNAsub) <- "niche0.3"
avg_expr <- AverageExpression(spRNAsub, assays = "SCT", slot = "data")[[1]]
avg_expr <- avg_expr[rowSums(avg_expr) > 0, ] # keep non-zero genes
expr_mat <- as.matrix(avg_expr)

# ---- Run GSVA Analysis ----
gsva_res <- gsva(expr_mat, hallmark_sets, method = "gsva", kcdf = "Gaussian")
rownames(gsva_res) <- str_replace(rownames(gsva_res), "HALLMARK_", "")

# ---- Plot Heatmap ----
p <- pheatmap(gsva_res,
              show_colnames = TRUE,
              scale = "row",
              angle_col = 45,
              cluster_cols = FALSE,
              color = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))

# ---- Save Plot ----
ggsave(file.path(output_dir, "GSVA.Hallmark.enrich.pheatmap.pdf"), p, width = 6, height = 9)

cat("GSVA Hallmark enrichment analysis completed.\n")

#### KEGG Pathway Enrichment Analysis (Celltype-specific) ####

# Load required libraries
library(clusterProfiler)
library(Seurat)
library(dplyr)
library(ggplot2)

# Identify cell-type-specific markers
sce.markers <- FindAllMarkers(
  object = spRNAsub,
  only.pos = TRUE,
  min.pct = 0.25,
  thresh.use = 0.25
)
write.csv(sce.markers, file = "LZY/spRNA/上皮异质性/sce.markers.clusters.csv")

# Filter markers with adjusted p-value < 0.05
markers <- read.csv("LZY/spRNA/上皮异质性/sce.markers.clusters.csv") |>
  group_by(cluster) |>
  filter(p_val_adj < 0.05) |>
  ungroup()

# Convert gene symbols to ENTREZ IDs
gid <- bitr(unique(markers$gene), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
markers <- full_join(markers, gid, by = c('gene' = 'SYMBOL'))

# Perform KEGG enrichment analysis by cell type
kegg_result <- compareCluster(ENTREZID ~ cluster, data = markers, fun = 'enrichKEGG')
write.csv(kegg_result@compareClusterResult, file = "LZY/spRNA/上皮异质性/spRNAsub.KEGG.markers.clusters.csv")

# Dotplot of KEGG results (optional quick view)
dotplot(kegg_result, label_format = 30) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Load manually curated significant KEGG terms for barplot
dt <- read.csv("LZY/spRNA/上皮异质性/selecet_KEGG.csv", header = TRUE)

# Format factors for plotting
celltype_levels <- c("CITED1+DTC","S100A1+DTC","DNAJA4+DTC","NPW+DTC",
                     "BMP8A+DHGTC","SIGLEC6+DHGTC","SLC34A21+DHGTC",
                     "SERPINE1+ATC","KRT5+ATC")
dt$Description <- factor(dt$Description, levels = rev(dt$Description))
dt$Cluster <- factor(dt$Cluster, levels = celltype_levels)

# Define color scheme for each celltype
mycol <- c("#4b6aa8", "#ece399", "#2d3462", "#c376a7", "#408444",
           "#61bada", "#b05545", "#df5734", "#9d3b62")

# Create KEGG enrichment barplot with labels
p <- ggplot(dt) +
  geom_bar(aes(x = -log10(pvalue), y = Description, fill = Cluster),
           stat = 'identity', width = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text.y = element_blank()) +
  geom_text(aes(x = 0.1, y = Description, label = Description),
            size = 4.5, hjust = 0) +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol)

# Save plot
ggsave("LZY/spRNA/上皮异质性/KEGG.celltype.barplot.pdf", p, width = 5, height = 6)


### GO enrichment analysis for epithelial subtypes
library(clusterProfiler)
library(ggplot2)
library(dplyr)

# Run GO (Biological Process) enrichment by cluster
GO_result <- compareCluster(
  ENTREZID ~ cluster,
  data = markers,
  fun = 'enrichGO',
  OrgDb = 'org.Hs.eg.db',
  keyType = 'ENTREZID',
  ont = 'BP'
)

# Save full GO result
write.csv(GO_result@compareClusterResult, "LZY/spRNA/上皮异质性/spRNAsub.GOBP.markers.clusters.csv")

# Basic dotplot
p_dot <- dotplot(GO_result) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
print(p_dot)

# Load filtered and selected GO terms for barplot
dt <- read.csv('LZY/spRNA/上皮异质性/select_GO.csv', header = TRUE)

# Set plotting factor levels
GO_levels <- c("CITED1+DTC","S100A1+DTC","DNAJA4+DTC","NPW+DTC",
               "BMP8A+DHGTC","SIGLEC6+DHGTC","SLC34A21+DHGTC",
               "SERPINE1+ATC","KRT5+ATC")

dt$Description <- factor(dt$Description, levels = rev(dt$Description))
dt$cluster <- factor(dt$cluster, levels = GO_levels)

# Create barplot
p <- ggplot(dt, aes(x = -log10(pvalue), y = Description, fill = cluster)) +
  geom_bar(stat = 'identity', width = 0.5) +
  theme_classic()

# Adjust x-axis origin
p1 <- p + scale_x_continuous(expand = c(0, 0))

# Add GO labels
p2 <- p1 +
  theme(axis.text.y = element_blank()) +
  geom_text(aes(x = 0.1, y = Description, label = Description), size = 4.5, hjust = 0)

# Define manual colors
mycol <- c("#4b6aa8", "#ece399", "#2d3462", "#c376a7",
           "#408444", "#61bada", "#b05545", "#df5734", "#9d3b62")

# Final colored barplot
p3 <- p2 +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol)

# Save final plot
ggsave("LZY/spRNA/上皮异质性/GOBP.celltype.barplot.pdf", p3, width = 5, height = 6)



############################################################
## GSVA-Based Metabolic Pathway Enrichment (Celltype-Level)
############################################################

# Load required libraries
library(Seurat)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(limma)
library(ggpubr)
options(stringsAsFactors = FALSE)

# Step 1: Load custom KEGG metabolism gene sets in GMT format
genesets <- read.gmt("LZY/spRNA/上皮异质性/KEGG_metabolism_nc.gmt")
genesets <- split(genesets$gene, genesets$term)  # Convert to list

# Step 2: Extract average expression matrix per group
Idents(spRNAsub) <- "niche0.3"
expr <- AverageExpression(spRNAsub, assays = "SCT", slot = "data")[[1]]
expr <- expr[rowSums(expr) > 0, ]  # Filter out zero-sum genes
expr <- as.matrix(expr)

# Step 3: Perform GSVA using Gaussian kernel
gsva.res <- gsva(expr, genesets, method = "gsva", kcdf = "Gaussian")

# Optional: Clean gene set names if needed
rownames(gsva.res) <- str_replace(rownames(gsva.res), "HALLMARK_", "")

# Step 4: Generate GSVA heatmap (not scaled)
p <- pheatmap::pheatmap(
  gsva.res,
  show_colnames = TRUE,
  scale = "none",
  angle_col = 45,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99)
)

# Save the heatmap to PDF
ggsave("LZY/spRNA/上皮异质性/GSVA.celltype.metabiolism.enrich.pheatmap.pdf", p, width = 8, height = 12)


# Load required packages
library(GeneNMF)
library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(dplyr)
library(tibble)

# Load Seurat object and subset target niches
spRNA <- readRDS("after_nicheidentified_spRNA.rds")
target_niches <- c("CITED1+DTC","S100A1+DTC","DNAJA4+DTC","NPW+DTC","BMP8A+DHGTC",
                   "SIGLEC6+DHGTC","SLC34A21+DHGTC","SERPINE1+ATC","KRT5+ATC")
spRNAsub <- subset(spRNA, niche0.3 %in% target_niches)

# Split object by sample and run multi-NMF
seu.list <- SplitObject(spRNAsub, split.by = "orig.ident")
geneNMF.programs <- multiNMF(seu.list, assay="SCT", k=4:10, min.exp=0.05, nfeatures=2000)

# Merge and summarize gene programs
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs, metric="cosine", weight.explained=0.5,
                                        nMP=9, max.genes=200)

# Save metaprogram object
save(geneNMF.metaprograms, file="spRNA.geneNMF.metaprograms.RData")

# Plot similarity heatmap
pdf("geneNMF.metaprograms.heatmap.pdf", width=9, height=8)
plotMetaPrograms(geneNMF.metaprograms, similarity.cutoff=c(0.1,1),
                 palette = c('white','#bfd3e6','#9ebcda','#8c96c6',
                             '#8c6bb1','#88419d','#810f7c','#4d004b'))
dev.off()

# Extract gene list from metaprograms
MP_genes <- unlist(geneNMF.metaprograms$metaprograms.genes)
MP_labels <- rep(names(geneNMF.metaprograms$metaprograms.genes), lengths(geneNMF.metaprograms$metaprograms.genes))
MP_genes_df <- data.frame(gene = MP_genes, MP = MP_labels)

# GO enrichment using clusterProfiler
gene_id_map <- bitr(MP_genes_df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
merged_genes <- merge(gene_id_map, MP_genes_df, by.x="SYMBOL", by.y="gene")

go_results <- compareCluster(ENTREZID ~ MP, data=merged_genes, fun="enrichGO",
                             OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=0.05)

# Simplify and visualize
go_simplified <- simplify(go_results, cutoff=0.7, by="p.adjust", select_fun=min)
ggplot(go_simplified, aes(Cluster, Description)) +
  geom_point(aes(size=GeneRatio, color=qvalue)) +
  scale_color_gradientn(colors=c('#CC3333','#FFFF99','#6699CC')) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("GOBP.geneNMF.metaprograms.pdf", width=9, height=8)

# Save results
write.csv(go_simplified@compareClusterResult, "GOBP.geneNMF.metaprograms.csv", row.names=TRUE)

# GSEA with MSigDB category C5:GO:BP
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(prog_genes) {
  runGSEA(prog_genes, universe=rownames(spRNAsub), category="C5", subcategory="GO:BP")
})
save(top_p, file="GOBP.geneNMF.metaprograms.top_p.RData")

# Module scoring per metaprogram
mp_genes <- geneNMF.metaprograms$metaprograms.genes
spRNAsub <- AddModuleScore(spRNAsub, features=mp_genes, name="", ncores=4)
VlnPlot(spRNAsub, features=1:9, group.by="niche0.3", pt.size=0, ncol=3) +
  scale_fill_manual(values=c("#4b6aa8","#ece399","#2d3462","#c376a7","#408444",
                             "#61bada","#b05545","#df5734","#9d3b62"))
ggsave("VlnPlot.geneNMF.metaprograms.pdf", width=9, height=8)
