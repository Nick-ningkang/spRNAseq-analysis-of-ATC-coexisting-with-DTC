######################################################################
###############Multi-sample integration analysis######################
######################################################################
#### Loading R packages ####
library(Seurat)
library(hdf5r)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggsci)
library(stringr)
library(RColorBrewer)
library(harmony)
library(forcats)
#### Function to visualize spatial clustering without visible labels
plot_spatial_dim <- function(spRNA, celltype_colors,group ,label = TRUE,label_size = 2) {
  p1 <- SpatialDimPlot(spRNA, label = label, label.size = label_size, images = "P12",pt.size.factor = 3.5,group.by = group) + 
    theme(legend.position = "none") + scale_fill_manual(values = celltype_colors) + ggtitle("P12")
  p2 <- SpatialDimPlot(spRNA, label = label, label.size = label_size, images = "P17",pt.size.factor = 4,group.by = group) + 
    theme(legend.position = "none") + scale_fill_manual(values = celltype_colors) + ggtitle("P17")
  p3 <- SpatialDimPlot(spRNA, label = label, label.size = label_size, images = "P26",pt.size.factor = 2,group.by = group) + 
    theme(legend.position = "none") + scale_fill_manual(values = celltype_colors) + ggtitle("P26")
  p4 <- SpatialDimPlot(spRNA, label = label, label.size = label_size, images = "P32",pt.size.factor = 4,group.by = group) + 
    theme(legend.position = "none") + scale_fill_manual(values = celltype_colors) + ggtitle("P32")
  p5 <- SpatialDimPlot(spRNA, label = label, label.size = label_size, images = "P44",pt.size.factor = 4,group.by = group) + 
    theme(legend.position = "none") + scale_fill_manual(values = celltype_colors) + ggtitle("P44")
  p7 <- SpatialDimPlot(spRNA, label = label, label.size = label_size, images = "P83",pt.size.factor = 2,group.by = group) + 
    theme(legend.position = "none") + scale_fill_manual(values = celltype_colors) + ggtitle("P83")
  p8 <- SpatialDimPlot(spRNA, label = label, label.size = label_size, images = "P98",pt.size.factor = 3.5,group.by = group) + 
    theme(legend.position = "none") + scale_fill_manual(values = celltype_colors) + ggtitle("P98")
  combined_plot <- (p1 | p2 | p3 | p4) / (p5 | p7 | p8)
  return(combined_plot)
}

##################################################
########## Harmony Integration Analysis ##########
##################################################
# Remove all objects from the current R environment and perform garbage collection
rm(list = ls())
gc()

# Load preprocessed spatial transcriptomics data (Seurat object containing multiple samples)
spRNA <- readRDS("data/spRNA_single_sample.rds")

# Check the number of spots per sample
table(spRNA$sample)

# Normalize data using SCTransform on the Spatial assay (recommended for multi-sample integration)
spRNA <- SCTransform(spRNA, assay = "Spatial", verbose = FALSE)

# Perform Principal Component Analysis on SCT-normalized data
spRNA <- RunPCA(spRNA, assay = "SCT", verbose = FALSE)

# Apply Harmony integration to correct for sample-specific batch effects
spRNA <- RunHarmony(spRNA, group.by.vars = "sample")

# Construct the nearest-neighbor graph using Harmony-reduced dimensions
spRNA <- FindNeighbors(spRNA, reduction = "harmony", dims = 1:10)

# Identify cell clusters using Louvain algorithm at multiple resolutions
spRNA <- FindClusters(spRNA, resolution = c(0.1, 0.2, 0.3, 0.5))

# Run UMAP for visualization using Harmony-integrated dimensions
spRNA <- RunUMAP(spRNA, reduction = "harmony", dims = 1:10)

# Extract image names (sample identifiers)
samples_name <- names(spRNA@images)

# Create output directory for integration results
dir.create("integration/harmony")


# Define marker genes for major cell types or functional states
cgs = list( 
  Differentiation = c("TG", "PAX8", "TSHR", "NKX2-1"),        # Thyroid differentiation markers
  Epithelial = c("CA9", "EPCAM", "KRT7", "KRT19", "CLDN4"),    # Epithelial cell markers
  SCC = c("KRT5", "KRT14", "TP63", "SOX2", "S100A9", "EGFR"),  # Squamous cell carcinoma markers
  Meyloid = c("CD68", "CD163", "LYZ", "AIF1", "MS4A6A"),       # Myeloid/macrophage markers
  T_cell = c("CD3D", "CD3E", "CD2", "CD8A"),                   # T cell markers
  NK = c("TYROBP", "XCL1"),                                   # NK cell markers
  B_cell = c("CD19", "MZB1", "MS4A1", "CD79A"),                # B cell markers
  TLS = c("CXCL13", "LAMP3", "CCL19", "CCL21"),                # Tertiary lymphoid structure markers
  Fibro = c("COL1A2", "COL3A1", "ACTA2", "PDPN", "COL1A1"),    # Fibroblast markers
  Endo = c("VWF", "FLT1", "RAMP2", "MYH11", "PECAM1")          # Endothelial cell markers
)

# Define color palette for visualizing different clusters
celltype_colors <- c(
  "#4b6aa8", "#ece399", "#df5734", "#c376a7", "#408444",
  "#b7deea", "#61bada", "#2d3462", "#9d3b62", "#83ab8e",
  "#d69a55", "#696a6c", "#d25774", "#cea5c7", "#8d689d",
  "#6c408e", "#efd2c9", "#a78982", "#e6b884", "#4490c4"
)

# Set clustering identity based on resolution 0.3
Idents(spRNA) <- "SCT_snn_res.0.3"

# Spatial plot of clusters on tissue sections
combined_plot <- plot_spatial_dim(spRNA, celltype_colors, "SCT_snn_res.0.3")
print(combined_plot)

# UMAP plot colored by cluster identity
p_cluster <- DimPlot(spRNA, reduction = "umap", label = TRUE, cols = celltype_colors)

# UMAP plot colored by sample origin
p_sample <- DimPlot(spRNA, reduction = "umap", label = TRUE, cols = celltype_colors, group.by = "sample")

# Dot plot showing marker gene expression across clusters
p_markgene <- DotPlot(spRNA, features = cgs, assay = "SCT") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1))

# Identify cluster-specific marker genes
cluster_markers <- FindAllMarkers(spRNA, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

# Save generated plots to PDF files
ggsave("integration/harmony/Clusters_image0.3.pdf", combined_plot, width = 10, height = 8)
ggsave("integration/harmony/Clusters_umap0.3.pdf", p_cluster, width = 8, height = 6)
ggsave("integration/harmony/sample_umap0.3.pdf", p_sample, width = 8, height = 6)
ggsave("integration/harmony/markgene0.3.pdf", p_markgene, width = 14, height = 6)

# Export identified cluster markers to CSV file
write.csv(cluster_markers, "integration/harmony/cluster_markers0.3.csv")

# Save the processed Seurat object after Harmony integration
saveRDS(spRNA, "data/after_harmony_spRNA.rds")


####################################################
########## Annotation of different niches ##########
####################################################
Idents(spRNA) <- "SCT_snn_res.0.3"

# Initialize a new metadata column for niche annotation
spRNA$niche0.3 <- ""

# Assign niche identities based on cluster ID
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="0"] <- "CITED1+DTC"
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="1"] <- 'SLC34A21+HGFCTC'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="2"] <- 'CCL9+Immune_cells'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="3"] <- 'SERPINE1+ATC'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="4"] <- 'SIGLEC6+HGFCTC'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="5"] <- 'KRT5+ATC'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="6"] <- 'SFRP4+Fibroblast'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="7"] <- 'DES+Myocyte'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="8"] <- 'BMP8A+HGFCTC'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="9"] <- 'FOSB+Immune_cells'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="10"] <- 'S100A1+DTC'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="11"] <- 'DNAJA4+DTC'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="12"] <- 'TTN+Myocyte'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="13"] <- 'MARCO+Immune_cells'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="14"] <- 'NPW+DTC'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="15"] <- 'KRT4+NOR'
spRNA$niche0.3[spRNA$SCT_snn_res.0.3=="16"] <- 'SLC34A21+HGFCTC'

# For regions manually annotated as "Squamous epithelium", overwrite with KRT4+NOR
spRNA$niche0.3[spRNA$region_single=="Squamous epithelium"] <- 'KRT4+NOR'

# Save the annotated Seurat object after niche classification
saveRDS(spRNA,"data/after_nicheidentified_spRNA.rds")


#### UMAP visualization of spatial niches
# Set the factor levels for consistent ordering
spRNA$niche0.3 <- factor(spRNA$niche0.3, levels = c(
  "CITED1+DTC", "S100A1+DTC", "DNAJA4+DTC", "NPW+DTC",      # DTC
  "BMP8A+HGFCTC", "SIGLEC6+HGFCTC", "SLC34A21+HGFCTC",         # HGFCTC
  "SERPINE1+ATC", "KRT5+ATC",                               # ATC
  "CCL9+Immune_cells", "FOSB+Immune_cells", "MARCO+Immune_cells",
  "SFRP4+Fibroblast", "DES+Myocyte", "TTN+Myocyte",         # Stromal
  "KRT4+NOR"                                                # Normal
))

# Assign the identity class
Idents(spRNA) <- "niche0.3"

# Define colors for each annotated spatial niche
celltype_colors <- c(
  'CITED1+DTC' = "#4b6aa8", 'S100A1+DTC' = "#ece399", 'DNAJA4+DTC' = "#2d3462", 'NPW+DTC' = "#c376a7",
  'BMP8A+HGFCTC' = "#408444", 'SIGLEC6+HGFCTC' = "#61bada", 'SLC34A21+HGFCTC' = '#b05545',
  'SERPINE1+ATC' = "#df5734", 'KRT5+ATC' = "#9d3b62",
  'CCL9+Immune_cells' = "#83ab8e", 'FOSB+Immune_cells' = "#d69a55", 'MARCO+Immune_cells' = "#d25774",
  'SFRP4+Fibroblast' = "#696a6c", 'DES+Myocyte' = "#a78982", 'TTN+Myocyte' = "#efd2c9",
  'KRT4+NOR' = "#6c408e"
)

# Generate spatial plots with and without labels
combined_plot <- plot_spatial_dim(spRNA, celltype_colors, "niche0.3"); print(combined_plot)
combined_plot_withoutlabel <- plot_spatial_dim(spRNA, label = FALSE, celltype_colors, "niche0.3"); print(combined_plot_withoutlabel)

# Save spatial plots to PDF
ggsave("Figure/niche_image0.3.pdf", combined_plot, width = 10, height = 8)
ggsave("Figure/niche_image0.3_withoutlabel.pdf", combined_plot_withoutlabel, width = 10, height = 8)

# UMAP visualization of annotated niches and samples
p_celltype <- DimPlot(spRNA, reduction = "umap", label = TRUE, cols = celltype_colors)
p_sample <- DimPlot(spRNA, reduction = "umap", label = TRUE, cols = celltype_colors, group.by = "sample")
ggsave("Figure/UMAP.pdf", p_celltype | p_sample, width = 12, height = 5)


#### DotPlot visualization of representative marker genes
# Define cell-type-specific marker genes
cgs = list( 
  Differentiation = c("TG", "PAX8", "TSHR", "NKX2-1"),
  Epithelial = c("EPCAM", "KRT7", "CLDN4"),
  Squamous = c("KRT5", "TP63"),
  T_cell = c("CD3D", "CD3E", "CD2"),
  B_cell = c("MZB1", "MS4A1", "CD79A", "AKNA", "CIITA"),
  Meyloid = c("CD68", "CD163", "LYZ"),
  Fibro = c("COL1A2", "COL3A1", "COL1A1"),
  Myocyte = c("TNC", "MYH11", "DES")
)

# Create DotPlot using SCT assay
p_markgene <- DotPlot(spRNA, features = cgs, assay = "SCT") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1))

# Extract data from DotPlot object for custom plotting
exp <- p_markgene$data
exp$features.plot <- as.factor(exp$features.plot)
exp$features.plot <- fct_inorder(exp$features.plot)

# Customize DotPlot appearance using ggplot2
p <- ggplot(exp, aes(x = id, y = features.plot)) +
  geom_point(aes(size = `pct.exp`, color = `avg.exp.scaled`)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('#6699CC', '#FFFF99', '#CC3333')
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  # Add dotted vertical lines to separate niche groups
  geom_vline(xintercept = c(4.5, 7.5, 9.5, 10.5, 13.5, 14.5), linetype = "dotted", size = 1) +
  # Highlight niche group regions using transparent rectangles
  geom_rect(aes(xmin = 1 - 0.4, xmax = 10 + 0.4, ymin = 1 - 0.4, ymax = 9 + 0.4), fill = "transparent", color = "black", size = 0.5) +
  geom_rect(aes(xmin = 11 - 0.4, xmax = 13 + 0.4, ymin = 10 - 0.4, ymax = 20 + 0.4), fill = "transparent", color = "black", size = 0.5) +
  geom_rect(aes(xmin = 14 - 0.4, xmax = 14 + 0.4, ymin = 21 - 0.4, ymax = 23 + 0.4), fill = "transparent", color = "black", size = 0.5) +
  geom_rect(aes(xmin = 15 - 0.4, xmax = 16 + 0.4, ymin = 24 - 0.4, ymax = 26 + 0.4), fill = "transparent", color = "black", size = 0.5) +
  theme(legend.direction = "horizontal", legend.position = "bottom")

# Save the final customized marker gene DotPlot
ggsave("Figure/markergene.pdf", p, width = 6, height = 8)

# Save the processed Seurat object after Annotation of different niches
saveRDS(spRNA, "data/after_nicheidentified_spRNA")



####################################################
############ Plotting and visualization ############
####################################################
#===========================================
# Gene Expression Visualization
# Define a minimal theme to remove axes, ticks, borders, and legend for a cleaner look
theme <- theme(
  axis.title = element_blank(),       # Remove axis titles
  axis.text = element_blank(),        # Remove axis text
  axis.ticks = element_blank(),       # Remove axis ticks
  axis.line = element_blank(),        # Remove axis lines
  panel.border = element_blank(),     # Remove plot borders
  panel.grid = element_blank(),       # Remove grid lines
  legend.position = "none"            # Remove legend
)

# Plot spatial expression of selected genes using FeaturePlot (UMAP-based)
p1 <- FeaturePlot(object = spRNA, features = c("TG"), cols = c('#6699CC', '#FFFF99', '#CC3333'), order = TRUE) + theme
p2 <- FeaturePlot(object = spRNA, features = c("TSHR"), cols = c('#6699CC', '#FFFF99', '#CC3333'), order = TRUE) + theme
p3 <- FeaturePlot(object = spRNA, features = c("PAX8"), cols = c('#6699CC', '#FFFF99', '#CC3333'), order = TRUE) + theme
p4 <- FeaturePlot(object = spRNA, features = c("NKX2-1"), cols = c('#6699CC', '#FFFF99', '#CC3333'), order = TRUE) + theme

# Combine and display the 4 gene plots in a single row
p1 | p2 | p3 | p4

# Save the combined FeaturePlot as a PDF
ggsave("Figure/分化基因可视化.pdf", p1 | p2 | p3 | p4, width = 12, height = 3)

# Create a violin plot showing expression of selected genes across niches
p <- VlnPlot(
  spRNA,
  features = c("TG", "TSHR", "PAX8", "NKX2-1"),
  pt.size = 0,
  stack = TRUE,
  flip = TRUE,
  cols = celltype_colors,
  split.by = "niche0.3"
) + theme(legend.position = "none")

# Save the violin plot
ggsave("Figure/分化基因小提琴图.pdf", p, width = 9, height = 4)


# Define a custom function for spatial feature plotting
# across multiple samples (images)
plot_spatial_feature <- function(spRNA, features) {
  # Create individual spatial plots for each tissue section
  p1 <- SpatialFeaturePlot(spRNA, features = features, images = "P12", pt.size.factor = 3.5) +
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) +
    ggtitle("P12")
  
  p2 <- SpatialFeaturePlot(spRNA, features = features, images = "P17", pt.size.factor = 4) +
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) +
    ggtitle("P17")
  
  p3 <- SpatialFeaturePlot(spRNA, features = features, images = "P26", pt.size.factor = 2) +
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) +
    ggtitle("P26")
  
  p4 <- SpatialFeaturePlot(spRNA, features = features, images = "P32", pt.size.factor = 4) +
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) +
    ggtitle("P32")
  
  p5 <- SpatialFeaturePlot(spRNA, features = features, images = "P44", pt.size.factor = 4) +
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) +
    ggtitle("P44")
  
  p6 <- SpatialFeaturePlot(spRNA, features = features, images = "P57", pt.size.factor = 2) +
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) +
    ggtitle("P57")
  
  p7 <- SpatialFeaturePlot(spRNA, features = features, images = "P83", pt.size.factor = 2) +
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) +
    ggtitle("P83")
  
  p8 <- SpatialFeaturePlot(spRNA, features = features, images = "P98", pt.size.factor = 3.5) +
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) +
    ggtitle("P98")
  
  # Arrange all plots in a 2-row layout
  combined_plot <- (p1 | p2 | p3 | p4) / (p5 | p6 | p7 | p8)
  return(combined_plot)
}


# Spatial expression plots of differentiation-related genes

p <- plot_spatial_feature(spRNA, features = "TG")
ggsave("Figure/TG_spatial.pdf", p, width = 10, height = 6)

p <- plot_spatial_feature(spRNA, features = "TSHR")
ggsave("Figure/TSHR_spatial.pdf", p, width = 10, height = 6)

p <- plot_spatial_feature(spRNA, features = "PAX8")
ggsave("Figure/PAX8_spatial.pdf", p, width = 10, height = 6)

p <- plot_spatial_feature(spRNA, features = "NKX2-1")
ggsave("Figure/NKX2-1_spatial.pdf", p, width = 10, height = 6)



#===========================================
# Cell composition analysis across samples
# 1. Analyze niche composition in ATC regions
# Subset data to include only cells from ATC regions
spRNAsub <- subset(spRNA, region_single == "ATC")

# Tabulate the number of cells from each sample (Var1) and each niche (Var2)
bar_data <- table(spRNAsub$sample, spRNAsub$niche0.3) %>% as.data.frame()

# Calculate relative proportions (%) within each sample
bar_per <- bar_data %>%
  group_by(Var1) %>%
  mutate(total = sum(Freq)) %>%
  mutate(percent = Freq / total)

# Optional: View specific niches of interest (e.g., ATC-related)
subset(bar_per, Var2 %in% c('SERPINE1+ATC', 'KRT5+ATC'))

# Set consistent factor order for bar plot
bar_per$Var2 <- factor(bar_per$Var2, levels = c(
  "CITED1+DTC", "S100A1+DTC", "DNAJA4+DTC", "NPW+DTC",             # DTC
  "BMP8A+HGFCTC", "SIGLEC6+HGFCTC", "SLC34A21+HGFCTC",                # HGFCTC
  "SERPINE1+ATC", "KRT5+ATC",                                      # ATC
  "CCL9+Immune_cells", "FOSB+Immune_cells", "MARCO+Immune_cells",  # Immune
  "SFRP4+Fibroblast", "DES+Myocyte", "TTN+Myocyte", "KRT4+NOR"     # Stromal & Normal
))

# Create a horizontal stacked bar plot of niche proportions by sample
p <- ggplot(bar_per, aes(x = Var1, y = percent)) +
  geom_bar(aes(fill = Var2), stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = celltype_colors) +
  theme(
    axis.ticks = element_line(linetype = "blank"),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA),
    plot.background = element_rect(colour = NA)
  ) +
  labs(x = NULL, y = "% Relative cell source", fill = NULL)

# Save the plot for ATC samples
ggsave("Figure/niche_percent1.pdf", p, width = 4, height = 4)


# 2. Analyze niche composition in HGFCTC and PTC regions

# Subset data to include only HGFCTC and PTC regions
spRNAsub <- subset(spRNA, region_single %in% c("HGFCTC", "PTC"))

# Tabulate and calculate proportions as above
bar_data <- table(spRNAsub$sample, spRNAsub$niche0.3) %>% as.data.frame()

bar_per <- bar_data %>%
  group_by(Var1) %>%
  mutate(total = sum(Freq)) %>%
  mutate(percent = Freq / total)

# Optional: Inspect selected DTC and HGFCTC niches
subset(bar_per, Var2 %in% c("CITED1+DTC", "S100A1+DTC", "DNAJA4+DTC", "NPW+DTC"))
subset(bar_per, Var2 %in% c("BMP8A+HGFCTC", "SIGLEC6+HGFCTC", "SLC34A21+HGFCTC"))

# Set factor level order for consistency
bar_per$Var2 <- factor(bar_per$Var2, levels = c(
  "CITED1+DTC", "S100A1+DTC", "DNAJA4+DTC", "NPW+DTC",             # DTC
  "BMP8A+HGFCTC", "SIGLEC6+HGFCTC", "SLC34A21+HGFCTC",                # HGFCTC
  "SERPINE1+ATC", "KRT5+ATC",                                      # ATC
  "CCL9+Immune_cells", "FOSB+Immune_cells", "MARCO+Immune_cells",  # Immune
  "SFRP4+Fibroblast", "DES+Myocyte", "TTN+Myocyte", "KRT4+NOR"     # Stromal & Normal
))

# Generate and save bar plot for HGFCTC/PTC regions
p <- ggplot(bar_per, aes(x = Var1, y = percent)) +
  geom_bar(aes(fill = Var2), stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = celltype_colors) +
  theme(
    axis.ticks = element_line(linetype = "blank"),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA),
    plot.background = element_rect(colour = NA)
  ) +
  labs(x = NULL, y = "% Relative cell source", fill = NULL)

ggsave("Figure/niche_percent2.pdf", p, width = 10, height = 6)









#===========================================
# Ro/e Analysis: Tissue Enrichment Visualization

# 1. Ro/e plot for all samples
# Extract metadata from Seurat object
metainfo <- spRNA@meta.data

# Load the custom function for Ro/e distribution analysis
source("/datadisk1/person/LZY/Gender/Code/distribution_Roe.R")

# Define custom colors for each spatial niche
cellcolors <- c(
  'CITED1+DTC' = "#4b6aa8", 'S100A1+DTC' = "#ece399", 'DNAJA4+DTC' = "#2d3462", 'NPW+DTC' = "#c376a7",
  'BMP8A+HGFCTC' = "#408444", 'SIGLEC6+HGFCTC' = "#61bada", 'SLC34A21+HGFCTC' = '#b05545',
  'SERPINE1+ATC' = "#df5734", 'KRT5+ATC' = "#9d3b62",
  'CCL9+Immune_cells' = "#83ab8e", 'FOSB+Immune_cells' = "#d69a55", 'MARCO+Immune_cells' = "#d25774",
  'SFRP4+Fibroblast' = "#696a6c", 'DES+Myocyte' = "#a78982", 'TTN+Myocyte' = "#efd2c9", 'KRT4+NOR' = "#6c408e"
)

# Run Ro/e (Ratio of observed vs expected) analysis across samples
distribution_Roe(
  meta_data = metainfo,                        # Metadata with cell annotations
  celltype_column = "niche0.3",                # Column indicating niche identity
  condition_column = "sample",                 # Grouping condition (e.g., by sample)
  celltype_level = c(                          # Predefined order of niche types
    "CITED1+DTC", "S100A1+DTC", "DNAJA4+DTC", "NPW+DTC",      # DTC
    "BMP8A+HGFCTC", "SIGLEC6+HGFCTC", "SLC34A21+HGFCTC",         # HGFCTC
    "SERPINE1+ATC", "KRT5+ATC",                                # ATC
    "CCL9+Immune_cells", "FOSB+Immune_cells", "MARCO+Immune_cells",
    "SFRP4+Fibroblast", "DES+Myocyte", "TTN+Myocyte", "KRT4+NOR"
  ),
  add_label = "sign",                          # Show significance markers (e.g., asterisks)
  celltype_color = cellcolors,                 # Use consistent colors
  relative_width = 0.2,                        # Relative width for significance bars
  tile_color = NA                              # Default background color
)

# Save the plot
ggsave("Figure/组织偏好性Roe图.pdf", width = 3, height = 6)


# 2. Ro/e plot for ATC samples only

# Subset metadata to include only ATC region cells
metainfosub <- subset(metainfo, metainfo$region_single == "ATC")

# Reload the Ro/e function (optional redundancy if environment is refreshed)
source("/datadisk1/person/LZY/Gender/Code/distribution_Roe.R")

# Run Ro/e analysis for ATC region only
distribution_Roe(
  meta_data = metainfosub,
  celltype_column = "niche0.3",
  condition_column = "sample",
  celltype_level = c(
    "CITED1+DTC", "S100A1+DTC", "DNAJA4+DTC", "NPW+DTC",
    "BMP8A+HGFCTC", "SIGLEC6+HGFCTC", "SLC34A21+HGFCTC",
    "SERPINE1+ATC", "KRT5+ATC",
    "CCL9+Immune_cells", "FOSB+Immune_cells", "MARCO+Immune_cells",
    "SFRP4+Fibroblast", "DES+Myocyte", "TTN+Myocyte", "KRT4+NOR"
  ),
  add_label = "sign",
  celltype_color = cellcolors,
  relative_width = 0.2,
  tile_color = NA
)











