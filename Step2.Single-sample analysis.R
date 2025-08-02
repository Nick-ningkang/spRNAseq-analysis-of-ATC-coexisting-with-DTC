######################################################################
######################Single-sample analysis##########################
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
library(ggrepel)
library(ggVennDiagram)
library(UpSetR)
###########################################
#### Single-sample Clustering Analysis ####
###########################################
# Load Spatial Transcriptomics Data 
rm(list = ls())         # Clear environment
gc()                    # Garbage collection
setwd('/datadisk1/person/NICK/spRNA')  # Set working directory
spRNA <- readRDS("data/after_QC_spRNA.rds")  # QC-passed Seurat object

sample_names <- unique(spRNA$sample)          # Extract sample identifiers
dir.create("pathological_segregation/cluster", recursive = TRUE)  # Create output directory
spRNAlist <- list()                           # Initialize result list
getPalette <- colorRampPalette(brewer.pal(8, "Set1"))  # Define color palette

#Marker genes serving as reference for downstream pathological annotation
cgs = list( 
  Differentiation = c("TG","PAX8",'TSHR','NKX2-1'),
  Epithelial = c("CA9","EPCAM","KRT7","KRT19","CLDN4"),
  SCC =c("KRT5","KRT14", "TP63", "SOX2", "S100A9","EGFR"),
  Meyloid = c("CD68","CD163","LYZ","AIF1","MS4A6A"), 
  T_cell = c("CD3D",'CD3E','CD2',"CD8A"),  
  NK=c("TYROBP","XCL1"),
  B_cell = c("CD19","MZB1","MS4A1","CD79A"), 
  TLS=c('CXCL13','LAMP3','CCL19','CCL21'),
  Fibro = c("COL1A2","COL3A1",'ACTA2','PDPN','COL1A1'),
  Endo = c("VWF","FLT1","RAMP2",'MYH11','PECAM1')
)

for (i in seq_along(sample_names)) {
  
  # Subset Seurat object by sample
  sample_id <- sample_names[i]
  spRNA_sub <- subset(spRNA, sample == sample_id)
  
  # Normalize & scale using SCTransform
  spRNA_sub <- SCTransform(spRNA_sub, assay = "Spatial", verbose = FALSE)
  
  # Dimensionality reduction & clustering
  spRNA_sub <- RunPCA(spRNA_sub, assay = "SCT", verbose = FALSE)
  spRNA_sub <- FindNeighbors(spRNA_sub, reduction = "pca", dims = 1:10)
  spRNA_sub <- FindClusters(spRNA_sub, resolution = 0.8, verbose = FALSE)
  spRNA_sub <- RunUMAP(spRNA_sub, reduction = "pca", dims = 1:10)
  
  # Define color palette based on number of clusters
  cluster_colors <- getPalette(length(unique(spRNA_sub$SCT_snn_res.0.8)) + 1)
  
  # Generate plots
  p1 <- DimPlot(spRNA_sub, reduction = "umap", label = TRUE, cols = cluster_colors)
  p2 <- SpatialDimPlot(
    spRNA_sub, label = TRUE, label.size = 3, 
    images = sample_id, pt.size.factor = 2
  ) +
    scale_fill_manual(values = cluster_colors) +
    theme(legend.position = "none")
  
  # DotPlot using pre-defined marker list 'cgs'
  p3 <- DotPlot(spRNA_sub, features = cgs, assay = "SCT") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1))
  
  # Identify marker genes
  markers <- FindAllMarkers(
    spRNA_sub, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25
  )
  
  # Save plots and marker genes
  ggsave(
    filename = paste0("pathological_segregation/cluster/", sample_id, "_UMAP.pdf"),
    plot = (p1 | p2) / p3, width = 15, height = 12
  )
  write.csv(markers, paste0("pathological_segregation/cluster/", sample_id, "_clustergene.csv"))
  
  # Store result
  spRNAlist[[i]] <- spRNA_sub
}

# Save all single-sample Seurat objects
saveRDS(spRNAlist, "data/spRNAlist.rds")


###########################################################################################
#### Manual annotation of pathological regions based on expert pathological assessment ####
###########################################################################################
rm(list = ls())         # Clear environment
gc()  
spRNAlist <- readRDS("data/spRNAlist.rds")  # QC-passed Seurat object

#Annotate pathological regions for each sample individually
#P12
i=1
spRNAlist[[i]]$region <- ""
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(2,3,8,12))] <- "DTC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(4,9,10))] <- "ATC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(1,5,11,13,7))] <- "Stroma"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(0,6))] <- "HGFCTC"
spRNAlist[[i]]$region <- factor(spRNAlist[[i]]$region,levels = c("ATC","DTC","Stroma","HGFCTC"))
celltype_colors <- c("#E64B35FF","#00A087FF","#3C5488FF","#F39B7FFF")
p1 <- SpatialDimPlot(spRNAlist[[i]], label = F, label.size = 2,images = "P12",pt.size.factor = 4,group.by ="region" )+scale_fill_manual(values=celltype_colors)+ theme(legend.position = "none") + ggtitle("P12")
ggsave('pathological_segregation/pathology/p_P12.pdf', p1, width = 10, height = 5)

#P17
i=2
spRNAlist[[i]]$region <- ""
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(0,3,7,8,12))] <- "DTC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(1,4,5,9))] <- "ATC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(2,6,10,11,13,14))] <- "Stroma"
spRNAlist[[i]]$region <- factor(spRNAlist[[i]]$region,levels = c("ATC","DTC","Stroma"))
celltype_colors <- c("#E64B35FF","#00A087FF","#3C5488FF")
p2 <- SpatialDimPlot(spRNAlist[[i]], label = F, label.size = 2,images = "P17",pt.size.factor = 4,group.by ="region" )+scale_fill_manual(values=celltype_colors)+ theme(legend.position = "none") + ggtitle("P17")
ggsave('pathological_segregation/pathology/p_P17.pdf', p2, width = 10, height = 5)

#P26
i=3     
spRNAlist[[i]]$region <- ""
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(0,2,4,7))] <- "DTC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(1,3,5,6,8))] <- "ATC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(9,10))] <- "Stroma"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(11,12))] <- "Tertiary"
spRNAlist[[i]]$region <- factor(spRNAlist[[i]]$region,levels = c("ATC","DTC","Stroma","Tertiary"))
celltype_colors <- c("#E64B35FF","#00A087FF","#3C5488FF","#8491B4FF")
p3 <- SpatialDimPlot(spRNAlist[[i]], label = F, label.size = 2,images = "P26",pt.size.factor = 2,group.by ="region" )+scale_fill_manual(values=celltype_colors)+ theme(legend.position = "none") + ggtitle("P26")
ggsave('pathological_segregation/pathology/p_P26.pdf', p3, width = 10, height = 5)

#P32
i=4     
spRNAlist[[i]]$region <- ""
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(2,4,5))] <- "DTC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(0,3,11))] <- "ATC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(1,6,8,10))] <- "Stroma"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(7,9))] <- "NSE"
spRNAlist[[i]]$region <- factor(spRNAlist[[i]]$region,levels = c("ATC","DTC","Stroma","NSE"))
celltype_colors <- c("#E64B35FF","#00A087FF","#3C5488FF","#B09C85FF")
p4 <- SpatialDimPlot(spRNAlist[[i]], label = F, label.size = 2,images = "P32",pt.size.factor = 4,group.by ="region" )+scale_fill_manual(values=celltype_colors)+ theme(legend.position = "none") + ggtitle("P32")
ggsave('pathological_segregation/pathology/p_P32.pdf', p4, width = 10, height = 5)

#P44
i=5     
spRNAlist[[i]]$region <- ""
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(5))] <- "DTC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(4,11))] <- "ATC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(6,9,10,12))] <- "Stroma"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(0,1,2,3,7,8))] <- "HGFCTC"
spRNAlist[[i]]$region <- factor(spRNAlist[[i]]$region,levels = c("ATC","DTC","Stroma","HGFCTC"))
celltype_colors <- c("#E64B35FF","#00A087FF","#3C5488FF","#F39B7FFF")
p5 <- SpatialDimPlot(spRNAlist[[i]], label = F, label.size = 2,images = "P44",pt.size.factor = 4,group.by ="region" )+scale_fill_manual(values=celltype_colors)+ theme(legend.position = "none") + ggtitle("P44")
ggsave('pathological_segregation/pathology/p_P44.pdf', p5, width = 10, height = 5)

#P83
i=6     
spRNAlist[[i]]$region <- ""
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(1,10,13))] <- "DTC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(0,2,8,9,14))] <- "ATC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(12,5,6,11,15))] <- "Stroma"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(3,4,7))] <- "HGFCTC"
spRNAlist[[i]]$region <- factor(spRNAlist[[i]]$region,levels = c("ATC","DTC","Stroma","HGFCTC"))
celltype_colors <- c("#E64B35FF","#00A087FF","#3C5488FF","#F39B7FFF")
p6 <- SpatialDimPlot(spRNAlist[[i]], label = F, label.size = 2,images = "P83",pt.size.factor = 2,group.by ="region" )+scale_fill_manual(values=celltype_colors)+ theme(legend.position = "none") + ggtitle("P83")
ggsave('pathological_segregation/pathology/p_P83.pdf', p6, width = 10, height = 5)

#P98
i=7     
spRNAlist[[i]]$region <- ""
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(1,2,4,5,7,10))] <- "DTC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(0,6,12,13))] <- "HGFCTC"
spRNAlist[[i]]$region[which(spRNAlist[[i]]$seurat_clusters %in% c(3,8,9,11))] <- "Stroma"
spRNAlist[[i]]$region <- factor(spRNAlist[[i]]$region,levels = c("HGFCTC","DTC","Stroma"))
celltype_colors <- c("#F39B7FFF","#00A087FF","#3C5488FF")
p7 <- SpatialDimPlot(spRNAlist[[i]], label = F, label.size = 2,images = "P98",pt.size.factor = 3.5,group.by ="region" )+scale_fill_manual(values=celltype_colors)+ theme(legend.position = "none") + ggtitle("P98")
ggsave('pathological_segregation/pathology/p_P98.pdf', p7, width = 10, height = 5)

#Merge all figures and save the output
combined_plot <- (p1 | p2 | p3 | p4) / (p5 | p6 | p7 )
combined_plot
ggsave('Figure/pathological_segregation.pdf', combined_plot, width = 10, height = 8)

#Save single-sample spatial transcriptomics analysis
saveRDS(spRNAlist,"data/spRNAlist_single_sample.rds")

########################################################################################################
#### Differential Gene Expression Analysis (Highly Differentiated vs Poorly Differentiated Regions) ####
########################################################################################################
rm(list = ls())         # Clear environment
gc()  
spRNAlist <- readRDS("data/spRNAlist_single_sample.rds") #loading data

# Classify each sample into highly and poorly differentiated regions:
# If ATC is present, HGFCTC is assigned to the highly differentiated group;
# otherwise, HGFCTC is assigned to the poorly differentiated group.
for (i in 1:7) {spRNAlist[[i]]$region1 <- spRNAlist[[i]]$region}
spRNAlist[[1]]$region1[spRNAlist[[1]]$region1=="HGFCTC"] <- 'DTC'
spRNAlist[[5]]$region1[spRNAlist[[5]]$region1=="HGFCTC"] <- 'DTC'
spRNAlist[[6]]$region1[spRNAlist[[6]]$region1=="HGFCTC"] <- 'DTC'
spRNAlist[[7]]$region1[spRNAlist[[7]]$region1=="HGFCTC"] <- 'ATC'
for (i in 1:7) {spRNAlist[[i]]$region1[spRNAlist[[i]]$region1=="DTC"] <- "highly_differentiated"}
for (i in 1:7) {spRNAlist[[i]]$region1[spRNAlist[[i]]$region1=="ATC"] <- "poorly_differentiated"}

# Define sample IDs
samples <- c('P12', 'P17', 'P26', 'P32', 'P44', 'P83', 'P98')

# Loop over samples
for (i in seq_along(samples)) {
  
  # Perform differential gene expression analysis between poorly and highly differentiated regions
  deg <- FindMarkers(
    spRNAlist[[i]],
    ident.1 = "poorly_differentiated",
    ident.2 = "highly_differentiated",
    group.by = "region1",               
    logfc.threshold = 0.25,
    min.pct = 0.05
  )
  
  # Subset significant upregulated and downregulated genes
  up_sigdeg <- subset(deg, p_val_adj < 0.05 & avg_log2FC > 1)
  down_sigdeg <- subset(deg, p_val_adj < 0.05 & avg_log2FC < -1)
  
  # Dynamically assign results to global environment
  assign(paste0("deg_", samples[i]), deg)
  assign(paste0("up_sigdeg_", samples[i]), up_sigdeg)
  assign(paste0("down_sigdeg_", samples[i]), down_sigdeg)
  
  # Save results as CSV files
  write.csv(deg, file = paste0("pathological_segregation/DEG/deg_", samples[i], ".csv"))
  write.csv(up_sigdeg, file = paste0("pathological_segregation/DEG/up_sigdeg_", samples[i], ".csv"))
  write.csv(down_sigdeg, file = paste0("pathological_segregation/DEG/down_sigdeg_", samples[i], ".csv"))
}


#### Upregulated DEGs (highly_differentiated vs poorly_differentiated)
# Step 1: Find the maximum number of upregulated genes across samples
max_length <- max(sapply(list(up_sigdeg_P12, up_sigdeg_P17, up_sigdeg_P26,
                              up_sigdeg_P32, up_sigdeg_P44, up_sigdeg_P83, up_sigdeg_P98),
                         function(x) nrow(x)))

# Step 2: Normalize gene list length by padding with NAs
P12 <- rownames(up_sigdeg_P12); length(P12) <- max_length
P17 <- rownames(up_sigdeg_P17); length(P17) <- max_length
P26 <- rownames(up_sigdeg_P26); length(P26) <- max_length
P32 <- rownames(up_sigdeg_P32); length(P32) <- max_length
P44 <- rownames(up_sigdeg_P44); length(P44) <- max_length
P83 <- rownames(up_sigdeg_P83); length(P83) <- max_length
P98 <- rownames(up_sigdeg_P98); length(P98) <- max_length

# Step 3: Create a dataframe from upregulated genes and save to CSV
df_up <- data.frame(P12, P17, P26, P32, P44, P83, P98)
write.csv(df_up, "pathological_segregation/DEG/upgene.csv")

# Step 4: Draw an upset plot for upregulated genes
cors <- pal_npg("nrc")(7)
pdf("Figure/upset_poorly_vs_highly_upgene.pdf", width = 10, height = 6)
p <- upset(fromList(df_up), nsets = 8, nintersects = 35, order.by = "freq", keep.order = TRUE,
           sets.bar.color = cors, point.size = 3, line.size = 1,
           sets = c('P12','P17','P26','P32','P44','P83','P98'),
           main.bar.color = "#33476B", matrix.color = "#33476B", mb.ratio = c(0.6, 0.4),
           queries = list(
             list(query = intersects, params = list('P12','P17','P26','P32','P44','P83','P98'), color = "#AD1F1F", active = TRUE),
             list(query = intersects, params = list('P12','P17','P26','P32','P44','P83'), color = "#AD1F1F", active = TRUE)
           ))
print(p)
dev.off()

# Step 5: Identify commonly upregulated genes
common_up_all <- Reduce(intersect, df_up); print(common_up_all)
common_up_ATC <- Reduce(intersect, df_up[1:6]); print(common_up_ATC)#The first six samples all contain ATC.

# Step 6: Combine all upregulated gene metadata and export
samples <- c("P12", "P17", "P26", "P32", "P44", "P83", "P98")
for (sample_id in samples) {
  obj <- get(paste0("up_sigdeg_", sample_id))
  obj$sample <- sample_id
  obj$gene <- rownames(obj)
  assign(paste0("up_sigdeg_", sample_id), obj)
}
up_sigdeg_all <- do.call(rbind, mget(paste0("up_sigdeg_", samples)))
write.csv(up_sigdeg_all, "Figure/upgene_poorly_vs_highly.csv")



#### Down DEGs (highly_differentiated vs poorly_differentiated)

# Step 1: Find the maximum number of downregulated genes across samples
max_length <- max(sapply(list(down_sigdeg_P12, down_sigdeg_P17, down_sigdeg_P26,
                              down_sigdeg_P32, down_sigdeg_P44, down_sigdeg_P83, down_sigdeg_P98),
                         function(x) nrow(x)))

# Step 2: Normalize gene list length by padding with NAs
P12 <- rownames(down_sigdeg_P12); length(P12) <- max_length
P17 <- rownames(down_sigdeg_P17); length(P17) <- max_length
P26 <- rownames(down_sigdeg_P26); length(P26) <- max_length
P32 <- rownames(down_sigdeg_P32); length(P32) <- max_length
P44 <- rownames(down_sigdeg_P44); length(P44) <- max_length
P83 <- rownames(down_sigdeg_P83); length(P83) <- max_length
P98 <- rownames(down_sigdeg_P98); length(P98) <- max_length

# Step 3: Create a dataframe from downregulated genes and save to CSV
df_down <- data.frame(P12, P17, P26, P32, P44, P83, P98)
write.csv(df_down, "pathological_segregation/DEG/downgene.csv")

# Step 4: Draw an upset plot for downregulated genes
pdf("Figure/upset_poorly_vs_highly_downgene.pdf", width = 10, height = 6)
p <- upset(fromList(df_down), nsets = 8, nintersects = 35, order.by = "freq", keep.order = TRUE,
           sets.bar.color = cors, point.size = 3, line.size = 1,
           sets = c('P12','P17','P26','P32','P44','P83','P98'),
           main.bar.color = "#33476B", matrix.color = "#33476B", mb.ratio = c(0.6, 0.4),
           queries = list(
             list(query = intersects, params = list('P12','P17','P26','P32','P44','P83','P98'), color = "#AD1F1F", active = TRUE),
             list(query = intersects, params = list('P12','P17','P26','P32','P44','P83'), color = "#AD1F1F", active = TRUE)
           ))
print(p)
dev.off()

# Step 5: Identify commonly downregulated genes
common_down_all <- Reduce(intersect, df_down); print(common_down_all)
common_down_ATC <- Reduce(intersect, df_down[1:6]); print(common_down_ATC)

# Step 6: Combine all downregulated gene metadata and export
for (sample_id in samples) {
  obj <- get(paste0("down_sigdeg_", sample_id))
  obj$sample <- sample_id
  obj$gene <- rownames(obj)
  assign(paste0("down_sigdeg_", sample_id), obj)
}
down_sigdeg_all <- do.call(rbind, mget(paste0("down_sigdeg_", samples)))
write.csv(down_sigdeg_all, "Figure/downgene_poorly_vs_highly.csv")


#### Volcano Plot for Differentially Expressed Genes (poorly vs highly)
### Combine differential expression data for all samples
sample_list <- list(deg_P12, deg_P17, deg_P26, deg_P32, deg_P44, deg_P83, deg_P98)
sample_names <- c('P12', 'P17', 'P26', 'P32', 'P44', 'P83', 'P98')

# Annotate each data frame with sample name and gene names
for (i in seq_along(sample_list)) {
  sample_list[[i]]$sample <- sample_names[i]
  sample_list[[i]]$gene <- rownames(sample_list[[i]])
}

# Merge all samples into one data frame
diff_sample <- do.call(rbind, sample_list)
diff_sample$sample <- factor(diff_sample$sample, levels = sample_names)

# Label based on log2 fold change direction
diff_sample$label <- ifelse(diff_sample$avg_log2FC < 0, "Highly differentiated", "Poorly differentiated")

# Extract top 10 DEGs by absolute log2FC for each sample
for (i in sample_names) {
  assign(paste0("top10_", i),
         dplyr::filter(diff_sample, sample == i) %>%
           dplyr::distinct(gene, .keep_all = TRUE) %>%
           dplyr::top_n(10, abs(avg_log2FC)))
}
top10 <- do.call(rbind, mget(paste0("top10_", sample_names)))

# Mark top10 genes for label sizing
diff_sample$size <- dplyr::case_when(
  !(diff_sample$gene %in% top10$gene) ~ 1,
  diff_sample$gene %in% top10$gene ~ 2
)

# Subset for plotting
top10$sample <- factor(top10$sample, levels = sample_names)
dt <- dplyr::filter(diff_sample, size == 1)
dt$label <- factor(dt$label, levels = c("Highly differentiated", "Poorly differentiated"))

### Basic volcano-like scatter plot 
p <- ggplot() +
  geom_jitter(data = dt, aes(x = sample, y = avg_log2FC, color = label), size = 0.85, width = 0.4) +
  geom_jitter(data = top10, aes(x = sample, y = avg_log2FC, color = label), size = 1, width = 0.4) +
  scale_color_manual(name = NULL, values = c("#AD1F1F", "#33476B"))

### Add neutral background bars 
dfbar <- data.frame(x = sample_names, y = c(4.4, 4.7, 5, 3.5, 4.2, 4.5, 4.1))
dfbar1 <- data.frame(x = sample_names, y = c(-4.5, -4.5, -3.5, -4.8, -4.5, -4.6, -4.6))
dfbar$x <- factor(dfbar$x, levels = sample_names)
dfbar1$x <- factor(dfbar1$x, levels = sample_names)

# Draw background bars
p1 <- ggplot() +
  geom_col(data = dfbar, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6) +
  geom_col(data = dfbar1, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6)

# Combine with scatter
p2 <- p1 +
  geom_jitter(data = dt, aes(x = sample, y = avg_log2FC, color = label), size = 0.85, width = 0.4) +
  geom_jitter(data = top10, aes(x = sample, y = avg_log2FC, color = label), size = 1, width = 0.4) +
  scale_color_manual(name = NULL, values = c("#AD1F1F", "#33476B"))

### Add colored x-axis tiles for each cluster ###
dfcol <- data.frame(x = sample_names, y = 0, label = sample_names)
dfcol$x <- factor(dfcol$x, levels = sample_names)
mycol <- pal_npg("nrc")(7)

p3 <- p2 +
  geom_tile(data = dfcol, aes(x = x, y = y), height = 0.8, color = "black", fill = mycol, alpha = 0.6, show.legend = FALSE)

### Add gene labels to top 10 genes ###
p4 <- p3 +
  geom_text_repel(data = top10, aes(x = sample, y = avg_log2FC, label = gene),
                  size = 3, arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last"))

### Final plot with axis and color labeling 
p5 <- p4 +
  labs(x = "Sample", y = "avg_log2FC") +
  geom_text(data = dfcol, aes(x = x, y = y, label = label), size = 3, color = "white")

# Refine theme and layout
p6 <- p5 +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 13, color = "black", face = "bold"),
    axis.line.y = element_line(color = "black", size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1, 0),
    legend.text = element_text(size = 10)
  )

### Save the final volcano plot 
ggsave("Figure/volcano_poorly_vs_highly.pdf", p6, width = 10, height = 8)
saveRDS(spRNAlist,"data/spRNAlist_single_sample.rds")


#############################
###### Pathway Analysis #####
#############################

## Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(stringr)

### GO Enrichment Analysis for Upregulated Genes ###
df <- read.csv("pathological_segregation/DEG/upgene.csv", row.names = 1)
symbols_list <- as.list(as.data.frame(df))

# Convert gene symbols to Entrez IDs for each sample
gcSample <- lapply(symbols_list, function(y) {
  y <- as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                  keys = y,
                                                  columns = 'ENTREZID',
                                                  keytype = 'SYMBOL')[, 2]))
  y
})

# Perform GO enrichment analysis (Biological Process category)
formula_res <- compareCluster(
  gcSample,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Dot plot for visualization
p <- dotplot(formula_res, showCategory = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
  scale_fill_gradientn(values = seq(0, 1, 0.2),
                       colours = c('#CC3333', '#FFFF99', '#6699CC'))

ggsave("Figure/upgene_GOenrich_poorly_vs_highly.pdf", p, width = 8, height = 8)

### GO Enrichment Analysis for Downregulated Genes ###
df <- read.csv("pathological_segregation/DEG/downgene.csv", row.names = 1)
symbols_list <- as.list(as.data.frame(df))

gcSample <- lapply(symbols_list, function(y) {
  y <- as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                  keys = y,
                                                  columns = 'ENTREZID',
                                                  keytype = 'SYMBOL')[, 2]))
  y
})

formula_res <- compareCluster(
  gcSample,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

p <- dotplot(formula_res, showCategory = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
  scale_fill_gradientn(values = seq(0, 1, 0.2),
                       colours = c('#CC3333', '#FFFF99', '#6699CC'))

ggsave("Figure/downgene_GOenrich_poorly_vs_highly.pdf", p, width = 8, height = 8)

# Save single-sample object for reuse
saveRDS(spRNAlist, "data/spRNAlist_single_sample.rds")


### GSVA Pathway Analysis Using HALLMARK Gene Sets
library(GSVA)
library(msigdbr)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)
library(ggpubr)

# Clean environment
rm(list = ls())
gc()

# Load spatial object list
spRNAlist <- readRDS("data/spRNAlist_single_sample.rds")

sample <- c("P12", "P17", "P26", "P32", "P44", "P57", "P83", "P98")
for (i in 1:8) {
  spRNAlist[[i]]@images[-which(names(spRNAlist[[i]]@images) == sample[i])] <- NULL
}

# Merge into a single Seurat object
spRNA <- merge(spRNAlist[[1]], spRNAlist[2:length(spRNAlist)])
saveRDS(spRNA, "data/spRNA_single_sample.rds")

# Construct combined region label
spRNA$region2 <- paste0(spRNA$sample, "_", spRNA$region1)

# Subset highly vs poorly differentiated tumor cells
# FIXED below with corrected syntax
spRNAsub <- subset(spRNA, region1 %in% c("highly_differentiated", "poorly_differentiated"))

# Load HALLMARK gene sets from .gmt file
genesets <- read.gmt("/datadisk1/person/NICK/spRNA/pathological_segregation/h.all.v7.5.1.symbols.gmt")
genesets <- split(genesets$gene, genesets$term)

# Set identity class by combined region
Idents(spRNAsub) <- "region2"

# Compute average expression matrix
expr <- AverageExpression(spRNAsub, assays = "Spatial", slot = "data")[[1]]
expr <- expr[rowSums(expr) > 0, ] %>% as.matrix()

# Perform GSVA
gsva.res <- gsva(expr, genesets, method = "gsva")

# Heatmap of GSVA results
pheatmap::pheatmap(gsva.res, show_colnames = T,
                   scale = "row", angle_col = "45", cluster_cols = F,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# Define sample grouping for comparison (P = poorly, A = highly differentiated)
group_list <- data.frame(
  sample = colnames(gsva.res),
  group = c('P', 'A', 'A', 'P', 'A', 'P', 'A', 'P', 'P', 'A', 'A', 'P', 'P', 'A', 'P', 'A')
)

# Differential pathway analysis using limma
design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva.res)

contrast.matrix <- makeContrasts(A - P, levels = design)
fit <- lmFit(gsva.res, design)
fit2 <- contrasts.fit(fit, contrast.matrix) %>% eBayes()
allDiff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

# Prepare data for visualization
pathway <- str_replace(rownames(allDiff), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = allDiff$t)
cutoff <- 0
df$group <- cut(df$score, breaks = c(-Inf, cutoff, Inf), labels = c(1, 2))

# Reorder for plotting
sortdf <- df[order(df$score), ]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)

# Bar plot with annotated text
pdf("Figure/GSVA_pathway_poorly_vs_highly.pdf")
ggplot(sortdf, aes(ID, score, fill = group)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_manual(values = c('#35496E', '#AD1F1F'), guide = "none") +
  geom_hline(yintercept = c(-1, 1), color = "white", linetype = 2, size = 0.3) +
  geom_text(data = subset(df, score < 0),
            aes(x = ID, y = 0.1, label = ID, color = group),
            size = 4, hjust = "inward") +
  geom_text(data = subset(df, score > 0),
            aes(x = ID, y = -0.1, label = paste0(" ", ID), color = group),
            size = 4, hjust = "outward") +
  scale_colour_manual(values = c("black", "black"), guide = "none") +
  labs(x = "", y = "t value of GSVA Score", title = "Malignant cell: ATC vs. DTC") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
dev.off()




#############################
##### Gene Visualization ####
#############################

# Clear environment and free memory
rm(list = ls())
gc()

# Load single-sample spatial transcriptomics data
spRNAlist <- readRDS("data/spRNAlist_single_sample.rds")

# Sample names
sample <- c("P12", "P17", "P26", "P32", "P44", "P83", "P98")

# Keep only the corresponding image for each sample
for (i in 1:7) {
  spRNAlist[[i]]@images[-which(names(spRNAlist[[i]]@images) == sample[i])] <- NULL
}

# Merge individual Seurat objects
spRNA <- merge(spRNAlist[[1]], spRNAlist[2:length(spRNAlist)])

# Define a function to generate spatial plots for a given gene
plot_spatial_feature <- function(spRNA, features) {
  p1 <- SpatialFeaturePlot(spRNA, features = features, images = "P12", pt.size.factor = 3.5) + 
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) + ggtitle("P12")
  p2 <- SpatialFeaturePlot(spRNA, features = features, images = "P17", pt.size.factor = 4) + 
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) + ggtitle("P17")
  p3 <- SpatialFeaturePlot(spRNA, features = features, images = "P26", pt.size.factor = 2) + 
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) + ggtitle("P26")
  p4 <- SpatialFeaturePlot(spRNA, features = features, images = "P32", pt.size.factor = 4) + 
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) + ggtitle("P32")
  p5 <- SpatialFeaturePlot(spRNA, features = features, images = "P44", pt.size.factor = 4) + 
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) + ggtitle("P44")
  p6 <- SpatialFeaturePlot(spRNA, features = features, images = "P83", pt.size.factor = 2) + 
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) + ggtitle("P83")
  p7 <- SpatialFeaturePlot(spRNA, features = features, images = "P98", pt.size.factor = 3.5) + 
    theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) + ggtitle("P98")
  
  combined_plot <- (p1 | p2 | p3 | p4) / (p5 | p6 | p7)
  return(combined_plot)
}

# Gene list to be visualized
gene <- c("ID4", "TG", "PDCD4", "NTRK2", "DAPK2", "PAX8", "S100A1", "ATP1B1", 
          "EPCAM", "SLC34A2", "CD24", "IYD", "TG", "PDCD4", "LMO3", "PEBP1", "DAPK2")

# Output spatial plots to a PDF
pdf('/datadisk1/person/NICK/spRNA/pathological_segregation/DEG/gene.pdf')
for (i in seq_along(gene)) {
  p <- plot_spatial_feature(spRNA, features = gene[i])
  print(p)
}
dev.off()


###### Violin Plots 

# Prepare expression matrix for selected genes
data <- cbind(spRNA@meta.data, t(spRNA@assays$SCT@data[c("TG", "PAX8", "TSHR", "NKX2-1"), ]))

# Filter out non-epithelial regions
data <- subset(data, data$region != "NSE" & data$region != "Tertiary")

# Set region factor levels
spRNA$region <- factor(spRNA$region, levels = c("ATC", "HGFCTC", "DTC", "Stroma"))

# Define color palette
celltype_colors <- c("#E64B35FF", "#F39B7FFF", "#00A087FF", "#3C5488FF")

# Gene list for plotting
gene <- c("TG", "PAX8", "TSHR", "NKX2-1")
p_list <- list()

# Generate violin plots for each gene
for (i in seq_along(gene)) {
  data$gene <- data[, gene[i]]
  p_list[[i]] <- ggplot(data, aes(sample, gene, fill = region)) +
    geom_violin(width = 2, position = position_dodge(width = 0.6)) +
    scale_fill_manual(values = celltype_colors) +
    labs(x = NULL, y = NULL) + theme_test() +
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          legend.position = "right",
          strip.background = element_blank(),
          strip.text.x = element_text(color = "black", face = "bold", size = 11)) +
    ggtitle(gene[i])
}

# Save violin plots to a PDF
ggsave("Figure/Differentiation_genes.pdf", p_list[[1]] / p_list[[2]] / p_list[[3]] / p_list[[4]], 
       width = 10, height = 8)







