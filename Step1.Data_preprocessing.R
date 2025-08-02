######################################################################
############Data preprocessing and quality control####################
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
rm(list=ls())
gc()

#### Load Spatial Transcriptomics Data ####
# Set working directory
setwd("/datadisk1/person/NICK/spRNA")

# List all sample folder names under 'rawdata' directory
samples_name <- list.files("./rawdata")

# Initialize an empty list to store individual Seurat objects
spRNAlist <- list()

# Loop through each sample and load the spatial data
for (i in 1:length(samples_name)) {
  sample_path <- paste0("./rawdata/", samples_name[i])
  sample_id <- paste0("P", str_sub(samples_name[i], -2, -1))  # Use last two characters as sample ID
  
  spRNAlist[[i]] <- Load10X_Spatial(
    data.dir = sample_path,
    filename = paste0(samples_name[i], ".h5"),
    slice = sample_id
  )
  
  spRNAlist[[i]]$orig.ident <- sample_id
}

# Merge all individual Seurat objects into one
spRNA <- merge(spRNAlist[[1]], spRNAlist[2:length(spRNAlist)])

# View metadata
View(spRNA@meta.data)


#### Clinical Metadata Annotation ####

# Copy sample ID
spRNA$sample <- spRNA$orig.ident

# Annotate gender
spRNA$gender <- "male"
spRNA$gender[spRNA$sample %in% c("P98", "P26", "P12", "P44", "P32")] <- "female"

# Annotate pathology types:
# PDTC: poorly differentiated thyroid carcinoma
# SC: squamous-type anaplastic thyroid carcinoma
# ATC: sarcomatoid-type anaplastic thyroid carcinoma (default)
spRNA$pathology <- "ATC"
spRNA$pathology[spRNA$sample %in% c("P17", "P83", "P26")] <- "SC"
spRNA$pathology[spRNA$sample == "P98"] <- "PDTC"

# Annotate spatial distribution between undifferentiated and papillary regions
# Types: "mixed" (intermingled), "separation" (clearly divided)
spRNA$distribution <- "separation"
spRNA$distribution[spRNA$sample %in% c("P83", "P44", "P32")] <- "mixed"

# Annotate survival time (in days): death date - surgery date
spRNA$survivaltime <- NA  # initialize with NA

spRNA$survivaltime[spRNA$sample == "P12"] <- 110
spRNA$survivaltime[spRNA$sample == "P17"] <- 216
spRNA$survivaltime[spRNA$sample == "P26"] <- 59
spRNA$survivaltime[spRNA$sample == "P32"] <- 58
spRNA$survivaltime[spRNA$sample == "P44"] <- 37
spRNA$survivaltime[spRNA$sample == "P83"] <- 470
spRNA$survivaltime[spRNA$sample == "P98"] <- 829


# Save raw object before QC
dir.create("data", showWarnings = FALSE)
saveRDS(spRNA, "data/before_QC_spRNA.rds")







#### Quality Control for Spatial Transcriptomics Data ####

# Clear workspace
rm(list = ls())
gc()

# Load the raw Seurat object (before QC)
spRNA <- readRDS("/datadisk1/person/NICK/spRNA/data/before_QC_spRNA.rds")

# Calculate the percentage of mitochondrial genes per spot
spRNA[["percent.mt"]] <- PercentageFeatureSet(spRNA, pattern = "^MT-")

# ---- QC Metrics Before Filtering ----
# Define custom colors
celltype_colors <- pal_npg("nrc")(8)

# Violin plots: before filtering
p1 <- VlnPlot(spRNA, features = "nCount_Spatial", group.by = "orig.ident", pt.size = 0) + 
  NoLegend() + scale_fill_manual(values = celltype_colors)
p2 <- VlnPlot(spRNA, features = "nFeature_Spatial", group.by = "orig.ident", pt.size = 0) + 
  NoLegend() + scale_fill_manual(values = celltype_colors)
p3 <- VlnPlot(spRNA, features = "percent.mt", group.by = "orig.ident", pt.size = 0) + 
  NoLegend() + scale_fill_manual(values = celltype_colors)

# Combine and save violin plots
dir.create("QC", showWarnings = FALSE)
(p1 | p2 | p3)
ggsave("QC/beforeQC.pdf", p1 | p2 | p3, width = 10, height = 4)

# ---- SpatialFeaturePlots: before QC ----

samples_name <- names(spRNA@images)

# Plot: nCount_Spatial
p1 <- SpatialFeaturePlot(spRNA, features = "nCount_Spatial", images = samples_name[1:4], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) +
  ggtitle("P12")
p2 <- SpatialFeaturePlot(spRNA, features = "nCount_Spatial", images = samples_name[1:4], pt.size.factor = 2) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))
p3 <- SpatialFeaturePlot(spRNA, features = "nCount_Spatial", images = samples_name[5:8], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))

p1 / p2 / p3
ggsave("QC/nCount_beforeQC.pdf", p1 / p2 / p3, width = 10, height = 10)

# Plot: nFeature_Spatial
p1 <- SpatialFeaturePlot(spRNA, features = "nFeature_Spatial", images = samples_name[1:4], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) +
  ggtitle("P12")
p2 <- SpatialFeaturePlot(spRNA, features = "nFeature_Spatial", images = samples_name[1:4], pt.size.factor = 2) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))
p3 <- SpatialFeaturePlot(spRNA, features = "nFeature_Spatial", images = samples_name[5:8], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))

p1 / p2 / p3
ggsave("QC/nFeature_beforeQC.pdf", p1 / p2 / p3, width = 10, height = 10)

# Plot: percent.mt
p1 <- SpatialFeaturePlot(spRNA, features = "percent.mt", images = samples_name[1:4], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8)) +
  ggtitle("P12")
p2 <- SpatialFeaturePlot(spRNA, features = "percent.mt", images = samples_name[1:4], pt.size.factor = 2) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))
p3 <- SpatialFeaturePlot(spRNA, features = "percent.mt", images = samples_name[5:8], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))

p1 / p2 / p3
ggsave("QC/percent.mt_beforeQC.pdf", p1 / p2 / p3, width = 10, height = 10)

# ---- Filtering: Apply QC thresholds ----

minCount <- 500
minFeature <- 300

spRNA2 <- subset(spRNA, nCount_Spatial > minCount & nFeature_Spatial > minFeature)

# ---- Post-QC Visualizations ----

samples_name <- names(spRNA2@images)

# nCount_Spatial
p1 <- SpatialFeaturePlot(spRNA2, features = "nCount_Spatial", images = samples_name[1:4], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))
p2 <- SpatialFeaturePlot(spRNA2, features = "nCount_Spatial", images = samples_name[1:4], pt.size.factor = 2) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))
p3 <- SpatialFeaturePlot(spRNA2, features = "nCount_Spatial", images = samples_name[5:8], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))

p1 / p2 / p3
ggsave("QC/nCount_afterQC.pdf", p1 / p2 / p3, width = 10, height = 10)

# nFeature_Spatial
p1 <- SpatialFeaturePlot(spRNA2, features = "nFeature_Spatial", images = samples_name[1:4], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))
p2 <- SpatialFeaturePlot(spRNA2, features = "nFeature_Spatial", images = samples_name[1:4], pt.size.factor = 2) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))
p3 <- SpatialFeaturePlot(spRNA2, features = "nFeature_Spatial", images = samples_name[5:8], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))

p1 / p2 / p3
ggsave("QC/nFeature_afterQC.pdf", p1 / p2 / p3, width = 10, height = 10)

# percent.mt
p1 <- SpatialFeaturePlot(spRNA2, features = "percent.mt", images = samples_name[1:4], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))
p2 <- SpatialFeaturePlot(spRNA2, features = "percent.mt", images = samples_name[1:4], pt.size.factor = 2) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))
p3 <- SpatialFeaturePlot(spRNA2, features = "percent.mt", images = samples_name[5:8], pt.size.factor = 3.5) +
  theme(legend.position = "top", legend.text = element_text(size = 5), legend.title = element_text(size = 8))

p1 / p2 / p3
ggsave("QC/percent.mt_afterQC.pdf", p1 / p2 / p3, width = 10, height = 10)

# Save filtered Seurat object
saveRDS(spRNA2, "data/after_QC_spRNA.rds")


















