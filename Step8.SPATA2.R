# =============================================
# Spatial gene expression and trajectory analysis using SPATA2
# =============================================

library(SPATA2)

# ---------------------------------------------------
# Function: Convert Seurat object to SPATA2 format
# ---------------------------------------------------
SPSeuratToSPATAData <- function(STSeurat, SampleName, aname) {
  STSeurat <- SPATA2::asSPATA2(
    STSeurat,
    sample_name = SampleName,
    platform = "VisiumSmall",
    assay_name = aname,
    assay_modality = "gene",
    img_name = SampleName,
    img_scale_fct = "lowres"
  )
  return(STSeurat)
}

# Convert Visium Seurat objects to SPATA2 format for each sample
P12StSPATA <- SPSeuratToSPATAData(P12STSeurat, "P12", "Spatial")
P17StSPATA <- SPSeuratToSPATAData(P17STSeurat, "P17", "Spatial")
P26StSPATA <- SPSeuratToSPATAData(P26STSeurat, "P26", "Spatial")
P32StSPATA <- SPSeuratToSPATAData(P32STSeurat, "P32", "Spatial")
P44StSPATA <- SPSeuratToSPATAData(P44STSeurat, "P44", "Spatial")
P83StSPATA <- SPSeuratToSPATAData(P83STSeurat, "P83", "Spatial")
P98StSPATA <- SPSeuratToSPATAData(P83STSeurat, "P98", "Spatial") # Note: input is P83STSeurat again – is this intentional?

# ---------------------------------------------------
# Create spatial trajectories for each SPATA2 object
# ---------------------------------------------------
P12StSPATA <- createSpatialTrajectories(P12StSPATA)
P17StSPATA <- createSpatialTrajectories(P17StSPATA)
P26StSPATA <- createSpatialTrajectories(P26StSPATA)
P32StSPATA <- createSpatialTrajectories(P32StSPATA)
P44StSPATA <- createSpatialTrajectories(P44StSPATA)
P83StSPATA <- createSpatialTrajectories(P83StSPATA)
P98StSPATA <- createSpatialTrajectories(P98StSPATA)

# ---------------------------------------------------
# Plot PDCD4 spatial expression along trajectories
# ---------------------------------------------------
# For each sample, overlay the defined trajectory and visualize PDCD4 expression
P12TraAddOn <- ggpLayerSpatialTrajectories(P12StSPATA, ids = "tra1")
plotSurfaceComparison(
  object = P12StSPATA,
  color_by = "PDCD4",
  pt_clrsp = "Reds 3",
  outline = TRUE,
  nrow = 1
) + P12TraAddOn

P17TraAddOn <- ggpLayerSpatialTrajectories(P17StSPATA, ids = "tra1")
plotSurfaceComparison(P17StSPATA, "PDCD4", pt_clrsp = "Reds 3", outline = TRUE, nrow = 1) + P17TraAddOn

P26TraAddOn <- ggpLayerSpatialTrajectories(P26StSPATA, ids = "tra1")
plotSurfaceComparison(P26StSPATA, "PDCD4", pt_clrsp = "Reds 3", outline = TRUE, nrow = 1) + P26TraAddOn

P32TraAddOn <- ggpLayerSpatialTrajectories(P32StSPATA, ids = "tra1")
plotSurfaceComparison(P32StSPATA, "PDCD4", pt_clrsp = "Reds 3", outline = TRUE, nrow = 1) + P32TraAddOn

P44TraAddOn <- ggpLayerSpatialTrajectories(P44StSPATA, ids = "tra1")
plotSurfaceComparison(P44StSPATA, "PDCD4", pt_clrsp = "Reds 3", outline = TRUE, nrow = 1) + P44TraAddOn

P83TraAddOn <- ggpLayerSpatialTrajectories(P83StSPATA, ids = "tra1")
plotSurfaceComparison(P83StSPATA, "PDCD4", pt_clrsp = "Reds 3", outline = TRUE, nrow = 1) + P83TraAddOn

P98TraAddOn <- ggpLayerSpatialTrajectories(P98StSPATA, ids = "tra1")
plotSurfaceComparison(P98StSPATA, "PDCD4", pt_clrsp = "Reds 3", outline = TRUE, nrow = 1) + P98TraAddOn

# ---------------------------------------------------
# Line plots: Expression trends of PDCD4 and TYMPTAM along trajectories
# ---------------------------------------------------
# For each sample, plot expression of PDCD4 and TYMPTAM along spatial trajectory ("tra1")
# Color adjusted manually for visual clarity; no smoothing errors shown

plotStsLineplot(P12StSPATA, variables = c("PDCD4", "TYMPTAM"), id = "tra1", clrp_adjust = c("red", "blue"), nrow = 1, display_facets = FALSE, smooth_se = FALSE, smooth_span = 0.5) + theme_cowplot()
plotStsLineplot(P17StSPATA, variables = c("PDCD4", "TYMPTAM"), id = "tra1", clrp_adjust = c("red", "blue"), nrow = 1, display_facets = FALSE, smooth_se = FALSE, smooth_span = 0.5) + theme_cowplot()
plotStsLineplot(P26StSPATA, variables = c("PDCD4", "TYMPTAM"), id = "tra1", clrp_adjust = c("red", "blue"), nrow = 1, display_facets = FALSE, smooth_se = FALSE, smooth_span = 0.5) + theme_cowplot()
plotStsLineplot(P32StSPATA, variables = c("PDCD4", "TYMPTAM"), id = "tra1", clrp_adjust = c("red", "blue"), nrow = 1, display_facets = FALSE, smooth_se = FALSE, smooth_span = 0.5) + theme_cowplot()
plotStsLineplot(P44StSPATA, variables = c("PDCD4", "TYMPTAM"), id = "tra1", clrp_adjust = c("red", "blue"), nrow = 1, display_facets = FALSE, smooth_se = FALSE, smooth_span = 0.6) + theme_cowplot()
plotStsLineplot(P83StSPATA, variables = c("PDCD4", "TYMPTAM"), id = "tra1", clrp_adjust = c("red", "blue"), nrow = 1, display_facets = FALSE, smooth_se = FALSE, smooth_span = 0.7) + theme_cowplot()
plotStsLineplot(P98StSPATA, variables = c("PDCD4", "TYMPTAM"), id = "tra1", clrp_adjust = c("red", "blue"), nrow = 1, display_facets = FALSE, smooth_se = FALSE, smooth_span = 0.8) + theme_cowplot()
