i# ============================================
# Cell2location pipeline for spatial mapping
# ============================================

import sys
import os
import scanpy as sc
import cell2location
import numpy as np
import matplotlib as mpl
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

# ---------------------------------------------------
# Function 1: Train regression model from scRNA-seq
# ---------------------------------------------------
def ModTrain(scfile, workdir):
    # Load single-cell data
    adata_sc = sc.read_h5ad(scfile)

    # Filter genes based on expression criteria
    selected = filter_genes(
        adata_sc,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12
    )
    adata_sc = adata_sc[:, selected].copy()

    # Set up AnnData for cell2location training
    RegressionModel.setup_anndata(
        adata=adata_sc,
        batch_key='Samples',
        labels_key='cell_type'
    )

    # Train the regression model to derive reference signatures
    mod = RegressionModel(adata_sc)
    mod.train(max_epochs=250, train_size=1)

    # Export posterior estimates
    adata_sc = mod.export_posterior(
        adata_sc,
        sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
    )
    adata_sc = mod.export_posterior(
        adata_sc,
        use_quantiles=True,
        add_to_varm=["q05", "q50", "q95", "q0001"],
        sample_kwargs={'batch_size': 2500}
    )

    # Save trained model and results
    moddir = f"{workdir}/reference_signatures/"
    h5file = f"{workdir}/sc_cell2location.h5ad"
    mod.save(moddir, overwrite=True)
    adata_sc.write(h5file)

# ---------------------------------------------------
# Function 2: Perform deconvolution on spatial data
# ---------------------------------------------------
def Convolution(stfile, sch5file, scmodsig, workdir):
    # Load spatial transcriptomics data
    adata_vis = sc.read_h5ad(stfile)

    # Remove mitochondrial genes
    adata_vis.var['MT_gene'] = [gene.startswith('mt-') for gene in adata_vis.var.index]
    adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

    # Load trained single-cell reference
    adata_sc = sc.read_h5ad(sch5file)

    # Extract inferred gene expression per cell type
    if 'means_per_cluster_mu_fg' in adata_sc.varm:
        inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[
            f'means_per_cluster_mu_fg_{i}' for i in adata_sc.uns['mod']['factor_names']
        ]].copy()
    else:
        inf_aver = adata_sc.var[[
            f'means_per_cluster_mu_fg_{i}' for i in adata_sc.uns['mod']['factor_names']
        ]].copy()

    inf_aver.columns = adata_sc.uns['mod']['factor_names']

    # Align genes between spatial and reference datasets
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    adata_sc = adata_sc[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # Set up and train spatial deconvolution model
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, layer="raw_count")
    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=30,
        detection_alpha=20
    )
    mod.view_anndata_setup()
    mod.train(max_epochs=20000, batch_size=None, train_size=1)
    mod.plot_history(1000)

    # Export posterior results
    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
    )

    # Save spatial results
    mod.plot_QC()
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    stmoddir = f"{workdir}/cell2location_map/"
    stfile_out = f"{workdir}/st_cell2location_res.h5ad"
    mod.save(stmoddir, overwrite=True)
    adata_vis.write(stfile_out)

    # Generate output figures
    picdir = f"{workdir}/picture"
    if not os.path.isdir(picdir):
        os.mkdir(picdir)

    for i in set(adata_sc.obs.cell_type):
        picfile = f"{picdir}/{i}.png"
        with mpl.rc_context({'axes.facecolor': 'black', 'figure.figsize': [4.5, 5]}):
            sc.pl.spatial(
                adata_vis,
                cmap='magma',
                color=[i],
                ncols=4,
                size=1.3,
                img_key='lowres',
                vmin=0,
                vmax='p99.2',
                save=picfile
            )

# ---------------------------------------------------
# Run the full pipeline
# ---------------------------------------------------

# Step 1 – Train model from single-cell data
scfile = "/data/total.singlecell.h5ad"
workdir = "/data/scmod"
ModTrain(scfile, workdir)

# Step 2 – Apply model to multiple spatial transcriptomics samples
st_files = [
    "/data/st/P1.st.h5ad", "/data/st/P2.st.h5ad", "/data/st/P3.st.h5ad", "/data/st/P4.st.h5ad",
    "/data/st/P5.st.h5ad", "/data/st/P6.st.h5ad", "/data/st/P7.st.h5ad", "/data/st/P8.st.h5ad"
]
workdirs = [
    "/data/convolution1", "/data/convolution2", "/data/convolution3", "/data/convolution4",
    "/data/convolution5", "/data/convolution6", "/data/convolution7", "/data/convolution8"
]
sch5file = "/data/scmod/sc_cell2location.h5ad"
scmodsig = "/data/scmod/reference_signatures/"

# Loop over spatial samples
for stfile, wd in zip(st_files, workdirs):
    Convolution(stfile, sch5file, scmodsig, wd)
