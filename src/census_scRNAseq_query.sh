#!/usr/bin/env python

import cellxgene_census
import scanpy as sc
import numpy as np
import scvi
from scipy.sparse import csr_matrix
#import tiledbsoma
import pickle
import pandas as pd
import subprocess
import scipy
import random
import os
import pickle
import skmisc
import sklearn
import matplotlib
import matplotlib.pyplot as plt
from xxhash import xxh3_64_intdigest
from base64 import b85encode

census = cellxgene_census.open_soma(census_version="latest")

file_path = "/staging/leuven/stg_00019/research/Kate/emseq/data/Ensemble_ENSG.txt"

with open(file_path, 'r') as file:
    lines = file.readlines()

lines = [line.strip() for line in lines]


#tabula sapiens
dataset_ids = ["53d208b0-2cfd-4366-9866-c3c6114081bc"]

#tabula sapiens + vento-tormo (2023) trophoblasts + sle
#dataset_ids = ["53d208b0-2cfd-4366-9866-c3c6114081bc", "ecf2e08e-2032-4a9e-b466-b65b395f4a02", "218acb0f-9f2f-4f76-b90b-15a4b7c7f629"]

#tabula sapiens + vento-tormo (2023) trophoblasts
#dataset_ids = ["53d208b0-2cfd-4366-9866-c3c6114081bc", "ecf2e08e-2032-4a9e-b466-b65b395f4a02"]

#query cesus for ensembl genes and datasets and no duplicates and not Smart-seq data
adata = cellxgene_census.get_anndata(
    census, organism="Homo sapiens", var_value_filter = f"feature_id in {lines}", obs_value_filter=f"dataset_id in {dataset_ids} and assay != 'Smart-seq2'"
)

#qc plots
#os.chdir('/staging/leuven/stg_00019/research/Kate/emseq/data/census')
#pd.set_option('display.max_columns', None)
#pd.set_option('display.max_rows', None)

#print the number of cells per assay
#print(adata.obs.assay.value_counts())
#print(adata.obs.cell_type.value_counts())
#print(adata.obs.dataset_id.value_counts())
#print(adata.obs.tissue_type.value_counts())
#print(adata.obs.disease.value_counts())
#print(adata.obs.donor_id.value_counts())


census.close()
del census

#ensemble gene filter removes MT genes
#adata.var["mt"] = adata.var["feature_name"].str.startswith("MT-")
#sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)

#plot qc metrics
#valid_dataset_ids = [dataset_id for dataset_id in dataset_ids if (adata.obs['dataset_id'] == dataset_id).sum() > 0]

#sc.pl.violin(
#    adata,
#    ["n_genes_by_counts", "total_counts"],
#    jitter=0.4,
#    multi_panel=True,
#    groupby = 'dataset_id',
#    rotation=90,
#    save='-qc.datasetid.png'
#)

#scatter
#for dataset_id in valid_dataset_ids:
#    subset_adata = adata[adata.obs['dataset_id'] == dataset_id]
#    sc.pl.scatter(
#        subset_adata,
#        x='total_counts',
#        y='n_genes_by_counts',
#        save=f'{dataset_id}-qc.datasetid.png'
#    )

#filter out normal control cells from sle dataset 
#adata = adata[~((adata.obs['disease'] == 'normal') & (adata.obs['dataset_id'] == '218acb0f-9f2f-4f76-b90b-15a4b7c7f629'))]
#subset to disease flare cases!
#mask_normal = adata.obs['disease'] == "normal"
#mask_sle_flare = (adata.obs['disease'] == "systemic lupus erythematosus") & (adata.obs['donor_id'].str.contains('FLARE'))

# Combine the masks using logical OR
#combined_mask = mask_normal | mask_sle_flare

# Subset the AnnData object based on the combined mask
#adata = adata[combined_mask, :]

#ANNOTATE WITH CELL TYPES AND BATCHES FROM VTTROPH ORIGINAL DATA OBJECT
#adatavento = sc.read_h5ad("/staging/leuven/stg_00019/research/Kate/emseq/data/vento-tormo/adata_all_donors_trophoblast_raw_counts_in_raw_normlog_counts_in_X_for_download_UPD_20230307.h5ad")
#print(adatavento)

#convert cell barcode ids to hashed ids in census objects
#got code to transform cell barcodes to hashed ids here: https://github.com/chanzuckerberg/single-cell-curation/blob/main/cellxgene_schema_cli/cellxgene_schema/utils.py
#def hash_identifier(identifier):
    # Convert identifier to hash digest
#    hash_digest = xxh3_64_intdigest(identifier.encode())
    # Encode hash digest into Base85
#    encoded_hash = b85encode(hash_digest.to_bytes(8, "big")).decode("ascii")
#    return encoded_hash

# Sample list of cell barcodes (replace with your actual list)
#cell_barcodes = adatavento.obs_names.tolist()

# Hash each cell barcode and create DataFrame
#hashed_ids = [hash_identifier(barcode) for barcode in cell_barcodes]
#df = pd.DataFrame({'observation_joinid': hashed_ids, 'Original_ID': cell_barcodes})

# Merge the hashed IDs DataFrame with adata.obs on 'observation_joinid'
#adata.obs = adata.obs.merge(df, on='observation_joinid', how='left')
#adatavento.obs['Original_ID'] = adatavento.obs.index
#adata.obs = adata.obs.merge(adatavento.obs, on='Original_ID', how='left')

#merge vento cell type column with census cell type column
#adata.obs['final_annot_troph'] = adata.obs['final_annot_troph'].astype(str)
#adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)

#adata.obs['final_annot_troph'] = adata.obs['final_annot_troph'].replace("nan", float("NaN"))

#adata.obs['final_annot_troph'] = adata.obs['final_annot_troph'].fillna(adata.obs['cell_type'])

#create batches for correction
#adata.obs['donor_id_assay'] = adata.obs['donor_id'].astype(str) + '_' + adata.obs['assay'].astype(str)
#adata.obs['sample'] = adata.obs['sample'].astype(str)
#adata.obs['donor_id_assay'] = adata.obs['donor_id_assay'].astype(str)

#adata.obs['sample'] = adata.obs['sample'].replace("nan", float("NaN"))
#adata.obs['sample'] = adata.obs['sample'].fillna(adata.obs['donor_id_assay'])

#donor batches with less than 100 cells
#exclude_donor_ids = [
#    "H2_fetus",
#    "P13_mother",
#    "P14_mother",
#    "7_mother",
#    "6_fetus",
#    "7_fetus",
#    "Hrv43_mother",
#    "Hrv46_mother",
#    "Hrv99_mother",
#    "P34_mother",
#    "12_mother",
#    "R0_mother",
#    "10_mother",
#    "8_mother"
#]

exclude_batch_ids = [
    "FCA7196225",
    "Pla_HDBR8624431",
    "Pla_HDBR10142864",
    "Pla_HDBR10142863",
    "Pla_HDBR10142865",
    "Pla_HDBR10142768",
    "FCA7167226",
    "FCA7167224",
    "FCA7167222",
    "Pla_HDBR10142767",
    "WSSS_PLA8811070",
    "FCA7474062",
    "Pla_HDBR10701668",
    "Pla_HDBR10701667",
    "WSSS_PLA8811069",
    "FCA7511882",
    "WSSS_PLA8811068",
    "FCA7196219",
    "FCA7167223",
    "FCA7167219",
    "Pla_HDBR8624430",
    "FCA7511881",
    "FCA7167221",
    "FCA7196224",
    "FCA7196218",
    "Pla_HDBR8715514",
    "Pla_HDBR8768477",
    "Pla_HDBR8715512"
]

# Subset the adata object to exclude donor_ids
#adata = adata[~adata.obs['sample'].isin(exclude_batch_ids)].copy()
#print(adata.obs['sample'].value_counts())

#basic normalization
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

#select variable genes
#sc.pp.highly_variable_genes(
#    adata,
#    n_top_genes=1000,
#    flavor="seurat_v3",
#    layer="counts",
#    batch_key="dataset_id",
#    subset=True,
#)

#visualize pre integration
#sc.tl.pca(adata)
#sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
#sc.tl.umap(adata)
#sc.pl.umap(adata, 
#    color="dataset_id",
#    legend_loc="on data",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.preintegration.tabulasapiens.only.dataset_id.png'
#)

#sc.pl.umap(
#    adata, 
#    color="tissue_general",
#    legend_loc="on data",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1, 
#    save='umap.preintegration.tbs.tissue_general.png'
#)

#perform integration
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="donor_id")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb", n_hidden=50)

vae.train(max_epochs=100)
adata.obsm["X_scVI"] = vae.get_latent_representation()
adata.layers["scVI_normalized"] = vae.get_normalized_expression()
adata.write('/staging/leuven/stg_00019/research/Kate/emseq/data/census/tbs_censusversion2024-04-29_norm_log_scale_integrated_on_donorid_allgenes_normalized_batchcorrected_geneexpression_in_scVI_normalized.h5ad')
#placed lower when outputting average expression#
adata.X = adata.layers["scVI_normalized"]
#adata=sc.read_h5ad('/staging/leuven/stg_00019/research/Kate/emseq/data/census/tbs_vt2023trophoblasts_censusversion2024-04-29_norm_log_scale_integrated_on_donorid_morethan100cells_nohighlyvariablegenes_allgenes_normalized_batchcorrected_geneexpression_in_scVI_normalized.h5ad')
#visualize datasets post integration
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color="donor_id",
    legend_loc="right margin",
    frameon=False,
    legend_fontsize=5,
    legend_fontoutline=1,
    save='umap.postintegration_on_donorid.tbs.donor_idpng'
)

sc.pl.umap(
    adata,
    color="tissue_general",
    legend_loc="right margin",
    frameon=False,
    legend_fontsize=5,
    legend_fontoutline=1,
    save='umap.postintegration_on_donorid.tbs.tissue_general.png'
)

#read in anndata object generated from above workflow
#adata = sc.read_h5ad("/staging/leuven/stg_00019/research/Kate/emseq/data/census/tbs_vt2023trophoblasts_sle_censusversion2024-04-29_norm_log_scale_integrated_nohighlyvariablegenes_allgenes.h5ad")

#adata = sc.read_h5ad("/staging/leuven/stg_00019/research/Kate/emseq/data/census/tbs_vt2023trophoblasts_sle_censusversion2024-04-29_norm_log_scale_integrated_nohighlyvariablegenes_allgenes.h5ad")
#print(adata)
#plt.figure(figsize=(50, 50))
#sc.pl.umap(
#    adata,
#    color="cell_type",
#    legend_loc="right margin",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.postintegration.tbs.vento.sle.cell_type.allgenes.png'
#)

#adata_placenta = adata[adata.obs['tissue_general'] == 'placenta'].copy()
#print(adata_placenta)
#print(adata_placenta.obs['cell_type'].dtype)
#print(adata_placenta.obs['tissue_general'].dtype)
#adata_placenta.obs['cell_type'] = pd.Categorical(adata_placenta.obs['cell_type'])
#adata_placenta.obs['donor_id'] = pd.Categorical(adata_placenta.obs['donor_id'])

# Regenerate neighbors and UMAP for the placenta subset

#sc.pp.neighbors(adata_placenta, use_rep="X_scVI")
#sc.tl.umap(adata_placenta)

# Plot the UMAP
#sc.pl.umap(
#    adata_placenta,
#    color="cell_type",
#    legend_loc="right margin",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.postintegration_on_donorid.morethan100cells.placentasubset.cell_type.allgenes.png'
#)

#sc.pl.umap(
#    adata_placenta,
#    color="donor_id",
#    legend_loc="right margin",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.postintegration_on_donorid.morethan100cells.placentasubset.donor_id.allgenes.png'
#)

#transfer cell type annotations
#adatavento = sc.read_h5ad("/staging/leuven/stg_00019/research/Kate/emseq/data/vento-tormo/adata_all_donors_trophoblast_raw_counts_in_raw_normlog_counts_in_X_for_download_UPD_20230307.h5ad")
#print(adatavento)

#convert cell barcode ids to hashed ids in census objects
#got code to transform cell barcodes to hashed ids here: https://github.com/chanzuckerberg/single-cell-curation/blob/main/cellxgene_schema_cli/cellxgene_schema/utils.py
#def hash_identifier(identifier):
    # Convert identifier to hash digest
#    hash_digest = xxh3_64_intdigest(identifier.encode())
    # Encode hash digest into Base85
#    encoded_hash = b85encode(hash_digest.to_bytes(8, "big")).decode("ascii")
#    return encoded_hash

# Sample list of cell barcodes (replace with your actual list)
#cell_barcodes = adatavento.obs_names.tolist()

# Hash each cell barcode and create DataFrame
#hashed_ids = [hash_identifier(barcode) for barcode in cell_barcodes]
#df = pd.DataFrame({'observation_joinid': hashed_ids, 'Original_ID': cell_barcodes})

# Merge the hashed IDs DataFrame with adata.obs on 'observation_joinid'
#adata.obs = adata.obs.merge(df, on='observation_joinid', how='left')
#adatavento.obs['Original_ID'] = adatavento.obs.index
#adata.obs = adata.obs.merge(adatavento.obs, on='Original_ID', how='left')

#merge vento cell type column with census cell type column
#adata.obs['final_annot_troph'] = adata.obs['final_annot_troph'].astype(str)
#adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)

#adata.obs['final_annot_troph'] = adata.obs['final_annot_troph'].replace("nan", float("NaN"))

#adata.obs['final_annot_troph'] = adata.obs['final_annot_troph'].fillna(adata.obs['cell_type'])

#create batches for correction
#adata.obs['donor_id_assay'] = adata.obs['donor_id'].astype(str) + '_' + adata.obs['assay'].astype(str)
#adata.obs['donor_id_assay'] = adata.obs['donor_id_assay'].astype(str)

#adata.obs['sample'] = adata.obs['sample'].replace("nan", float("NaN"))
#adata.obs['sample'] = adata.obs['sample'].fillna(adata.obs['donor_id_assay'])

#adata.write('/staging/leuven/stg_00019/research/Kate/emseq/data/census/tbs_vt2023trophoblasts_sleflare_censusversion2024-04-29_norm_log_scale_integrated_on_donorid_morethan100cells_nohighlyvariablegenes_allgenes_vtannotationsadded.h5ad')

#adata = sc.read_h5ad("/staging/leuven/stg_00019/research/Kate/emseq/data/census/tbs_vt2023trophoblasts_censusversion2024-04-29_norm_log_scale_integrated_on_donorid_morethan100cells_nohighlyvariablegenes_allgenes_vtannotationsadded.h5ad")

#adata_placenta = adata[adata.obs['tissue_general'] == 'placenta'].copy()
#print(adata_placenta)
#print(adata_placenta.obs['cell_type'].dtype)
#print(adata_placenta.obs['tissue_general'].dtype)
#adata_placenta.obs['cell_type'] = pd.Categorical(adata_placenta.obs['cell_type'])
#adata_placenta.obs['donor_id'] = pd.Categorical(adata_placenta.obs['donor_id'])

# Regenerate neighbors and UMAP for the placenta subset

#sc.pp.neighbors(adata_placenta, use_rep="X_scVI")
#sc.tl.umap(adata_placenta)

# Plot the UMAP
#sc.pl.umap(
#    adata_placenta,
#    color="final_annot_troph",
#    legend_loc="right margin",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.postintegration_on_sequencingbatch.morethan100cells.placentasubset.tbs.vento.final_annot_troph.allgenes.png'
#)

#sc.pl.umap(
#    adata_placenta,
#    color="donor_id",
#    legend_loc="right margin",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.postintegration_on_sequencingbatch.morethan100cells.bloodsubset.tbs.vento.donor_id.allgenes.png'
#)

#sc.pl.umap(
#    adata_placenta,
#    color="donor_id",
#    legend_loc="right margin",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.postintegration_on_sequencingbatch.morethan100cells.bloodsubset.tbs.vento.sequencing_batch.allgenes.png'
#)


# Display adata.obs with the added 'Original_ID' column
#print(adata.obs)

#print(adatavento.obs.index.head())
#print(adata_placenta.obs.index.head())
#print(adata.obs.index.head())

#sc.pl.umap(
#    adata_placenta,
#    color="assay",
#    legend_loc="on data",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.postintegration.placentasubset.assay.allgenes.png'
#)

#sc.pl.umap(
#    adata_placenta,
#    color="donor_id",
#    legend_loc="on data",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.postintegration.placentasubset.donor_id.allgenes.png'
#)

#sc.pl.umap(
#    adata_placenta,
#    color="development_stage",
#    legend_loc="on data",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.postintegration.placentasubset.dev_stage.allgenes.png'
#)

#Regenerate plots for Blood subset
#adata_blood = adata[adata.obs['tissue_general'] == 'blood'].copy()
#print(adata_blood)
#adata_blood.obs['cell_type'] = pd.Categorical(adata_blood.obs['cell_type'])

#sc.pp.neighbors(adata_blood, use_rep="X_scVI")
#sc.tl.umap(adata_blood)

#sc.pl.umap(
#    adata_blood,
#    color="cell_type",
#    legend_loc="right margin",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.postintegration_on_donorid.morethan100cells.bloodsubset.tbs.vento.sleflare.cell_type.allgenes.png'
#)

#sc.pl.umap(
#    adata_blood,
#    color="dataset_id",
#    legend_loc="right margin",
#    frameon=False,
#    legend_fontsize=5,
#    legend_fontoutline=1,
#    save='umap.postintegration_on_donorid.morethan100cells.bloodsubset.tbs.vento.sleflare.dataset_id.allgenes.png'
#)

#option to add grouped dev age category
#data = {
#    'donor_id': ['10_fetus', '11_fetus', '12_fetus', '8_fetus', '9_fetus', 'Hrv100_fetus', 'Hrv43_fetus', 'Hrv46_fetus', 'Hrv98_fetus', 'Hrv99_fetus', 'P13_fetus', 'P34_fetus', 'P14_fetus'],
#    'grouped_developmental_stage': ['4-7PCW', '4-7PCW', '8-10PCW', '4-7PCW', '8-10PCW', '4-7PCW', '4-7PCW', '8-10PCW', '8-10PCW', '8-10PCW', '8-10PCW', '11-13PCW', '8-10PCW']
#}

#write out gene expression matrix with gene names and cell type names
expr_matrix = adata.X.T
if isinstance(expr_matrix, scipy.sparse.spmatrix):
    expr_matrix = expr_matrix.todense()

#set X to the normalized batch-corrected gene expression values
#adata.layers["scVI_normalized"] = model.get_normalized_expression()
#adata.X = adata.layers["scVI_normalized"]

#read in non batch correction expression values for tbs + vt + sle dataset
#adata=sc.read_h5ad("/staging/leuven/stg_00019/research/Kate/emseq/data/census/tbs_censusversion2024-04-29_norm_log_scale_integrated_on_donorid.h5ad")

adata.var['original_var_names'] = adata.var_names
adata.var_names = adata.var['feature_id']
adata.var_names

#get average expression per tissue cell type + disease
#adata.obs['tissue_cell_type'] = adata.obs['tissue_general'].astype(str) + ' ' + adata.obs['final_annot_troph'].astype(str) + ' ' + adata.obs['dataset_id'].astype(str)
adata.obs['tissue_cell_type'] = adata.obs['tissue_general'].astype(str) + ' ' + adata.obs['cell_type'].astype(str)

cell_types = adata.obs['tissue_cell_type'].unique()

avg_expression = {}
for cell_type in cell_types:
    subset = adata[adata.obs['tissue_cell_type'] == cell_type, :]
    avg_expression[cell_type] = subset.X.mean(axis=0)

#for key, value in avg_expression.items():
#    avg_expression[key] = np.array(value).ravel()

avg_expression_df = pd.DataFrame(avg_expression, index=adata.var_names)
avg_expression_df.to_csv('/staging/leuven/stg_00019/research/Kate/emseq/data/census/tabsap_censusversion2024-04-15_normalized_log1p_scale_integrated_on_donor_id_tissue_general_cell_types_data.csv')

#get average expression per broad cell type + disease
#adata.obs['broad_cell_type'] = adata.obs['final_annot_troph'].astype(str) + ' ' + adata.obs['dataset_id'].astype(str)
#cell_types = adata.obs['broad_dataset_cell_type'].unique()
cell_types = adata.obs['cell_type'].unique()
avg_expression = {}
for cell_type in cell_types:
    subset = adata[adata.obs['cell_type'] == cell_type, :]
    avg_expression[cell_type] = subset.X.mean(axis=0)

avg_expression_df = pd.DataFrame(avg_expression, index=adata.var_names)
avg_expression_df.to_csv('/staging/leuven/stg_00019/research/Kate/emseq/data/census/tabsap_censusversion2024-04-15_normalized_log1p_scale_integrated_on_donor_id_broad_cell_types_data.csv')

