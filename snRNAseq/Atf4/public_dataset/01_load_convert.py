#%%
"""
Load and convert Linnarsson dataset to AnnData format.

"""
#%%

import scanpy as sc


loom_path = "../data/linnarsson_adolescent.loom"

adata = sc.read_loom(loom_path)
print(adata)
#%%
import numpy as np

adata.obsm['X_umap'] = np.column_stack((adata.obs['_X'], adata.obs['_Y']))
# sc.pl.umap(adata, color=["Class"], frameon=False)
adata.obsm['X_tsne'] = np.column_stack((adata.obs['_tSNE1'], adata.obs['_tSNE2']))
# adata.obs['X_PC1'] = np.column_stack((adata.obs['_PC1'], adata.obs['_PC2']))
adata.obs['X_PCA1'] = adata.obs['_PC1']
adata.obs['X_PCA2'] = adata.obs['_PC2']
#%%
# adata.obsm['X_umap'] = adata.obsm['UMAP']
# adata.obsm['X_tsne'] = adata.obsm['TSNE']
# sc.pl.tsne(adata, color=["Label"], frameon=False)

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# adata.layers["counts"] = adata.layers["matrix"].copy()
for tkey in ['matrix', 'expected', 'pooled', 'spliced', 'unspliced']:
    adata.layers.pop(tkey, None)  # Remove specific layers if they exist
# adata.layers.pop("counts", None)  # Remove counts layer if it exists
#%%

adata.var_names_make_unique()
adata.obs_names_make_unique()
adata.write_h5ad("../h5ad_data/linnarsson_adolescent.h5ad")
#%%
import scanpy as sc
adata = sc.read_h5ad("../h5ad_data/linnarsson_adolescent.h5ad")
print(adata)



#%%

import scanpy as sc


loom_path = "../data/linnarsson_dev_all.loom"

adata = sc.read_loom(loom_path)
print(adata)
#%%
adata.obsm['X_umap'] = adata.obsm['UMAP']
adata.obsm['X_tsne'] = adata.obsm['TSNE']
sc.pl.tsne(adata, color=["Label"], frameon=False)


#%%
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# adata.layers["counts"] = adata.layers["matrix"].copy()
for tkey in ['matrix', 'expected', 'pooled', 'spliced', 'unspliced']:
    adata.layers.pop(tkey, None)  # Remove specific layers if they exist
# adata.layers.pop("counts", None)  # Remove counts layer if it exists
#%%
adata.var_names_make_unique()
adata.write_h5ad("../h5ad_data/linnarsson_dev_all.h5ad")
