#%%
import scanpy as sc
import glob
from desc_pytorch import train
import os 
import sys
sys.path.append('/Users/xuz3/single_cell_python/Pipeline_ref')
from map_single_file import map_adata
#%%

save_intermediate = False

# version = 'v01'
# n_top_genes = 2048
# ae_dims = [n_top_genes, 256, 64]
# louvain_resolution = [0.6,1.0,1.2]
# tneighbor = 80

# version = 'v02'
# n_top_genes = 2048
# ae_dims = [n_top_genes, 256, 64]
# louvain_resolution = [0.6,1.0,1.2]
# tneighbor = 32

version = 'v03'
n_top_genes = 1024
ae_dims = [n_top_genes, 256, 64]
louvain_resolution = [0.6,1.0,1.2]
tneighbor = 20


root_dir = '../Atf4/data/'
flist= glob.glob(root_dir + '*/*_filtered.h5')


tdict= {'72_97':'Con_Cux2_WT_WT','75_99':'Con_Cux2_WT_WT' ,'74_98':'Con_Cux2_WT_WT',
        '25_88':'Con_Cux2_Cre_Cre' ,'19_84': 'Con_Cux2_Cre_Cre' ,'11_81': 'Con_Cux2_Cre_Cre',
        '50_94': 'Cux2_Cre_WT_Atf4_ff','9_93': 'Cux2_Cre_WT_Atf4_ff','51_95': 'Cux2_Cre_WT_Atf4_ff',
        '6_4': 'Cux2_Cre_Cre_Atf4_ff','10_92': 'Cux2_Cre_Cre_Atf4_ff','8_90': 'Cux2_Cre_Cre_Atf4_ff'}


#%% # load data
if os.path.exists(root_dir + 'processed/11_preprocess.h5ad') is False:
    sc_list = []

    for tf in flist:
        tdata = sc.read_10x_h5(tf)
        tsample =  '_'.join(tf.split('/')[-1].split('_')[:2])
        tgroup = tdict[tsample]
        tdata.obs['Condition'] = tgroup
        tdata.obs['Channel'] = tsample
        # Add channel name to obs_names and make them unique
        tdata.obs_names = [f"{tsample}_{x}" for x in tdata.obs_names]
        # Make var_names unique
        tdata.var_names_make_unique()
        tdata.obs_names_make_unique()
        sc_list.append(tdata)
    #%% # concatenate data
    adata = sc.concat(sc_list)
    adata.obs_names_make_unique()

    print(adata)
    if save_intermediate:   
        adata.write(root_dir + 'processed/10_load.h5ad')


    #%% QC and filtering

    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("Mt-")  # "MT-" for human, "Mt-" for mouse
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(( "Rps", "Rpl"))  # "RPS" and "RPL" for human, "Rps" and "Rpl" for mouse
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^Hb[^(p)]")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True)

    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata = adata[adata.obs['pct_counts_mt'] < 10].copy()

    sc.pp.scrublet(adata, batch_key="Channel")

    adata.layers["counts"] = adata.X.copy()

    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data:
    sc.pp.log1p(adata)
    #if save_intermediate:       
    map_adata(adata)
    adata.write(root_dir+"processed/11_preprocess.h5ad")


#%% # DESC clustering
#check if adata exists
try:
    adata
except NameError:
    # adata = None    
    adata = sc.read(root_dir+"processed/11_preprocess.h5ad")


adata_org = adata.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes,\
                                          flavor="seurat_v3", subset=True,layer="counts")
sc.pp.scale(adata, max_value=3)

adata_desc=train(data=adata,dims=ae_dims,n_neighbors=tneighbor,\
                 louvain_resolution = louvain_resolution,\
                    do_tsne=False,do_umap=True,\
                        save_dir='results/%s'%(version),)

# adata_org = sc.read("/Users/xuz3/single_cell_python/Caffeine/data/11_preprocess.h5ad")

for tres in louvain_resolution:
    res = str(tres).replace('.', '_')
    umap_key = f'umap_{res}'
    desc_key = f'desc_{res}'
    adata_org.obsm[f'X_umap_{res}'] = adata_desc.obsm[umap_key]
    adata_org.obs[desc_key] = adata_desc.obs[desc_key]
    
    adata_org.obsm['X_umap'] = adata_desc.obsm[umap_key]
    sc.pl.umap(adata_org, color=[desc_key], frameon=False,save='_%s_%s.pdf'%(desc_key,version),show=False)
#%%

adata_org.write(root_dir+"processed/12_desc_%s.h5ad"%version)

# %%  Where and how to use the cell type mapping code ?
if __name__ == '__main__':
    map_adata(adata_org)
    adata_org.write(root_dir+"processed/12_desc_%s.h5ad"%version)