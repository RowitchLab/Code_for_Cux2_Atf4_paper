#%%
import scanpy as sc
import pandas as pd

#%%
# adolescent_file = "/Users/xuz3/single_cell_python/h5ad_data/linnarsson_adolescent.h5ad"
dev_file = "/Users/xuz3/single_cell_python/h5ad_data/linnarsson_dev_all.h5ad"

# adata = sc.read_h5ad(adolescent_file)
adata = sc.read_h5ad(dev_file)

#%%

adata = adata[adata.obs.Age.isin(["e17.0", "e18.0","e17.5"])]
print(adata.obs.Age.value_counts())

adata = adata[adata.obs['Subclass'].isin(["Cerebelllum glutamatergic",\
                                          "Cerebelllum GABAergic",\
                                        "Forebrain glutamatergic",\
                                        "Forebrain GABAergic",\
                                        "Midbrain glutamatergic",\
                                        "Midbrain GABAergic",\
                                        "Cortical or hippocampal glutamatergic"])]
print(adata.obs['Subclass'].value_counts())

#%%

adolescent_file = "/Users/xuz3/single_cell_python/h5ad_data/linnarsson_adolescent.h5ad"

adata = sc.read_h5ad(adolescent_file)
print(adata)



adata = adata[adata.obs.Age.isin(["p21", "p22","p23", "p25","p26","p27","p28", "p29","p60"])]
print(adata.obs.Age.value_counts())

adata = adata[adata.obs['TaxonomyRank4'].isin(['Telencephalon inhibitory interneurons','Telencephalon projecting excitatory neurons', ])]
print(adata.obs['TaxonomyRank4'].value_counts())

print(adata.X.max())


#%%
import os
os.environ["R_HOME"]="/Library/Frameworks/R.framework/Resources"
#export R_LD_LIBRARY_PATH="/Library/Frameworks/R.framework/Resources/lib"
os.environ["R_LD_LIBRARY_PATH"]="/Library/Frameworks/R.framework/Resources/lib"
os.environ["RSCRIPT_PATH"]="/usr/local/bin/Rscript"


import anndata
import scipy.sparse
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import StrVector,FloatMatrix, ListVector,FloatVector



def convert_h5ad_to_rds(adata, rds_path: str, save_coords=False,normalize=False):
    """
    AnnData -> Seurat .rds
    1. 优先用 raw counts (adata.raw 或 layers['counts'])
    2. 携带全部 obs 矩阵  &  obsm 中的嵌入坐标（UMAP/t-SNE等）
    3. 基因名零重复、无下划线
    """
    # ---------- 1.  R 包 ----------
    base   = importr('base')
    seurat = importr('Seurat')
    Matrix = importr('Matrix')

    # ---------- 2.  读 h5ad ----------
    # adata = anndata.read_h5ad(h5ad_path)
    print(f'original shape: {adata.shape}')

    # ---------- 3.  取 raw counts ----------
    if adata.raw is not None:
        # 用 .raw 里的 counts
        counts = adata.raw.X.copy()
        var_names = adata.raw.var_names
        print('Use adata.raw.X as raw counts')
    elif 'counts' in adata.layers:
        counts = adata.layers['counts'].copy()
        var_names = adata.var_names
        print('Use adata.layers["counts"] as raw counts')
    else:
        counts = adata.X.copy()
        var_names = adata.var_names
        print('WARNING: no raw counts found -> fall back to adata.X')

    # 保证稀疏
    if not scipy.sparse.issparse(counts):
        counts = scipy.sparse.csr_matrix(counts)

    # ---------- 4.  基因名去重 ----------
    var_names = pd.Index(var_names).astype(str)
    if var_names.duplicated().any():
        var_names = var_names.make_unique()   # pandas 自带
    counts = counts[:, ~var_names.duplicated()].copy()  # 丢掉重复列
    adata = adata[:, var_names].copy()        # 让 adata 与 counts 对齐
    print(f'after dup-remove: {adata.shape}')

    # ---------- 5.  构造 R sparseMatrix ----------
    counts_coo = counts.T.tocoo()               # Seurat 要 基因×细胞
    genes = var_names.tolist()
    cells = adata.obs_names.astype(str).tolist()

    r_counts = Matrix.sparseMatrix(
        i=ro.IntVector(counts_coo.row + 1),
        j=ro.IntVector(counts_coo.col + 1),
        x=ro.FloatVector(counts_coo.data),
        dims=ro.IntVector(counts_coo.shape),
        dimnames=ro.r.list(genes=StrVector(genes), cells=StrVector(cells))
    )

    # ---------- 6.  清洗行名（下划线->dash，零重复） ----------
    ro.r.assign('counts', r_counts)
    ro.r('''
    rn <- (rownames(counts))
    rn[rn == ""] <- NA
    rn <- gsub("_", "-", rn)
    rn[is.na(rn)] <- paste0("GENE-", seq_len(sum(is.na(rn))))
    rownames(counts) <- make.unique(rn)
    ''')
    r_counts_clean = ro.r('counts')
    ro.r('stopifnot(anyDuplicated(rownames(counts)) == 0)')

    # ---------- 7.  meta.data ----------
    with localconverter(default_converter + pandas2ri.converter):
        meta_r = pandas2ri.py2rpy(adata.obs)

    # ---------- 8.  创建 Seurat 对象 ----------
    seurat_obj = ro.r['CreateSeuratObject'](counts=r_counts_clean, meta_data=meta_r)

    if save_coords:
        # ---------- 9.  把 obsm 嵌入坐标写进去 ----------
        # 支持的维度名
        # 如果 reductions 不存在，先初始化为空 list()
        # 初始化 reductions（仅第一次循环前即可）
        ro.r.assign('seurat_obj', seurat_obj)
        ro.r('if (!"reductions" %in% slotNames(seurat_obj)) \
                seurat_obj@reductions <- list()')
        embed_map = {'X_umap': 'umap', 'X_tsne': 'tsne'}
        for key, seurat_name in embed_map.items():
            if key in adata.obsm:
                coord = adata.obsm[key][:, :2]
                cells = adata.obs_names.astype(str).tolist()
                coord_r = ro.r.matrix(FloatVector(coord.T.ravel()),
                                    nrow=coord.shape[0],
                                    ncol=2,
                                    dimnames=ro.r.list(StrVector(cells),
                                                        StrVector(['Dim1', 'Dim2'])))

                # 导变量
                ro.r.assign('seurat_obj', seurat_obj)
                ro.r.assign('coord_r', coord_r)
                ro.r.assign('name', seurat_name)

                # 官方写法：创建 + 直接写 slot
                ro.r('seurat_obj@reductions[[name]] <- CreateDimReducObject(\
                        embeddings = coord_r,\
                        key        = paste0(name, "_"),\
                        assay      = DefaultAssay(seurat_obj))')
                print(f'  -- added {seurat_name} reduction')

                seurat_obj = ro.r('seurat_obj')
    
    if normalize:
        ro.r.assign('seurat_obj', seurat_obj)
        ro.r('seurat_obj<-NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000,verbose = FALSE)')
        # # ---------- 9.5  Normalization ----------
        #seurat_obj = ro.r['NormalizeData'](seurat_obj, normalization_method = "LogNormalize", scale_factor = 10000)
        print('✅ NormalizeData done (log-normalized)')
        seurat_obj = ro.r('seurat_obj')
        print(ro.r('print(seurat_obj)'))

    # ---------- 10.  保存 ----------
    base.saveRDS(seurat_obj, file=rds_path)
    print(f'✅ RDS saved -> {rds_path}')

#%%
# rds_path = "linnarsson_adolescent_filtered.rds"
rds_path = "linnarsson_dev_all_filtered.rds"
convert_h5ad_to_rds(adata, rds_path, save_coords= True,normalize=True)