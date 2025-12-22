
#%%
import pandas as pd
import scanpy as sc
import numpy as np

## ------ Parameters ---------

fname = '../../h5ad_data/Lucas_org.h5ad'
tcluster = 'anno3'
fout_prefix ='DEGbox_Lucas_withDDR'
num_cells_sampled = 200
tpval_cutoff=0.01
tlog2fc=0.5
sratio = 1.5
ddr_file = '../utils/DDR_list_final.csv'
# ddr_file= None

tgroup = 'diagnosis'
#'WT' 'Null'Mut
tgroups =  ['Control','MS'] # Have to be in this order


n_repeats = 50

output_file = './figures/%s_%s_vs_%s_boxplot2.pdf'%(fout_prefix,tgroups[1],tgroups[0])


## ------ Load data---------

adata = sc.read(fname)
if ddr_file is not None:
    ddr_list = pd.read_csv(ddr_file,header=None)[0].tolist()


# ddr_list = get_ddr_list(adata)

def get_degs_per_cluster(tadata,tgroupby,tgroups,tpval_cutoff,tlog2fc):
    """
    tgroups = [reference_data,exp_data]
    tadata 
    return up,down,all
    """
    sc.tl.rank_genes_groups(tadata,groupby=tgroupby,\
                            reference=tgroups[0],\
                            use_raw=False,\
                             method='wilcoxon')

    tpdf_up= sc.get.rank_genes_groups_df(tadata,group=tgroups[1],\
                                    pval_cutoff=tpval_cutoff,log2fc_min=tlog2fc)
    tpdf_down= sc.get.rank_genes_groups_df(tadata,group=tgroups[1],\
                                    pval_cutoff=tpval_cutoff,log2fc_max=tlog2fc*(-1))
    tpdf=pd.concat([tpdf_down,tpdf_up],ignore_index=True)
    return tpdf_up,tpdf_down,tpdf

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

#%%

tdict = {}
if ddr_file is not None:
    tdict_ddr = {}
for tkey in adata.obs[tcluster].unique().tolist():
    ## sample 1000 cell per cluster 
    sadata = adata[adata.obs[tcluster]==tkey]
    
    print("working on %s, %d cells"%(tkey,len(sadata)))
    
    if len(sadata)<num_cells_sampled*sratio:    
        print("Skip %s, only %d cells"%(tkey,len(sadata)))
        continue
    tdict[tkey] = []
    if ddr_file is not None:
        tdict_ddr[tkey] = []
    for tindex in range(n_repeats):
        
        #%% can't use fraction here, use n_obs instead
        # tdata = sc.pp.subsample(sadata,fraction=0.3,copy=True,random_state=np.random.randint(10000))
        try:
            tdata = sc.pp.subsample(sadata,n_obs=num_cells_sampled,copy=True,
                                    random_state=np.random.randint(10000))

            # print(tdata.obs[tgroup].value_counts())
            _,_,tall = get_degs_per_cluster(tdata,tgroupby=tgroup,\
                                                tgroups=tgroups,tpval_cutoff=tpval_cutoff,tlog2fc=tlog2fc)
            # print(tkey,len(tup),len(tdown),len(tall))
            ## TODO add ddr list version
            if ddr_file is not None:
                tall_ddr = intersection(tall.names.tolist(),ddr_list)
                tdict_ddr[tkey].append(len(tall_ddr))
            # tall_ddr = intersection(tall.names.tolist(),ddr_list)
            tdict[tkey].append(len(tall))
            # tdict_ddr[tkey].append(len(tall_ddr))
        except Exception as e:
            print("Error in subsampling %s, %d cells"%(tkey,len(sadata)))
            if tkey in tdict.keys():
                tdict.pop(tkey)
            if ddr_file is not None and tkey in tdict_ddr.keys():
                tdict_ddr.pop(tkey)
            print(e)
    # print(tdict)

tpdf=pd.DataFrame.from_dict(tdict)
print(tpdf.head())

#%%
# Call R function to generate boxplot
from rpy2.robjects import r, globalenv, StrVector, FloatVector # Ensure FloatVector is imported
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.rinterface import NULL as R_NULL

#%%
r_base = importr('base')
r.source("../utils/box_plot_base.R")
#%%
tpdf = tpdf.reset_index().rename(columns={'index': 'X'})
with localconverter(ro.default_converter + pandas2ri.converter):
    r_boxplot_func = r['generate_boxplot']
    r_boxplot_func(tpdf, output_file,fig_width=8, fig_height=6)



