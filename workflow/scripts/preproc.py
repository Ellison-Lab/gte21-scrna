import scanpy
import numpy as np
import pandas as pd
import re
import anndata
import sys
import scrublet
import matplotlib.pyplot as plt
import random

random.seed(snakemake.params['seed'])
np.random.seed(snakemake.params['seed'])

rand_state=snakemake.params['seed']

THREADS = snakemake.threads
print('AAAA')
print(THREADS)
sys.path.append('scripts')
from funcs import *

# h5_1 = "../../data/larval-testes-01.10x.h5"
# h5_2 = "../../data/larval-testes-02.10x.h5"
# h5_3 = "../../data/larval-testes-03.10x.h5"
# scrub_thresh = [0.35,0.3,0.35]
# fbgn_f = "/home/mlawlor/work/TestisTpn/workflow/all-testis/ftp.flybase.org/releases/FB2020_01/precomputed_files/genes/fbgn_annotation_ID_fb_2020_01.tsv.gz"
# mito_f = "/home/mlawlor/work/TestisTpn/results2/external/mito.tsv"
# cc_f = "/home/mlawlor/work/TestisTpn/results2/external/cell_cycle.tsv"
# vars_to_regress = ['S_score','G2M_score'] #,'percent_mito']
# mito_cutoff= 0.05
# max_value = 100
# gene_cut = 5000

h5_1 = snakemake.input[0]
h5_2 = snakemake.input[1]
h5_3 = snakemake.input[2]

scrub_thresh = snakemake.params['scrub_thresh']
max_value = snakemake.params['max_value']

mito_cutoff = snakemake.params['mito_cutoff'] # 0.05
gene_cut = snakemake.params['gene_cut'] # 5000
mito_f = snakemake.input['mito_f']
gene_min = snakemake.params['gene_min'] #250
gene_min_count = snakemake.params['gene_min_count'] # 1
gene_min_cells = snakemake.params['gene_min_cells'] # 3
leiden_resolution = snakemake.params['leiden_resolution'] # 3
n_neighbors = snakemake.params['n_neighbors']
n_pcs_for_neighbors = snakemake.params['n_pcs_for_neighbors']
umap_min_dist = snakemake.params['umap_min_dist']
umap_spread = snakemake.params['umap_spread']

cc_f = snakemake.input['cc_f']
vars_to_regress = snakemake.params['vars']

fbgn_f = snakemake.input["fbgn_f"]


ads = [scanpy.read_10x_h5(h5) for h5 in [h5_1, h5_2, h5_3]]

ads = [merge_ltr(a,fbgn_f) for a in ads]

ads = [call_doublets(a,st) for a,st in zip(ads,scrub_thresh)]

ads = [qc_and_lognorm(a, mito_f, cc_f, mito_cutoff=mito_cutoff, gene_cut=gene_cut, gene_min=gene_min, gene_min_count=gene_min_count, gene_min_cells=gene_min_cells) for a in ads]

ads = [regress_out(a, mito_f, cc_f, max_value, vars_to_regress, n_jobs=THREADS) for a in ads]

var_names = ads[0].var_names.intersection(ads[1].var_names).intersection(ads[2].var_names)

adata_1 = ads[0][:,var_names]
adata_2 = ads[1][:,var_names]
adata_3 = ads[2][:,var_names]

for a in [adata_1, adata_2, adata_3]:
  a.var.loc[~a.var_names.str.match("FBgn"),'highly_variable'] = False
  sc.tl.pca(a, n_comps= n_pcs_for_neighbors*2, use_highly_variable=True)
  sc.pp.neighbors(a, n_neighbors=n_neighbors, n_pcs=n_pcs_for_neighbors)

sc.tl.leiden(adata_3, key_added='clusters', resolution=leiden_resolution, random_state=rand_state)

sc.tl.umap(adata_3, random_state=rand_state,  min_dist=umap_min_dist, spread=umap_spread)

for a in [adata_1,adata_2,]:
  sc.tl.ingest(a, adata_3, obs=['clusters'])

adata = adata_3.concatenate(adata_1,adata_2, batch_categories=['larval03', 'larval01','larval02'])

adata.write(snakemake.output['h5ad'])

ax = sc.pl.umap(adata, color=['clusters','Copia2'], legend_loc='on data', gene_symbols='gene_symbol', use_raw=True, show=False)

# DEBUG = plt.savefig("test.png")
plt.savefig(snakemake.output["clusters"])
