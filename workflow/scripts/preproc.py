import scanpy
import numpy as np
import pandas as pd
import re
import anndata
import sys
import scrublet
import matplotlib.pyplot as plt
import random
import itertools
import functools

random.seed(snakemake.params['seed'])
np.random.seed(snakemake.params['seed'])

rand_state=snakemake.params['seed']

THREADS = snakemake.threads

sys.path.append('scripts')
from funcs import *

# h5s = ["/home/mlawlor/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/cellranger/counts-larval-testes-01/outs/filtered_feature_bc_matrix.h5","/home/mlawlor/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/cellranger/counts-larval-testes-02/outs/filtered_feature_bc_matrix.h5","/home/mlawlor/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/cellranger/counts-larval-testes-03/outs/filtered_feature_bc_matrix.h5"]
# scrub_thresh = [0.35,0.3,0.35]
# fbgn_f = "/home/mlawlor/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/ftp.flybase.org/releases/FB2020_01/precomputed_files/genes/fbgn_annotation_ID_fb_2020_01.tsv.gz"
# mito_f = "/home/mlawlor/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/tables/mito.tsv"
# cc_f = "/home/mlawlor/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/tables/cell_cycle.tsv"
# vars_to_regress = ['S_score','G2M_score'] #,'percent_mito']
# mito_cutoff= 0.05
# max_value = 100
# gene_cut = 5000

h5s = list(itertools.chain(*snakemake.input['h5s']))

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

ads = [scanpy.read_10x_h5(h5) for h5 in h5s]

ads = [merge_ltr(a,fbgn_f) for a in ads]

ads = [call_doublets(a,st) for a,st in zip(ads,scrub_thresh)]

ads = [qc_and_lognorm(a, mito_f, cc_f, mito_cutoff=mito_cutoff, gene_cut=gene_cut, gene_min=gene_min, gene_min_count=gene_min_count, gene_min_cells=gene_min_cells) for a in ads]

ads = [regress_out(a, mito_f, cc_f, max_value, vars_to_regress, n_jobs=THREADS) for a in ads]

var_names = list(functools.reduce(lambda x,y: x.intersection(y), [z.var_names for z in ads]))

ads = [a[:,var_names].copy() for a in ads]

for a in ads:
  a.var.loc[~a.var_names.str.match("FBgn"),'highly_variable'] = False
  sc.tl.pca(a, n_comps= n_pcs_for_neighbors*2, use_highly_variable=True)
  sc.pp.neighbors(a, n_neighbors=n_neighbors, n_pcs=n_pcs_for_neighbors)

sc.tl.leiden(ads[-1], key_added='clusters', resolution=leiden_resolution, random_state=rand_state)

sc.tl.umap(ads[-1], random_state=rand_state,  min_dist=umap_min_dist, spread=umap_spread)

if len(ads) > 1:
    for a in ads[:-1]:
        sc.tl.ingest(a, ads[-1], obs=['clusters'])
    adata = ads[-1].concatenate(ads[:-1])
else:
    adata = ads[0]

adata.write(snakemake.output['h5ad'])

ax = sc.pl.umap(adata, color=['clusters','Copia2'], legend_loc='on data', gene_symbols='gene_symbol', use_raw=True, show=False)

# DEBUG = plt.savefig("test.png")
plt.savefig(snakemake.output["clusters"])
