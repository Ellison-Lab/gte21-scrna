import scanpy as sc
import anndata
import pandas as pd
import matplotlib.pyplot as plt

#adata = sc.read_h5ad('/home/mlawlor/work/TestisTpn/out/preproc/larval-testes.ingest-integrated.h5ad')
adata = sc.read_h5ad(snakemake.input['ad'])

#celltypes = pd.read_csv("exploratory/200711_dev_cellassign/anno/community-calls.csv")
celltypes = pd.read_csv(snakemake.input['ct'])

celltypes = {str(celltypes.iloc[x,0]):celltypes.iloc[x,4] for x in celltypes.index}

adata.obs = adata.obs.replace({'clusters':celltypes}, inplace=False)

ax = sc.pl.umap(adata, color=['clusters'], legend_loc='on data', gene_symbols='gene_symbol', use_raw=True, show=False)

# DEBUG = plt.savefig("test.png")
plt.savefig(snakemake.output["clusters"])

adata.write(snakemake.output['ad'])
