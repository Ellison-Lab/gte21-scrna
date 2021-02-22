import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sci
import random
import matplotlib.pyplot as plt

# DEBUG adata = sc.read_h5ad('out/preproc/larval-testes.celltypes.h5ad')
adata = sc.read_h5ad(snakemake.input['ad'])

expression_raw = pd.DataFrame(adata.raw.X,index = adata.raw.obs_names, columns=adata.raw.var_names.str.decode('utf-8'))
expression = pd.DataFrame(adata.X,index = adata.obs_names, columns=adata.var_names.str.decode('utf-8'))

expression_raw.reset_index().to_csv(snakemake.output['raw'], index=False)
expression.reset_index().to_csv(snakemake.output['exp'], index=False)
