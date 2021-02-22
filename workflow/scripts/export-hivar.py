import scanpy as sc
import pandas as pd

ad = sc.read(snakemake.input['adata'])

hi_var = list(ad.var_names[ad.var.loc[:,ad.var.columns.str.contains('highly_variable')].any(axis=1)].str.decode('utf-8'))

tes = list(ad.var_names[~ad.var_names.str.decode('utf-8').str.contains("FBgn")].str.decode('utf-8'))

genes = hi_var + tes

df = pd.DataFrame({'genes': genes})

df.to_csv(snakemake.output['genes'],index=False)
