import anndata
import scanpy as sc

#ad = sc.read_h5ad("wf/preproc-testis/out/larval-testes.ingest-integrated.h5ad")
ad = sc.read_h5ad(snakemake.input[0])

ad.write_csvs(snakemake.output['anno'])

if snakemake.params['type'] == 'raw':
    ad = ad.raw.to_adata()

ad.write_loom(snakemake.output['loom'])
