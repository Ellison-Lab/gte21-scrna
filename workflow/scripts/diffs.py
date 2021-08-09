import scanpy
import pandas

#h5ad_fl = "results/scanpy/larval-w1118-testes/celltypes.h5ad"
h5ad_fl = snakemake.input['ad']

ad = scanpy.read_h5ad(h5ad_fl)

scanpy.tl.rank_genes_groups(ad, "clusters")

unique_clusters = ad.obs.clusters.cat.categories.to_list()

res_dict_list = [scanpy.get.rank_genes_groups_df(ad, group=str(x)) for x in unique_clusters]

res = pandas.concat(res_dict_list, keys=unique_clusters, names=['id','index'])

res = res.reset_index(level=0)

res['id'] = res['id'].str.decode('utf-8')
res['names'] = res['names'].str.decode('utf-8')

res.to_csv(snakemake.output['ad'])
