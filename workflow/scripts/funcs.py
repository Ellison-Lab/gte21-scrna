import scanpy as sc
import numpy as np
import pandas as pd
import re
import anndata
import scrublet

def merge_ltr(adata_1, fbgn_f):

    # get set of TE features from the h5ad
    te_df = adata_1.var[~adata_1.var.gene_ids.str.contains('FBgn')]

    # get set of LTR TEs - 1 df each for the LTR features for the internal component
    ltr_df = te_df[te_df.gene_ids.str.contains('LTR')]
    internal_df = te_df[te_df.gene_ids.str.contains(r'[-_]I')]

    # get tuples in format: (TE shortened name with no _I or _LTR, feature name)
    ltrs = [(re.split(r'[_ | -]LTR',x)[0],x) for x in ltr_df.gene_ids]
    internal = [(re.split(r'[_ | -]I',x)[0],x) for x in internal_df.gene_ids]

    # get unique short names and feature names
    ltr_tes = set([x[0] for x in ltrs] + [x[0] for x in internal])
    ltr_tes_txs = set([x[1] for x in ltrs] + [x[1] for x in internal])

    # loop through to get dictionary of features associated with
    # each shorted TE name (only LTR family TEs)
    ltr_tes_names = {x:[] for x in ltr_tes}
    for l in ltr_tes:
        for a in ltrs:
            if a[0] == l:
                ltr_tes_names[l] += [a[1]]
        for b in internal:
            if b[0] == l:
                ltr_tes_names[l] += [b[1]]

    # DF of non-LTR TEs and LTR tes
    te_df = pd.DataFrame(index=[x for x in te_df.gene_ids.to_list() if x not in ltr_df.gene_ids.to_list() and x not in internal_df.gene_ids.to_list()] + list(ltr_tes_names.keys()))
    te_df['gene_symbol'] = te_df.index.to_list()

    # get raw count dict, then loop through the name lookup df from above.
    # on each loop I extract the two cols (_I and _LTR) for each LTR TE.
    # I then sum the counts and add back to the summed_te_df.
    m = pd.DataFrame(adata_1.X.todense(),columns=adata_1.var.index, index=adata_1.obs.index)
    summed_te_df = pd.DataFrame()
    for t in ltr_tes_names.keys():
        ix = [x for x in ltr_tes_names[t] if x != None]
        summed_te_df[t] = m[ix].sum(axis=1)

    # list of non-LTR Tes and regular genes
    non_tes = [x for x in m.columns.to_list() if not x in ltr_tes_txs]

    # get gene name/symbol data
    gene_df = pd.read_table(fbgn_f, skiprows=4, index_col='primary_FBgn#')

    # make new var df with ids from my original set, joined with full
    # set of flybase nomenclature
    var_df = pd.DataFrame(index=non_tes)
    var_df = var_df.merge(gene_df, left_index=True, right_index=True, how='left')
    var_df.columns = var_df.columns.str.replace('#','')
    var_df.columns = var_df.columns.str.replace(r'[()]','')

    var_df = var_df.append(te_df)

    var_df.loc[var_df.gene_symbol.isna(),'gene_symbol'] = var_df.loc[var_df.gene_symbol.isna()].index.to_list()

    # get rid of dup tes (because I made orig df from `non_tes`, and am appending the `te_df`
    var_df = var_df.drop_duplicates()

    # count dict from summed LTR TEs, non-LTR tes, and regular genes
    df = pd.concat([m[non_tes],summed_te_df], axis=1)

    # reorder + subset to only contain genes from my matrix - just a sanity check
    genes_to_use = df.columns.to_list()
    ix2 = [x in genes_to_use for x in var_df.index.to_list()]
    var_df = var_df.loc[ix2,:]

    # note that anndata constructor checks that indexes/columns match
    ad = anndata.AnnData(df,var=var_df)

    return ad


def call_doublets(ad, scrub_thresh, rand_state=2020):
    adata = ad.copy()

    scrub = scrublet.Scrublet(adata.X, random_state=rand_state)

    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    predicted_doublets = scrub.call_doublets(threshold=scrub_thresh)

    adata.obs['doublet_score'] = doublet_scores

    adata.obs['predicted_doublet'] = predicted_doublets

    adata = adata[~predicted_doublets]

    return adata


def qc_and_lognorm(ad, mito_f, cc_f, mito_cutoff, gene_cut, gene_min, gene_min_cells, gene_min_count):
    adata = ad.copy()

    # DEBUG mito_genes = [x.strip() for x in open("out/external/mito.tsv")]
    mito_genes = [x.strip() for x in open(mito_f)]

    # DEBUG cell_cycle_genes = pd.read_table("out/external/cell_cycle.tsv")
    cell_cycle_genes = pd.read_table(cc_f)
    g2m_genes = cell_cycle_genes[cell_cycle_genes.phase_group == 'g2m.genes'].gene_id.to_list()
    s_genes = cell_cycle_genes[cell_cycle_genes.phase_group == 's.genes'].gene_id.to_list()

    # ----
    # Basic filtering
    # ----

    sc.pp.filter_genes(adata, min_cells=gene_min_cells)
    sc.pp.filter_cells(adata, min_genes=gene_min)
    sc.pp.filter_genes(adata, min_counts=gene_min_count)

    # ----
    # Populate qc metrics
    # ----

    sc.pp.calculate_qc_metrics(adata, inplace=True)
    mg = [x for x in adata.var_names.to_list() if x in mito_genes]
    adata.obs['percent_mito'] = np.sum(adata[:, mg].X, axis=1) / np.sum(adata.X, axis=1)

    # -----
    # mito and n_gene filter
    # -----

    adata = adata[adata.obs.percent_mito < mito_cutoff, :]
    adata = adata[adata.obs.n_genes < gene_cut, :]

    g2mg = [x for x in adata.var_names.to_list() if x in g2m_genes]
    sg = [x for x in adata.var_names.to_list() if x in s_genes]

    # ----
    # Norm-transform
    # ----

    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)

    return adata


def regress_out(ad, mito_f, cc_f, max_value, vars_to_regress=None, rand_state=2020,n_jobs=1):
    adata = ad.copy()

    # DEBUG mito_genes = [x.strip() for x in open("out/external/mito.tsv")]
    mito_genes = [x.strip() for x in open(mito_f)]

    # DEBUG cell_cycle_genes = pd.read_table("out/external/cell_cycle.tsv")
    cell_cycle_genes = pd.read_table(cc_f)
    g2m_genes = cell_cycle_genes[cell_cycle_genes.phase_group == 'g2m.genes'].gene_id.to_list()
    s_genes = cell_cycle_genes[cell_cycle_genes.phase_group == 's.genes'].gene_id.to_list()

    g2mg = [x for x in adata.var_names.to_list() if x in g2m_genes]
    sg = [x for x in adata.var_names.to_list() if x in s_genes]

    adata.raw = adata

    sc.pp.scale(adata)

    # ----
    # Get scores for cell cycle
    # ----
    sc.tl.score_genes_cell_cycle(adata, s_genes=sg, g2m_genes=g2mg, random_state = rand_state)

    # ----
    # regress out
    # ----
    if vars_to_regress:
        sc.pp.regress_out(adata, vars_to_regress, n_jobs=n_jobs)
        sc.pp.scale(adata, max_value = max_value) # rescale per vignette

    return adata
