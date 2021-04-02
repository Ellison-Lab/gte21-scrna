from os.path import realpath
from os.path import split as pathsplit
import subprocess
import sys
import random

envvars:
    "PYTHONHASHSEED"

CUTOFF_APPROACH_ICA = config.get("CUTOFF_APPROACH_ICA", 1)
QVAL_ICA = config.get("QVAL_ICA",0.005)
COMPS_ICA = config.get("COMPS_ICA", 140)
REPS_ICA = config.get("REPS_ICA",10)
INDIV_RANDOM_SEEDS_ICA = random.sample(range(0,100000,1),k=REPS_ICA)
OUTLIER_FILT_KNN_ICA = config.get("OUTLIER_FILT_KNN_ICA",5)
OUTLIER_MAX_DIST_ICA = config.get("OUTLIER_MAX_DIST_ICA",250)
LEIDEN_RESOLUTION = config.get("LEIDEN_RESOLUTION",0.4)
ONTS = config.get("ONTS",["BP","MF","CC"])


rule get_cc_mito_genes:
    input:
        txome = custom_genome("results/custom-genome/combined.fixed.gtf"),
    output:
        cc = "results/tables/cell_cycle.tsv",
        mito = "results/tables/mito.tsv"
    conda:
        "../envs/bioc-general.yaml"
    script:
        "../scripts/get-cc-mito-genes.R"

rule export_expression:
    input:
        ad = "results/scanpy/{group}/tmp.h5ad",
    output:
        raw = "results/scanpy/{group}/lognorm-expression.csv.gz",
        exp = "results/scanpy/{group}/scaled-expression.csv.gz"
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/export-expression.py"


# rule export_umap_coords:
#     input:
#         larv_final = "results/scanpy/{group}/celltypes.h5ad",
#     output:
#         umap = "results/scanpy/{group}/umap.csv.gz",
#     conda:
#         "../envs/scanpy.yaml"
#     script:
#         "../scripts/export-coords.py"

rule export_hi_var_genes:
    input:
        adata = "results/scanpy/{group}/celltypes.h5ad"
    output:
        genes = "results/scanpy/{group}/hivar-and-tes.csv",
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/export-hivar.py"

rule get_fbgns:
    output:
        "results/fbgns/fbgns.tsv.gz"
    params:
        uri = config.get("FBGNS")
    shell:
        """
        wget -O {output} {params.uri}
        """

rule ingest:
    """
    One of several replicate integration strategies.
    Integration approach with ingest,
    """
    input:
        lambda wc: expand("results/cellranger/counts-{s}/outs/filtered_feature_bc_matrix.h5",s=config.get("SCRNA_GROUPS").get(wc.group)),
        #h5_1 = "../../data/larval-testes-01.10x.h5",
        #h5_2 = "../../data/larval-testes-02.10x.h5",
        #h5_3 = "../../data/larval-testes-03.10x.h5",
        fbgn_f = rules.get_fbgns.output, #determine_resource(config.get("FBGNS",None)),
        mito_f = rules.get_cc_mito_genes.output.mito,
        cc_f = rules.get_cc_mito_genes.output.cc,
    resources:
        time=20,
        mem=12000,
        cpus=8
    threads:
        8
    output:
        h5ad = temp("results/scanpy/{group}/tmp.h5ad"),
        clusters = 'results/scanpy/{group}/figs/communities.pdf'
    params:
        scrub_thresh = lambda wc: [float(pep.get_sample(x).doublet_threshold[0]) for x in config.get("SCRNA_GROUPS").get(wc.group)],
        vars = config.get("REGRESS_OUT",None),
        max_value = config.get("MAX_SCALE"),
        mito_cutoff = config.get('MITO_CUTOFF'),
        gene_cut = config.get('MAX_GENES'),
        seed = config.get('RANDOM_SEED'),
        gene_min = config.get('GENE_MIN'), #250
        gene_min_count = config.get('GENE_MIN_COUNT'), # 1
        gene_min_cells = config.get('GENE_MIN_CELLS'), # 3
        n_neighbors = config.get('N_NEIGHBORS'),
        n_pcs_for_neighbors = config.get('N_PCS_FOR_NEIGHBORS'),
        leiden_resolution = LEIDEN_RESOLUTION,
        umap_min_dist = config.get('UMAP_MIN_DIST'),
        umap_spread = config.get('UMAP_SPREAD'),
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/preproc.py"

rule export_raw_loom:
    input:
        "results/scanpy/{group}/tmp.h5ad"
    output:
        loom = temp("results/scanpy/{group}/raw.loom"),
        anno = directory('results/scanpy/{group}/anno/')
    params:
        type = "raw",
    resources:
        time=20,
        mem=12000,
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/export_loom.py"

rule get_anno:
    input:
        'results/scanpy/{group}/anno/'
    output:
        var='results/scanpy/{group}/var.csv',
        obs='results/scanpy/{group}/obs.csv',
        obsm='results/scanpy/{group}/obsm.csv',
    params:
        dir = lambda wc: 'results/scanpy/{g}/anno/'.format(g=wc.group)
    shell:
        """
        cp {input}/var.csv {output.var} &&
        cp {input}/obs.csv {output.obs} &&
        cp {input}/obsm.csv {output.obsm}
        """

rule get_raw_umis:
    input:
        h5 = lambda wc: expand("results/cellranger/counts-{s}/outs/filtered_feature_bc_matrix.h5",s=config.get("SCRNA_GROUPS").get(wc.group)),
        obs = 'results/scanpy/{group}/obs.csv',
        var = 'results/scanpy/{group}/var.csv',
    output:
        csv="results/scanpy/{group}/umis.csv.gz"
    resources:
        time=20,
        mem=20000,
        cpus=2
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/get_raw_umis.R"


rule larval_testis_community_to_celltypes:
    """
    Assign likely cell types from the larval clusters with 'garnett' R
    package.
    """
    input:
        larval = "results/scanpy/{group}/raw.loom",
        markers = config.get('LARVAL_TESTIS_MARKERS'),
        csvs = 'results/scanpy/{group}/anno/',
    threads:
        24
    resources:
        time=60,
        mem=12000,
        cpus=24
    output:
        ofl = "results/scanpy/{group}/clusters-to-celltypes.csv",
        #g_table = "results/scanpy/{group}/figs/community-to-celltype-table.pdf",
        #g_garnett_results = "results/scanpy/{group}/figs/community-to-celltype-barchart.pdf",
        #g_marker_check = "results/scanpy/{group}/figs/celltype-marker-check.pdf",
    #conda:
    #    "../envs/monocle3.yaml"
    shell:
        "Rscript workflow/scripts/clusters-to-celltypes.R {input.larval} {input.markers} {input.csvs} {output.ofl} {threads}"

rule add_larval_celltype:
    """
    Add larval cell types to testis anndata.
    """
    input:
        ad = "results/scanpy/{group}/tmp.h5ad",
        ct = "results/scanpy/{group}/clusters-to-celltypes.csv"
    output:
        ad = "results/scanpy/{group}/celltypes.h5ad",
        clusters = 'results/scanpy/{group}/figs/clusters.pdf'
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/add-larval-celltypes.py"
