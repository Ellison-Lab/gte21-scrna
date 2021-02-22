.libPaths(
	c("/cache/home/mal456/miniconda3/lib/R/library",
	"/cache/home/mal456/R/4.0.2/lib64/R/library")
)

library(LoomExperiment)
library(monocle3)
library(tidyverse)
library(org.Dm.eg.db)
library(garnett)

args = commandArgs(trailingOnly=TRUE)
fl = args[1]
markers_fl = args[2]
obs_fl = paste0(args[3],'/obs.csv')
var_fl = paste0(args[3],'/var.csv')
ofl = args[4]
cores = args[5]

sessionInfo()
set.seed(2020)
#set.seed(snakemake@params[['seed']])
# fl <- "~/amarel-scratch/TestisTEs2021/results/scanpy/larval-w1118-testes/raw.loom"; obs_fl <- "~/amarel-scratch/TestisTEs2021/results/scanpy/larval-w1118-testes/anno/obs.csv"; markers_fl <- "~/amarel-scratch/TestisTEs2021/resources/larval-testis-markers.txt"; var_fl <- "~/amarel-scratch/TestisTEs2021/results/scanpy/larval-w1118-testes/anno/var.csv";cores <- 8
#fl <- snakemake@input[['larval']]
#obs_fl <- paste0(snakemake@input[['csvs']],'/obs.csv')
#var_fl <- paste0(snakemake@input[['csvs']],'/var.csv')
#markers_fl <- snakemake@input[['markers']]
#ofl <- snakemake@output[["ofl"]]
#cores <- snakemake@threads[[1]]

obs <- read_csv(obs_fl) %>%
	mutate_if(is.character,~str_extract(.,regex("(?<=b').+(?=')"))) %>%
	column_to_rownames("X1")

var <- read_csv(var_fl) %>%
	mutate_if(is.character,~str_extract(.,regex("(?<=b').+(?=')"))) %>%
	mutate(gene_symbol = ifelse(is.na(gene_symbol),annotation_ID,gene_symbol)) %>%
  mutate(gene_short_name = `gene_symbol`) %>%
  column_to_rownames('gene_symbol')

scle <- import(fl, type="SingleCellLoomExperiment")

X <- assay(scle,'matrix')
X <- as.matrix(X)

colnames(X) <- rownames(obs)
rownames(X) <- rownames(var)

#pd <- new("AnnotatedDataFrame", data = obs)
#fd <- new("AnnotatedDataFrame", data = var)

message('make cds')
message(class(X))

#cds <- newCellDataSet(as(X, "dgCMatrix"),
#                             phenoData = pd,
#                             featureData = fd)
cds <- new_cell_data_set(X,
                             cell_metadata = obs,
                             gene_metadata = var)

message('size factors')
cds <- estimate_size_factors(cds)

message('marker check')
marker_check <- check_markers(cds, markers_fl,
                              db=org.Dm.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

g_marker_check <- plot_markers(marker_check)

classifier <- train_cell_classifier(cds = cds,cores=cores,
                                         marker_file = markers_fl,
                                         db=org.Dm.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 500,
                                         marker_file_gene_id_type = "SYMBOL")

cds <- classify_cells(cds, classifier,
                           db = org.Dm.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "SYMBOL")

saveRDS(cds, '~/test.rds')

message('printing result')
cds@colData


res <- cds@colData %>%
  as_tibble(rownames = 'X1') %>%
  dplyr::select(X1,clusters, garnett_cluster,cell_type, cluster_ext_type) %>%
  mutate_at(vars('clusters','garnett_cluster'), as.character)

# output final summary table
res2 <- res %>%
  group_by(clusters,cell_type) %>%
  summarize(n=n()) %>%
  group_by(clusters) %>%
  mutate(pct = n/sum(n)) %>%
  filter(cell_type != "Unknown") %>%
  top_n(1,n) %>%
  ungroup() %>%
  mutate(clusters2 = paste(clusters,cell_type,sep="/"))

message('printing res2')
message(res2)

write_csv(res2, ofl)
