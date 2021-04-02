library(Seurat)
library(tidyverse)

#h5s_fl <- Sys.glob("~/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/cellranger/counts-larval-testes-*/outs/filtered_feature_bc_matrix.h5")
#obs_fl <- "~/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/scanpy/larval-w1118-testes/obs.csv"
#var_fl <-"~/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/scanpy/larval-w1118-testes/var.csv" 

h5s_fl <- snakemake@input[["h5"]]
obs_fl <- snakemake@input[["obs"]]
var_fl <- snakemake@input[["var"]]

obs <- read_csv(obs_fl) %>%
  mutate_if(is.character, ~str_remove(.,"^b(?=')")) %>%
  mutate_if(is.character, ~str_remove_all(.,"'"))

var <- read_csv(var_fl) %>%
  mutate_if(is.character, ~str_remove(.,"^b(?=')")) %>%
  mutate_if(is.character, ~str_remove_all(.,"'"))

grp <- ifelse(str_detect(snakemake@wildcards[["group"]],"larval"),"larval","batch")

# get list of filtered cells in each rep.
filtered_cells <- obs %>%
  dplyr::select(batch, X1) %>%
  mutate(X1=str_remove(X1,"(?<=-1)-.+")) %>%
  split(.,.$batch) %>%
  map(pull, X1)

# read in raw umis for post filtered cells
dat <- h5s_fl %>% set_names(.,paste0(grp,str_extract(.,"(?<=testes-)\\d+(?=\\/)"))) %>%
  map(Seurat::Read10X_h5) %>%
  imap(~{
    .x[which(rownames(.x) %in% var$X1),which(colnames(.x) %in% filtered_cells[[.y]])]
    #.x[which(!str_detect(rownames(.x),"FBgn")),which(colnames(.x) %in% filtered_cells[[.y]])]
  }) %>%
  map(Matrix::as.matrix) %>%
  map(as_tibble, rownames = "gene") %>%
  map_df(~gather(.,cell,UMIs,-gene), .id = "batch") %>%
  mutate(cell = paste(cell,batch, sep="-"))

stopifnot(all(obs$X1 %in% dat$cell))

write_csv(dat,snakemake@output[["csv"]])

# quick check
#dat %>%
#  filter(!str_detect(gene,"FBgn")) %>%
#  left_join(dplyr::select(obs, X1, clusters), by=c(cell='X1')) %>%
#  left_join(read_csv("~/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/scanpy/larval-w1118-testes/clusters-to-celltypes.csv",col_types = "ccddc" )) %>%
#  group_by(clusters2) %>%
#  summarise(UMIs = sum(UMIs)/n()) %>%
#  ggplot(aes(clusters2,UMIs)) +
#  geom_col() +
#  coord_flip()
