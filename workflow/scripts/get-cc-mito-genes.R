library(Seurat)
library(tidyverse)
library(biomaRt)

txome <- snakemake@input[["txome"]]
ofl1 <- snakemake@output[["cc"]]
ofl2 <- snakemake@output[["mito"]]

# inspiration: https://support.bioconductor.org/p/96955/

message('getting mart')
# get mappings
mart1 <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

message('getting homologs')
homologs <- getBM(mart=mart1,
             attributes = c("ensembl_gene_id","external_gene_name",
                            "dmelanogaster_homolog_associated_gene_name")) %>%
  as_tibble() %>%
  mutate(ms_gene = toupper(external_gene_name))

message('importing txome')
refgenes <- rtracklayer::import(txome) %>%
  as_tibble() %>%
  mutate(dm_gene = toupper(gene_symbol))

cell_cycle <- cc.genes %>%
  enframe() %>%
  unnest(cols=value) %>%
  dplyr::rename(phase_group=name, ms_gene=value)


cc2 <- cell_cycle %>% left_join(homologs, by="ms_gene") %>%
  mutate(dmelanogaster_homolog_associated_gene_name=str_to_upper(dmelanogaster_homolog_associated_gene_name)) %>%
  left_join(refgenes, by=c(dmelanogaster_homolog_associated_gene_name="dm_gene")) %>%
  dplyr::select(phase_group, ms_gene, dmelanogaster_homolog_associated_gene_name, gene_id) %>%
  distinct() %>%
  dplyr::filter(!is.na(gene_id))


write_tsv(cc2, ofl1)

# get set of genes on mt dna
mt_tx <- refgenes %>%
  dplyr::filter(seqnames == "mitochondrion_genome") %>%
  dplyr::select(gene_id) %>%
  distinct()

write_tsv(mt_tx, ofl2,col_names = F)
