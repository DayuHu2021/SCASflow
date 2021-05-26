library(Seurat)
library(tidyverse)
library(tibble)
library(dplyr)
library(cluster)
`%>%` <- magrittr::`%>%`

sil <- readRDS(snakemake@input[[1]])
res<- snakemake@wildcards[["res"]]
pc<- snakemake@wildcards[["pc"]]

cluster_sil <- function(sil, res) {
  sil <- tibble::as_tibble(sil[, ]) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise("avg_sil" = mean(sil_width)) %>%
    tibble::add_column("res" = res)
  return(sil)
}

sil <- cluster_sil(sil,res)

saveRDS(sil, file = paste0("results/ClusterSilhouette","_resolution_", res, "_PC_", pc, ".rds"))