#轮廓指数
library(Seurat)
library(tidyverse)
library(tibble)
library(dplyr)
library(cluster)


assay <- "SCT"
reduction <- "pca"
res<- snakemake@wildcards[["res"]]
pc<- snakemake@wildcards[["pc"]]
mtchs <- readRDS(snakemake@input[[1]])
seurat_obj <- readRDS("rawdata/pbmc3k.rds")

find_clusters <- function(
                          obj,
                          reduction = "pca",
                          npcs = 100,
                          assay = "SCT",
                          features = NULL,
                          resolution = 0.8,
                          verbose = FALSE) {
  obj <- Seurat::FindNeighbors(
    obj,
    reduction = reduction,
    dims = 1:npcs,
    assay = assay,    ###changed from original repository
    features = features,
    verbose = verbose,
    graph.name = paste(reduction, assay, sep = ".")
  )
  obj <- Seurat::FindClusters(
    obj,
    resolution = resolution,
    graph.name = paste(reduction, assay, sep = "."),
    verbose = verbose
  )
  return(obj)
}


obj<- eval(parse(text=paste("find_clusters", "(", "seurat_obj",  ",", "npcs=", pc, ",",
                                   "res=", res, ")")))
clusters <-obj[[glue::glue("{reduction}.{assay}_res.{res}")]]

sil <- cluster::silhouette(
  x = as.numeric(as.character(unlist(clusters))),
  dmatrix = (1 - as.matrix(mtchs))
)

saveRDS(sil, file = paste0("results/CellSilhouette","_resolution_", res, "_PC_", pc, ".rds"))