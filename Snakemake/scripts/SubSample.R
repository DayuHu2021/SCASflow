library(Seurat)
library(tidyverse)
library(dplyr)


res<- snakemake@wildcards[["res"]]
pc<- snakemake@wildcards[["pc"]]
seurat_obj <- readRDS("rawdata/sample.rds")

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

SubSample <- function(
  obj,
  n = 100,
  size = 0.8,
  npcs = 100,
  res = 1.2,
  reduction = "pca",
  assay = "SCT") {
  
  # Initialise tibble for data
  clusters <- as_tibble(Cells(obj))
  clusters <- rename(clusters, "cell" = value)
  
  # Get samples
  samples<- replicate(
      n,
      sample(Cells(obj),as.integer(length(Cells(obj)) * size),replace = FALSE),
      simplify = FALSE
    )
  # Repeated clusters
  j <- 1
  for (idx in samples) {
    message(paste0("\tClustering ", j, "..."))
    small_obj <- obj[, idx]
    small_obj <- find_clusters(
      small_obj,
      reduction = reduction,
      npcs = npcs,
      resolution = res,
      assay = assay
    )
    clusters <- dplyr::left_join(
      clusters,
      dplyr::as_tibble(Seurat::Idents(small_obj), rownames = "cell"),
      by = "cell"
    )
    j <- j + 1
  }
  return(clusters)
}


Clusters<- eval(parse(text=paste("SubSample", "(", "seurat_obj",  ",", "npcs=", pc, ",",
                                   "res=", res, ")")))
saveRDS(Clusters, file = paste0("results/SubSample","_resolution_", res, "_PC_", pc, ".rds"))
