#Cocluster
library(Seurat)
library(tidyverse)
library(dplyr)


res<- snakemake@wildcards[["res"]]
pc<- snakemake@wildcards[["pc"]]
results <- readRDS(snakemake@input[[1]])

find_matches <- function(col, df) {
  mtchs <- outer(df[[col]], df[[col]], "==")
  mtchs[is.na(mtchs)] <- 1i
  return(mtchs)
}

percent_match <- function(x, n = 100) {
  return(Re(x) / (n - Im(x)))
}

columns <- colnames(select(results, -cell))
mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
i <- 1 # Counter
for (col in columns) {
  mtchs <- Reduce("+", list(
    mtchs,
    find_matches(col, df = results)
  ))
  i <- i + 1
}

mtchs <- mutate_all(
  as_tibble(mtchs),
  function(x) if_else(Re(x) > 0, percent_match(x), 0)
)

saveRDS(mtchs, file = paste0("results/CoClustering","_resolution_", res, "_PC_", pc, ".rds"))