

#the list of function
DrawFile <-function(results_path,res,pc){
  scores <- purrr::map(
    paste0(results_path, "ClusterSilhouette_resolution_", res,"_PC_",pc, ".rds"),
    readRDS
  )
  
  scores <- dplyr::bind_rows(scores) %>%
    dplyr::group_by(res) %>%
    dplyr::mutate("n_clusters" = dplyr::n()) %>%
    dplyr::ungroup()
  return(scores)
}

Median <- function(x, interval = 0.95, R = 25000, type = "bca") {
  # Define median to take data and indices for use with boot::
  med <- function(data, indices) {
    resample <- data[indices]
    return(median(resample))
  }
  
  # Calculate intervals
  boot_data <- boot::boot(data = x, statistic = med, R = R)
  boot_ci <- boot::boot.ci(boot_data, conf = interval, type = type)
  
  # Extract desired statistics
  ci <- list(
    low_med = boot_ci$bca[4],
    med = boot_ci$t0,
    high_med = boot_ci$bca[5]
  )
  return(ci)
}

AddMedian <- function(x){
  scores <- x
  meds <- scores %>%
    dplyr::group_by(res) %>%
    dplyr::summarise(
      "boot" = list(Median(avg_sil)),
      "n_clusters" = mean(n_clusters)
    ) %>%
    tidyr::unnest_wider(boot)
  return(meds)
}

FindBestRes <- function(x){
  scores <- x
  BestRes<-scores %>%
    filter(avg_sil >= 0.7)%>%
    group_by(res) %>%
    count %>% ungroup %>%
    right_join(scores,by="res")%>%mutate("percent" = n/n_clusters)%>%
    mutate("med" = n/n_clusters)%>%
    filter(percent >= 0.5)%>%
    dplyr::arrange(percent) %>%
    head(n = 1) %>%
    dplyr::pull(res)
  print(paste0("the best res is"," ",as.character(BestRes)))
  return(as.character(BestRes))
}
