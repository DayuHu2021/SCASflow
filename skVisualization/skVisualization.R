library(Seurat)
library(ggplot2)
library(tidyverse)
library(tibble)
library(dplyr)
library(ggsci)
`%>%` <- magrittr::`%>%`
source("Functions.R")


results_path <- ("~/Desktop/resultspmbc/")
pc <- 30
res <- c(0.1,0.5,1, 1.5, 2, 4,8,16)



scores<-  DrawFile(results_path,res,pc)
scores_meds <- AddMedian(scores)
scores_meds$res <- as.numeric(scores_meds$res)
scores_meds <- arrange(scores_meds,scores_meds$res)

#Set threshold and the best resolution 
threshold <- 0.7
best_res <- FindBestRes(scores)





####Select the best res####
ggplot(scores_meds, aes(factor(res), med)) +
  geom_crossbar(
    aes(ymin = low_med, ymax = high_med),
    fill = alpha("#FFA500",0.5),
    size = 0.25
  ) +
  geom_hline(aes(yintercept = threshold), colour = "blue") +
  geom_vline(aes(xintercept = best_res), colour = "red") +
  geom_jitter(
    data = scores,
    aes(factor(res), avg_sil),
    size = 0.35,
    width = 0.15
  ) +
  scale_x_discrete("Resolution") +
  scale_y_continuous(
    "Silhouette Score",
    expand = c(0, 0),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_hgrid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  )

ggsave(
  filename = paste0(results_path,"silhouette_of_each_res_pc_",pc,".png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)


####Stable ratio####
stable_cluster<-scores %>%
  filter(avg_sil >= 0.7)%>%
  group_by(res) %>%
  count %>% ungroup %>%
  right_join(scores,by="res")%>%mutate("percent" = n/n_clusters)%>%
  mutate("med" = n/n_clusters)%>%
  dplyr::group_by(res)%>%
  dplyr::summarise(
    "percent" = mean(percent)
  )
stable_cluster$res <- as.numeric(stable_cluster$res)

ggplot(stable_cluster, aes(factor(res), percent)) +
  geom_crossbar(
    aes(ymin = 0, ymax = percent),
    fill = alpha("#FFA500",0.5),
    size = 0.25
  ) +geom_text(aes(label=paste0(100*round(a$percent, 2),"%")),position=position_dodge(width=0.9),
               vjust=0,colour="red",size=4)+
  geom_hline(aes(yintercept = 0.5), colour = "blue")+
  scale_x_discrete("Resolution") +
  scale_y_continuous(
    "The proportion of stable clusters to the population",
    expand = c(0, 0),
    limits = c(0, 1.1),
    breaks = seq(0, 1, 0.1),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_hgrid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  )+scale_fill_lancet()

ggsave(
  filename = paste0(results_path, "stable_ratio.png"),
  dpi = 600,
  height = 3.5,
  width = 3.5,
  units = "in"
)


####Select the best PCs####
pcs_num <- c(10,20,30,40,50)
res <- c(0.1,0.5,1,1.5,2,4,8,16)

for (n in pcs_num){
  scores<-  DrawFile(results_path,res,n)
  scores<- mutate(scores,pcs = n)
  scores_meds <- AddMedian(scores)
  scores_meds$res <- as.numeric(scores_meds$res)
  scores_meds <- arrange(scores_meds,scores_meds$res)
  assign(paste0("score",n) ,scores)
  assign(paste0("score_meds",n) ,scores_meds)
}
  
all<-rbind(score10,score20,score30,score40,score50)
all$res <- as.numeric(all$res)


b1 <- subset(all,res==16)
b1_meds <- b1 %>%
  dplyr::group_by(pcs) %>%
  dplyr::summarise(
    "boot" = list(boot_median(avg_sil)),
    "n_clusters" = mean(n_clusters),
    "pcs"=mean(pcs)
  ) %>%
  tidyr::unnest_wider(boot)



ggplot(b1_meds, aes(factor(pcs), med)) +
  geom_crossbar(
    aes(ymin = low_med, ymax = high_med),
    fill = alpha("#FFA500",0.5),
    size = 0.25
  )  +
  geom_jitter(
    data = b1,
    aes(factor(pcs), avg_sil),
    size = 0.35,
    width = 0.15
  ) +
  scale_x_discrete("Number of PCs") +
  scale_y_continuous(
    "Silhouette Score",
    expand = c(0, 0),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_hgrid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  )

ggsave(
  filename = paste0(results_path,"silhouette_of_each_res_16.png"),
  dpi = 600,
  height = 3.5,
  width = 3.5,
  units = "in"
)



