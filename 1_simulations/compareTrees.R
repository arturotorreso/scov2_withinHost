library(treespace)
library(tidyverse)
library(plotly)
library(phangorn)
library(RColorBrewer)
library(here)

setwd(here()) # should set working directory to project folder - if not, do manually

# 3 lists of 100 trees, where indices correspond. So real_tree[[1]] is a simulation, and cons_tree[[1]], hets_tree[[1]] are the corresponding inferred trees from the associated sequences
load("data/sim_trees.RData")

l <- length(real_tree) # 100

all_trees <- c(cons_tree,hets_tree,real_tree)
class(all_trees) <- "multiPhylo"

# compare the trees as unrooted. Metrics: KF, RF, wRF, SP, wSP

distances_no_scaling <- tibble(.rows=l)
distances_no_scaling$tree_no <- 1:l

distances_no_scaling <- distances_no_scaling %>%
  mutate("KF consensus" = sapply(1:l, function(x) KF.dist(real_tree[[x]], cons_tree[[x]], rooted=FALSE)) ) %>%
  mutate("RF consensus" = sapply(1:l, function(x) RF.dist(real_tree[[x]], cons_tree[[x]], rooted=FALSE))) %>%
  mutate("wRF consensus" = sapply(1:l, function(x) wRF.dist(real_tree[[x]], cons_tree[[x]], rooted=FALSE))) %>%
  mutate("SP consensus" = sapply(1:l, function(x) path.dist(real_tree[[x]], cons_tree[[x]], use.weight=FALSE)) ) %>%
  mutate("wSP consensus" = sapply(1:l, function(x) path.dist(real_tree[[x]], cons_tree[[x]], use.weight=TRUE)))  %>%
  mutate("KF hets" = sapply(1:l, function(x) KF.dist(real_tree[[x]], hets_tree[[x]], rooted=FALSE))) %>%
  mutate("RF hets" = sapply(1:l, function(x) RF.dist(real_tree[[x]], hets_tree[[x]], rooted=FALSE))) %>%
  mutate("wRF hets" = sapply(1:l, function(x) wRF.dist(real_tree[[x]], hets_tree[[x]], rooted=FALSE)) ) %>%
  mutate("SP hets" = sapply(1:l, function(x) path.dist(real_tree[[x]], hets_tree[[x]], use.weight=FALSE)) ) %>%
  mutate("wSP hets" = sapply(1:l, function(x) path.dist(real_tree[[x]], hets_tree[[x]], use.weight=TRUE)) )

# go wide to long
distances_no_scaling_long <- distances_no_scaling %>%
  gather(key="metric and tree type", value="value",
         colnames(distances_no_scaling)[2:11])
distances_no_scaling_long <- distances_no_scaling_long %>%
  mutate("metric" = c(rep("KF",l),rep("RF",l),rep("wRF",l),rep("SP",l),rep("wSP",l),
                      rep("KF",l),rep("RF",l),rep("wRF",l),rep("SP",l),rep("wSP",l))) %>%
  mutate("tree type" = c(rep("Consensus",l*5),rep("Hets",l*5)))



# some convenient functions for tidy plotting:
ylimit <- c(0.00425, 180, 0.05, 2100, 0.32)
ytitles <- c("Distance from real tree", rep("",4))
metric_short_names <- unique(distances_no_scaling_long$metric)
metric_long_names <- c("Kuhner Felsenstein", 
                       "Robinson Foulds",
                       "Weighted Robinson Foulds",
                       "Steel and Penny",
                       "Weighted Steel and Penny")

# plotting aesthetics
f1 <- list(
  family = "Arial, sans-serif",
  size = 28,
  color = "black"
)
f2 <- list(
  family = "Arial, sans-serif",
  size = 22,
  color = "black"
)
f3 <- list(
  family = "Arial, sans-serif",
  size = 16,
  color = "black"
)

plot_list_no_jitter <- lapply(1:5, function(m) {
  plot_ly(distances_no_scaling_long %>% filter(metric == metric_short_names[[m]]), 
          type="box", 
          # boxpoints = "all", jitter=1,
          # pointpos = 0, marker = list(size = 2),
          y=~value, color=~`metric and tree type`,
          colors=brewer.pal(10, "Paired")[c((2*m-1),2*m)]) %>%
    layout(
      xaxis=list(
        tickvals=list(paste0(metric_short_names[[m]]," consensus"), paste0(metric_short_names[[m]]," hets")),
        ticktext=list("Consensus","Within-host"),
        title=metric_long_names[[m]],
        titlefont = f2,
        tickfont = f3
      ),
      yaxis=list(
        range=c(0,ylimit[[m]]),
        tickfont = f3,
        title=ytitles[[m]],
        titlefont=f1
      ),
      showlegend = FALSE
    )
})

final_no_jitter <- subplot(plot_list_no_jitter, titleX = TRUE, titleY=TRUE)

orca(final_no_jitter, "plots/within-host_figure_4.pdf",
     width="1600")
orca(final_no_jitter, "plots/within-host_figure_4.png",
     width="1600")










