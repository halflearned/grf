library(causalbenchmark)
library(grf)
library(data.table)
library(tidyverse)
Macintosh HD⁩ ▸ ⁨Users⁩ ▸ ⁨vitorh⁩ ▸ ⁨Downloads⁩ ▸ ⁨centering_step_sims_results⁩
setwd("~/Downloads/centering_step_sims_results/results*.csv")
files <- list.files("~/Downloads/centering_step_sims_results/")

DT <- fread(files[2])

for (i in c(2:length(files))){
  DT <- rbind(DT, fread(files[i]), fill = TRUE)
}

DT <- DT[rule_num == 1, rule := "num_trees"]
DT <- DT[rule_num == 2, rule := "max(100,num_trees/10)"]
DT <- DT[rule_num == 3, rule := "max(100, sqrt(num_trees))"]
DT <- DT[rule_num == 4, rule := "max(100, num_trees/4)"]


plot_function <- function( N_val){
  DT[N == N_val,] %>% ggplot(aes(as.factor(num_trees),
                                 loss1,
                                 fill= as.factor(rule),
                                       as.factor(rule))) +
    stat_boxplot(geom ='errorbar', width = 0.6) +
    geom_boxplot(width = 0.6, lwd=.1, outlier.size = .1) +
    facet_wrap(~ df_name, scales = "free") +
    theme_bw() +
    theme(legend.position="bottom") +
    theme(text = element_text(size=8),
          axis.text.x = element_text(angle=0, hjust=1)) +
    theme(legend.title = element_blank()) +
    xlab("num_trees") +
    ylab("mse") +
    ggtitle(paste0("N = ", N_val)) +
    ggsave(paste0("~/Desktop/centering", N_val, ".pdf"), width = 7, height = 7)
}

plot_function(1000)
plot_function(2000)
plot_function(5000)
