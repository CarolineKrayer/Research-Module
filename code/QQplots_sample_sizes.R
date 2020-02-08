################################################################################
# This file generates the quantile-quantile plots in Figure 4.2 and Figure D.3 #
################################################################################

# Create quantile-quantile-plots for different sample sizes.

rm(list=ls())

library(ggplot2)
library(rdd)
library(rddtools)

source("simulation_functions.R")
source("manipulation_test.R")

# Simulation parameters
n_vec <- c(500, 1000, 5000, 10000, 20000)
M <- 1000
select_frac <- 0.2

for(n in n_vec){
  set.seed(12345)
  z_scores <- sapply(rep(1, M), function(x) simul_fun(n=n, distribution="stdnorm",
                                                      select_frac=select_frac, 
                                                      h=NULL, b=NULL, out="z"))
  
  z_scores_df <- data.frame(z_scores)
  final = ggplot(z_scores_df, aes(sample = z_scores)) +
    stat_qq(shape=21, alpha=0.5, size=3, color="blue") +
    stat_qq_line(size=0.7) +
    theme_minimal() +
    theme(text = element_text(size = 14, hjust=0.5)) +
    xlab("Quantile of standard normal distribution") + 
    ylab("Quantile of test distribution")
    ggsave(final, file=paste0("../Figures/QQplot_", n, "obs_", M, "rep.png"),
            height=6, width=6)
}




