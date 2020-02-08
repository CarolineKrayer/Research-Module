################################################################################
######## This file generates the quantile-quantile plots in Figure D.4 #########
################################################################################

# Create quantile-quantile-plots for different degrees of undersmoothing.

rm(list=ls())

library(ggplot2)
library(rdd)
library(rddtools)

source("simulation_functions.R")
source("manipulation_test.R")

# Simulation parameters
n <- 20000
M <- 1000
select_frac <- 0.2
undersmoothing <- c(0.9, 0.75, 0.5, 0.25)

for (us in undersmoothing){
  set.seed(12345)
  z_scores <- sapply(rep(1, M), function(x) simul_fun(n=n, distribution="stdnorm",
                                                      select_frac=select_frac, 
                                                      h=us, b=NULL, out="z"))
  
  z_scores_df <- data.frame(z_scores)
  final = ggplot(z_scores_df, aes(sample = z_scores)) +
    stat_qq(shape=21, alpha=0.5, size=3, color="blue") +
    stat_qq_line(size=0.7) +
    theme_minimal() +
    theme(text = element_text(size = 14, hjust=0.5)) +
    xlab("Quantile of standard normal distribution") + 
    ylab("Quantile of test distribution")
  ggsave(final, file=paste0("../Figures/QQplot_", n, "obs_undersmoothing", us, "_", M, "rep.png"),
         height=6, width=6)
}

