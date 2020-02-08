################################################################################
########### This file generates the distribution plot in Figure D.5 ############
################################################################################


# Plot the test statistic's distributions for fixed
# number of observations, share of selectors,
# underlying standard normal distribution of running variable
# and varying degree of undersmoothing.

rm(list=ls())

library(ggplot2)
library(reshape2)

source("simulation_functions.R")
source("manipulation_test.R")


# Housekeeping
out_path <- "../Figures/"


# Simulation parameters
n <- 1000
M <- 5000
distr <- "stdnorm"

# RDD specific parameters
select_frac <- .2
undersmoothing <- .5


# Get data.
set.seed(12345)


# Compute test statistic for default bandwidth.
z_scores <- sapply(rep(1, M), function(x) simul_fun(n=n, distribution=distr,
                                                    select_frac=select_frac,
                                                    h=NULL, b=NULL, out="z"))

# Compute test statistic for 50% undersmoothing of default bandwidth.
z_scores_us <- sapply(rep(1, M), function(x) simul_fun(n=n, distribution=distr,
                                                       select_frac=select_frac,
                                                       h=undersmoothing, b=NULL, out="z"))

std_norm <- rnorm(n, mean=0, sd=1)

df_z_scores <- cbind(z_scores, z_scores_us, std_norm)
df_z_scores <- data.frame(df_z_scores)

# Stack data on top, required by ggplot.
df_m <- melt(df_z_scores)

# Plot densities.
p <- ggplot(df_m, aes(x=value, group=variable)) +
     geom_density(aes(fill=variable, color=variable), alpha=.65, size=1.05) +
     geom_hline(yintercept=0, color="#ebebeb", size=1.05) +
     xlim(-6.5, 6.5) +
     scale_fill_manual(values=c("#f8766c", "#00bec4", "#c67bff"), #7bad00
                       name="", labels=c("Pilot bandwidth",
                                         paste0(undersmoothing*100, "% undersmoothing"),
                                         "Limiting distribution")) +
     scale_color_manual(values=c("#f8766c", "#00bec4", "#c67bff")) +
     xlab("Test statistic") +
     ylab("Density") +
     theme_minimal() +
     theme(plot.margin=margin(0.4, 0.4, 0.4, 0.4, "cm"),
           text=element_text(size=18))

# Removes default legend.
final <- p + guides(color=FALSE)

# Save plot.
ggsave(final, file=paste0(out_path, "consistency_centering_", distr, "_", n, "obs_", M, "rep.png"),
        height=6, width=10)
