######################################################################################
########### This file generates the power functions in Figure 4.1 and D.2 ############
######################################################################################


# Simulate the power of the density test 
# for different distributions and sample sizes
# when increasing the share of perfect 
# manipulators. 
# Uncomment the distribution of interest,
# default is standard normal.

rm(list=ls())

library(ggplot2)
library(rdrobust)
library(sjstats)
library(parallel)

source("simulation_functions.R")
source("manipulation_test.R")


# Housekeeping
out_path <- "../Figures/"

# Define function.
generate_powerfunction <- function(cl, n_vec, M, select_frac, h, b,
                                   distribution, alpha, out="p"){
  matrix <- matrix(select_frac)
  for(n in n_vec){
    set.seed(12345)
    # Compute p-values of the test.
    mout_p_val <- parSapply(cl, rep(1, M), function(x) vsimul_fun_select_frac(
                      n, distribution, select_frac, h, b, out="p"))
    # Compute power of the test in rejecting H0.
    power <- powerfunction(alpha, mout_p_val)
    matrix = cbind(matrix, power)
  }
  return(matrix)
}


# Call function
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterCall(cl, function() { source("simulation_functions.R") })
clusterCall(cl, function() { source("manipulation_test.R") })

matrix <- generate_powerfunction(cl=cl, 
                                 n_vec=c(500,1000,5000),
                                 distribution="stdnorm", 
#                                distribution="uniform",
#                                distribution="uniform_large_var",
#                                distribution="f_distr",
                                 M=5000, 
                                 select_frac=c(seq(0, .5, .005)), 
                                 alpha=.05, 
                                 h=NULL, 
                                 b=NULL)

stopCluster(cl)

# Create dataframe with results.
colnames(matrix) <- c("select_frac_vec", "n500", "n1000", "n5000")
df <- data.frame(matrix)

# Plot the powerfunction.
final = ggplot(df, aes(x=select_frac_vec)) + 
  geom_point(aes(y=n500, color="n = 500"), alpha=.7, size=1.4) +
  geom_point(aes(y=n1000, color="n = 1000"), alpha=.7, size=1.4) +
  geom_point(aes(y=n5000, color="n = 5000"), alpha=.7,  size=1.4) +
  labs(
         x = "Share of perfect manipulators",
         y = (expression(beta))) + 
  geom_hline(yintercept = 0.05, lty=2) +
  geom_text(aes(x=0.4, y=0.08, label='alpha==0.05'), 
            parse=TRUE, size=6) +
  theme_minimal() +
         scale_color_discrete(name = "Sample size") +
  theme(plot.margin=margin(0.4, 0.4, 0.4, 0.4, "cm"),
        text = element_text(size=18))

# Save plot.
ggsave(final, file=paste0(out_path, "powerfunction.png"),
       height=6, width=8)
