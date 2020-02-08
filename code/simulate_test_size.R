################################################################################
######## This file generates the data for the test's size in Table D.1 #########
################################################################################


# Simulate the test's size for varying degrees
# undersmoothing the pilot bandwidth.

rm(list=ls())


library(pbapply)

source("simulation_functions.R")
source("manipulation_test.R")


#######################################

# Simulation parameters
M <- 1000
n <- 20000

# RDD specific parameters
select_frac <- 0

# Test specific parameters
alpha <- .05

undersmooth.vec <- c(seq(.25, 1, .05))

#######################################



# Draw running variable values from continuous density
# and run local linear density estimation with selection
# share of zero (i.e. no discontinuity aside via chance).
# Run for different degrees of undersmoothing to see
# if size of test is influenced.
set.seed(12345)
p_output <- pbsapply(rep(1, M), function(x) vsimul_fun_bandwidth(n=n,
                                                                 distribution="stdnorm",
                                                                 select_frac=select_frac,
                                                                 h=undersmooth.vec,
                                                                 b=NULL, out="p"))
share_rejections <- powerfunction(alpha, p_output)
df_test_size_stdnorm <- data.frame(t(share_rejections))
names(df_test_size_stdnorm) <- undersmooth.vec


###
set.seed(12345)
p_output <- pbsapply(rep(1, M), function(x) vsimul_fun_bandwidth(n=n,
                                                                 distribution="uniform",
                                                                 select_frac=select_frac,
                                                                 h=undersmooth.vec,
                                                                 b=NULL, out="p"))
share_rejections <- powerfunction(alpha, p_output)
df_test_size_uniform <- data.frame(t(share_rejections))
names(df_test_size_uniform) <- undersmooth.vec


###
set.seed(12345)
p_output <- pbsapply(rep(1, M), function(x) vsimul_fun_bandwidth(n=n,
                                                                 distribution="f_distr",
                                                                 select_frac=select_frac,
                                                                 h=undersmooth.vec,
                                                                 b=NULL, out="p"))
share_rejections <- powerfunction(alpha, p_output)
df_test_size_fdistr <- data.frame(t(share_rejections))
names(df_test_size_fdistr) <- undersmooth.vec


# Combine results in one dataframe.
df_test_size <- rbind(df_test_size_stdnorm,
                      df_test_size_uniform,
                      df_test_size_fdistr)


# Make table.
#table = xtable(output_bw, digits=3)
#print.xtable(table, type = "latex")
