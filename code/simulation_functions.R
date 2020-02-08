
# Implement the data generating process and
# build the simulation framework. This module
# consists only of functions that are called in
# other modules.

library(rdd)
library(rdrobust)
library(boot)


draw_runvar <- function(n, distribution="stdnorm"){
    # Draw n observations of running variable.

    # Args:
    #     n (int): Number of observations.
    #     distribution (str): Choose distribution to draw
    #                          running variable from.
    #                          Cutoff and impact regions are
    #                          defined to ensure same amount of
    #                          individuals move over cutoff.
    #                          Default is N(0, 1).

    if (distribution == "stdnorm"){
        # Benchmark distribution.
        R <- rnorm(n, mean=0, sd=1)
        cutoff <- -.25
        impact_region <- .5
    }

    if (distribution == "uniform"){
        R <- runif(n, min=0, max=1)
        cutoff <- .5
        impact_region <- .17466
    }

    if (distribution == "uniform_large_var"){
        R <- runif(n, min=0, max=3.5)
        cutoff <- 2
        impact_region <- 0.61131
    }

    if (distribution == "f_distr"){
        R <- rf(n, 10, 20)
        cutoff <- 1.3
        impact_region <- .292 # For F(10, 20) at cutoff = 1.3
        #cutoff <- .6
        #impact_region <- .2923 # For F(10, 20) at cutoff = .6
    }
    return(list(R=R,
                cutoff=cutoff,
                impact_region=impact_region))
}


select_fun <- function(r, cutoff, select_frac, impact_region){
    # Define selection/manipulation procedure.

    # Individuals with values in impact region below cutoff
    # perfectly continuously move over cutoff with probability based
    # on fraction of individuals that can perfectly manipulate.
    if (cutoff-impact_region < r && r < cutoff){
        selection <- rbinom(1, size=1, prob=select_frac)
        r_manipu <- r + abs(impact_region) * selection
        return(r_manipu)
    }
    else{return(r)}
}


simul_fun <- function(n, distribution="stdnorm",
                      select_frac, h=NULL, b=NULL, out="theta"){
    # Simulate discontinuity gap at cutoff
    # through local linear density estimation.

    # Args:
    #       distribution (str): Specifies distribution running variable
    #                            is drawn from. Default is "stdnorm",
    #                            "uniform", and "f_distr" also possible.
    #       n (int): Number of observations.
    #       select_frac (float): Fraction of perfect manipulators.
    #       h (float): 1 or NULL takes default McCrary bandwidth. A value
    #                   smaller than one denotes the degree of undersmoothing,
    #                   i.e. percentage of default
    #                   McCrary bandwidth.
    #       b (float): Binsize. Default is McCrary's procedure.
    #       out (str): Select return value of function. Default is "theta".
    #                  Options are: "theta", "se", "z", "p", "data"

    # Returns:
    # Depending on Argument.
    #       theta (float): Discontinuity estimate transformed by ln().
    #       se (float): Standard error of theta.
    #       z (float): Z-statistic of the test.
    #       p (float): P-value of the test.
    #       data (dataframe): Dataframe for binning of the histogram.

    out_run_var <- draw_runvar(n, distribution=distribution)
    R <- out_run_var$R
    cutoff <- out_run_var$cutoff
    impact_region <- out_run_var$impact_region


    R_manipu <- sapply(R, function(r) select_fun(r=r, cutoff=cutoff,
                                                 select_frac=select_frac,
                                                 impact_region=impact_region))

    if (is.null(h) || h == 1){
        # Perform McCrary manipulation test.
        #dens_test_object <- DCdensity(runvar=R_manipu, cutpoint=cutoff,
        #                              bin=b, bw=NULL, ext.out=TRUE, plot=FALSE)
        dens_test_object <- manipulation_test(runvar=R_manipu, cutpoint=cutoff,
                                              bw=NULL, bin=NULL,
                                              binning_proced="mccrary_2008")

        # Collect relevant variables for calculating power.
        theta <- dens_test_object$theta
        se <- dens_test_object$se
        z <- dens_test_object$z
        p <- dens_test_object$p
        b <- dens_test_object$binsize
        h <- dens_test_object$bw
        df_binned_data <- dens_test_object$data

        if (out=="theta"){return(theta)}
        else if (out=="se"){return(se)}
        else if (out=="z"){return(z)}
        else if (out=="p"){return(p)}
        else if (out=="b"){return(b)}
        else if (out=="h"){return(h)}
        else if (out=="data"){return(df_binned_data)}
    }

    else if (h > 0 && h<1){
        # Mc Crary's standard bandwidth choice.
        bw_mcc <- manipulation_test(runvar=R_manipu, cutpoint=cutoff,
                                    bin=b, bw=NULL,
                                    binning_proced="mccrary_2008")$bw

        # Undersmooth the default McCrary bandwidth.
        bw_undersmoothed <- h * bw_mcc

        dens_test_object <- manipulation_test(runvar=R_manipu, cutpoint=cutoff,
                                              bw=bw_undersmoothed, bin=b,
                                              binning_proced="mccrary_2008")

        # Collect relevant variables for calculating power.
        theta <- dens_test_object$theta
        se <- dens_test_object$se
        z <- dens_test_object$z
        p <- dens_test_object$p
        b <- dens_test_object$binsize
        h <- dens_test_object$bw
        df_binned_data <- dens_test_object$data

        if (out=="theta"){return(theta)}
        else if (out=="se"){return(se)}
        else if (out=="z"){return(z)}
        else if (out=="p"){return(p)}
        else if (out=="b"){return(b)}
        else if (out=="h"){return(h)}
        else if (out=="data"){return(df_binned_data)}
    }

    else{stop("Stop: specified degree of undersmoothing is incorrect.")}

}

# Layer 1. Vectorize by fraction of perfect manipulators.
vsimul_fun_select_frac <- Vectorize(FUN=simul_fun, vectorize.args="select_frac")

# Layer 2. Vary by number of observations.
vsimul_fun_nobs <- Vectorize(FUN=simul_fun, vectorize.args="n")

# Layer 3. Vary by bandwidth selection procedures.
vsimul_fun_bandwidth <- Vectorize(FUN=simul_fun, vectorize.args="h")


count_reject <- function(alpha, entry){
    # Entry (float) evaluates to 1 if smaller than significance level alpha.
    if (entry < alpha){return(1)}
    else{return(0)}
}
vcount_reject <- Vectorize(FUN=count_reject, vectorize.args="entry")


powerfunction <- function(alpha, p_output){
    # Simulate the power of the density test in rejecting H0.

    # Args:
    #       alpha (float): Significance level used for evaluation of the
    #                       test decision.
    #       p_output (list): Multi-dimensional List of p-values computed
    #                         for different runs and shares of manipulators.

    # Returns:
    #       num_rejections (list): Averaged number of rejections of the test for
    #                               every share of manipulators specified.

    share_rejections <- c(rep(0, dim(p_output)[1]))
    for (i in 1:dim(p_output)[1]){
        num_rejections <- vcount_reject(entry=p_output[i,], alpha=alpha)
        share_rejections[i] <- mean(num_rejections)
    }
    return(share_rejections)
}
