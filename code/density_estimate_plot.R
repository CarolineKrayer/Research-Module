################################################################################
######### This file generates the density estimate plot in Figure D.1 ##########
################################################################################

# Implement the local linear density estimation
# and plot the density estimate with the
# respective first step histogram.


rm(list=ls())


library(ggplot2)

source("simulation_functions.R")
source("manipulation_test.R")

# Housekeeping
out_path <- "../Figures/"


set.seed(12345)

# Simulation parameters
n <- 5000
select_frac <- .2



loc_lin_density_estimate <- function(R, cutpoint, bin=NULL, bw=NULL,
                                     binning_proced="mccrary_2008"){

    # Perform discontinuity estimation by means of
    # local linear density estimation for poverty scores for each year.
    # Data from Camacho and Conover (2011).


    # Args:
    #       R (vector of floats): Data on running variable.
    #       cutpoint (float): RDD cutoff.
    #       bw (float): Bandwidth. Default is NULL and McCrary's procedure
    #                              is performed.
    #       bin (float): Binsize. Default is McCrary's procedure.
    #       binning_proced (str): Specify binning procedure, either
    #                             McCrary (2008) (default) or
    #                             Camacho and Conover (2011) ("camacho_conover_2011")
    #       out (str): Either return dataframe, binsize or bandwidth.

    # Out:
    #     Return list of dataframe (containing data to plot
    #     density estimates), binsize or bandwidth.


    # Perform manipulation test based on local linear density estimation.
    dens_test_object <- manipulation_test(runvar=R, cutpoint=cutpoint,
                                          bin=bin, bw=bw,
                                          binning_proced=binning_proced)

    # Retrieve relevant outputs to reconstruct density estimation
    # for plotting purposes.
    bw <- dens_test_object$bw
    bin <- dens_test_object$bin
    nbins <- dens_test_object$nbins

    cellmp <- dens_test_object$data$cellmp
    cellval <- dens_test_object$data$cellval

    # Estimate density to either side of the cutpoint using a triangular kernel.
    # Left of cutoff.
    d.l <- data.frame(cellmp=cellmp[cellmp <= cutpoint],
                      cellval=cellval[cellmp <= cutpoint],
                      dist=NA, est=NA, lwr=NA, upr=NA)

    # Add one additional observation such that plot includes
    # estimate and confidence at cutoff.
    if (binning_proced == "camacho_conover_2011"){
      values_at_c <- data.frame(cellmp=cutpoint+bin/2,
                                cellval=cellval[cutpoint+1],
                                dist=NA, est=NA, lwr=NA, upr=NA)
      d.l <- rbind(d.l, values_at_c)
    }

    for (i in 1:nrow(d.l)){
      # Fill dataframe with data on cell midpoints,
      # density function estimate at cell midpoints,
      # and lower and upper confidence.

      # Distance to cutoff.
      d.l$dist <- d.l$cellmp - d.l[i,"cellmp"]

      w <- kernelwts(d.l$dist, 0, bw, kernel="triangular")

      # Predict density function's value at zero distance to cutoff,
      # i.e. at cutoff.
      newd <- data.frame(dist=0)

      # Predict density function's value at zero distance to cutoff,
      # based on weighted least squares estimation.
      pred <- predict(lm(cellval ~ dist, weights=w, data=d.l),
                      interval="confidence", newdata=newd)

      # Density point estimate at bin midpoint.
      d.l$est[i] <- pred[1]

      # Lower confidence interval.
      d.l$lwr[i] <- pred[2]

      if (pred[2] < 0){
        # Deal with special case with confidence interval
        # entering negative values. Cannot be in density
        # estimation, and for purpose of plotting set zero.
        d.l$lwr[i] <- 0
      }
      # Upper confidence interval.
      d.l$upr[i] <- pred[3]
    }

    # Right of cutoff.

    d.r <- data.frame(cellmp=cellmp[cellmp > cutpoint],
                      cellval=cellval[cellmp > cutpoint],
                      dist=NA, est=NA, lwr=NA, upr=NA)

    # Add one additional observation such that plot includes
    # estimate and confidence at cutoff.
    if (binning_proced == "camacho_conover_2011"){
      values_at_c <- data.frame(cellmp=cutpoint+bin/2+.2,
                                cellval=cellval[cutpoint+1],
                                dist=NA, est=NA, lwr=NA, upr=NA)
      d.r <- rbind(values_at_c, d.r)
    }

    for (i in 1:nrow(d.r)){
      d.r$dist <- d.r$cellmp - d.r[i,"cellmp"]
      w <- kernelwts(d.r$dist, 0, bw, kernel="triangular")
      newd <- data.frame(dist=0)
      pred <- predict(lm(cellval ~ dist, weights=w, data=d.r),
                      interval="confidence", newdata=newd)
      d.r$est[i] <- pred[1]
      d.r$lwr[i] <- pred[2]
      d.r$upr[i] <- pred[3]
    }

    # Prepare data for plotting.
    df <- rbind(d.l, d.r)
    df$lest <- df$est
    df$rest <- df$est
    df$llwr <- df$lwr
    df$rlwr <- df$lwr
    df$lupr <- df$upr
    df$rupr <- df$upr

    for (i in (dim(d.l)[1]+1):dim(df)[1]){
      # Fill estimate and confidence interval of left density with NA
      # for values right to cutoff.
      df$lest[i] <- NA
      df$llwr[i] <- NA
      df$lupr[i] <- NA
    }

    for (i in 1:dim(d.l)[1]){
      df$rest[i] <- NA
      df$rlwr[i] <- NA
      df$rupr[i] <- NA
    }

    return(list(bw=bw,
                bin=bin,
                nbins=nbins,
                df=df))
}


loc_lin_density_plot <- function(df, xlim=NULL, ylim=NULL, cutpoint=0, bin,
                                 xtitle="", ytitle="", hist_only=FALSE,
                                 n, nbins, bw,
                                 title="Local Linear Density Estimation",
                                 camacho=FALSE,
                                 year=1900){
    # Plot running variable's histogram and local linear
    # density estimates on both sides of RDD cutoff.

    # Args:
    #       df (data.frame): Data containing loc. lin. dens. estimates
    #                         and 95% confidence interval.
    #                         Dataframe must be output of loc_lin_density_estimate.
    #       cutpoint (float): RDD cutoff that is drawn in plot.
    #       xlim (vector): Vector of length two containing interval on x-axis.
    #       ylim (vector): Vector of length two containing interval on y-axis.
    #       xtitle (str): Title of x-axis.
    #       ytitle (str): Title of y-axis.
    #       hist_only (boolean): If TRUE plot histogram only.
    #       n (int): Number of observations. Required for labeling plot only.
    #       nbins (int): Number of bins. Required for labeling plot only.
    #       bw (float): Bandwidth. Required for labeling plot only.
    #       camacho (boolean): Allows specific arrangements required
    #                          for plots to replicate Camacho and Conover (2011)
    #       year (int): Declares year of data from Camacho and Conover (2011),
    #                    required for saving plots appropriately.
    #       bin (float): Binsize of underlying histogram required for plot.

    # Out:
    #     Save respective plot in folder.


    # New design
    p <- ggplot(data=df, aes(x=df$cellmp)) +
          geom_col(aes(y=df$cellval), width=.85*bin,
                   fill="#5e5ef7", alpha=1) +
          labs(title=title) +
          labs(x=xtitle, y=ytitle) +
          theme_minimal() +
          theme(plot.margin=margin(0.4, 0.4, 0.4, 0.4, "cm"),
                text = element_text(size=18)) +
          xlim(xlim[1], xlim[2]) +
          ylim(ylim[1], ylim[2])

    vec_temp <- df$cellmp <= cutpoint
    binnum_c <- length(vec_temp[vec_temp==TRUE])

    if (hist_only == FALSE){
        p <- p + geom_line(aes(y=df$lest), size=1, colour="red", alpha=.8) +
                 geom_line(aes(y=df$rest), size=1, colour="red", alpha=.8) +
                 geom_ribbon(data=subset(df, df$cellmp[1:binnum_c] <= cutpoint),
                                         aes(ymin=df$llwr, ymax=df$lupr),
                                         fill="red", alpha="0.35") +
                 geom_ribbon(data=subset(df, df$cellmp[(binnum_c+1):dim(df)[1]] > cutpoint),
                                         aes(ymin=df$rlwr, ymax=df$upr),
                                         fill="red", alpha="0.35")
                 #labs(subtitle=paste0(n, " observations in ", nbins, " bins",
                 #                   "\nCutoff at ", cutpoint, " and share of manipulators is ",
                 #                     select_frac, "\nBandwidth: ", round(bw, 3)))
    }
    if (hist_only == TRUE){
        p <- p + labs(subtitle=paste0(n, " observations in ", nbins, " bins",
                                      "\nCutoff at ", cutpoint, " and share of manipulators is ",
                                      select_frac, "\n"))
    }

    # Save final plot.
    if (camacho == TRUE){
        final <- p + geom_vline(xintercept=(cutpoint+.5),
                                color="#262626", size=.9,
                                linetype="longdash")
        ggsave(final,
               file=paste0(out_path, "density_estimate_plot.png"),
               height=6, width=10)
    }

    if (camacho == FALSE){
        final <- p + geom_vline(xintercept=(cutpoint),
                                    color="#262626", size=.9,
                                    linetype="longdash")
        ggsave(final,
               file=paste0(out_path, "density_estimate_plot.png"),
               height=6, width=10)
    }
}



# Actually call functions and plot.
# All parameters of the Data Generating Process
# are specified at the script's beginning.


# Standard normal.
out <- draw_runvar(n, distribution="stdnorm")
R <- out$R
cutoff <- out$cutoff
impact_region <- out$impact_region

R_manipu <- sapply(R, function(r) select_fun(r=r, cutoff=cutoff,
                                             select_frac=select_frac,
                                             impact_region=impact_region))

out <- loc_lin_density_estimate(R_manipu, cutpoint=cutoff)
df <- out$df
bin <- out$bin
nbins <- out$nbins
bw <- out$bw

loc_lin_density_plot(df, xlim=c(-4, 4), ylim=c(0, .65),
                     cutpoint=cutoff, bin=bin,
                     n=n, nbins=nbins, bw=bw,
                     hist_only=FALSE,
                     title="",
                     xtitle=paste0("Running variable"))
