
# Implement McCrary (2008)'s manipulation test and
# use source code by Drew Dimmery <drewd@nyu.edu>.
# The purpose of this script is to follow implementation,
# understand code better, increase robustness in practical use,
# and comment sections. Further, it is adjusted to allow
# local linear density estimation based on a different binning technique
# required by a real data application (Camacho and Conover, 2011).

# This code heavily builds on Drew Dimmery's code in 'rdd' package 
# version 0.57, March 14, 2016. All mistakes are ours.
# Github: https://github.com/cran/rdd
# Package release: https://cran.r-project.org/web/packages/rdd/rdd.pdf



manipulation_test <- function(runvar, cutpoint, bw=NULL, bin=NULL,
                              binning_proced="mccrary_2008"){
    # Perform McCrary (2008) manipulation test
    # in regression discontinuity designs.

    # Args:
    #       runvar (vector of floats): Data on running variable.
    #       cutpoint (float): RDD cutoff.
    #       bw (float): Bandwidth. Default is NULL and McCrary's procedure
    #                              is performed.
    #       bin (float): Binsize. Default is McCrary's procedure.
    #       binning_proced (str): Specify binning procedure, either
    #                             McCrary (2008) (default) or specific for
    #                             Camacho and Conover (2011) ("camacho_conover_2011")

    # Returns:
    #       theta (float): Discontinuity estimate transformed by ln().
    #       se (float): Standard error of theta.
    #       z (float): Z-statistic of the test.
    #       p (float): P-value of the test.
    #       binsize (float): Binsize.
    #       bw (float): Bandwidth.
    #       cutpoint (float): RDD cutoff.
    #       data (dataframe): Dataframe for binning of the histogram.


    # Drop missing values.
    runvar <- runvar[complete.cases(runvar)]

    # Summary statistics frequently used.
    rn <- length(runvar)
    rsd <- sd(runvar)
    rmin <- min(runvar)
    rmax <- max(runvar)

    if (cutpoint<=rmin | cutpoint>=rmax){
        stop("Cutpoint must lie within range of runvar")
    }

    if (is.null(bin)){
        # Default binsize stems from McCrary (2008)'s
        # automatic bandwidth selection, step 1.
        bin <- 2*rsd*rn^(-1/2)

        # Manual Binning technique for Camacho and Conover 2011.
        if (binning_proced=="camacho_conover_2011"){
            bin <- 1
        }
    }


    if (binning_proced=="mccrary_2008"){
        # Original McCrary Binning technique.
        l <- floor((rmin - cutpoint)/bin)*bin + bin/2 + cutpoint # Midpoint of lowest bin
        r <- floor((rmax - cutpoint)/bin)*bin + bin/2 + cutpoint # Midpoint of highest bin
        lc <- cutpoint - (bin/2) # Midpoint of bin just left of cutoff
        rc <- cutpoint + (bin/2) # Midpoint of bin just right of cutoff
        j <- floor((rmax - rmin)/bin) + 2 # Total number of bins.
    }

    if (binning_proced=="camacho_conover_2011"){
        # Manual Binning technique for Camacho and Conover 2011.
        # Their data are already discretisized and binning has to
        # acknowledge the underlying data generating process.
        # Poverty scores range from 0 to 100.
        l <- 0
        r <- 98
        lc <- 47
        rc <- 48
        j <- 99
    }


    # Assign each value of running variable its bin (returns vector of bin numbers).
    binnum <- round((((floor((runvar - cutpoint)/bin)*bin + bin/2 + cutpoint) - l)/bin) + 1)

    # With Camacho and Conover (2011) data no need to assign bins as categorical data
    # already correspond to bin.
    if (binning_proced=="camacho_conover_2011"){binnum <- runvar}

    # Create histogram.
    cellval <- rep(0, j)
    for (i in seq(1, rn)){
        # Count observations of running variable in each bin.
        cnum <- binnum[i]

        # Increment index cnum by 1 since R indexing starts
        # with 1. E.g. individual has score zero, hence cnum
        # will be zero but cellval for zero cannot be accessed
        # by index zero. In fact, it corresponds to index 1.
        if (binning_proced=="camacho_conover_2011"){cnum <- cnum+1}

        cellval[cnum] <- cellval[cnum]+1
    }
    # Normalize bin height.
    cellval <- (cellval/rn) / bin

    # Calculate bin midpoints.
    cellmp <- seq(from=1, to=j, by=1)
    cellmp <- floor(((l + (cellmp - 1)*bin ) - cutpoint)/bin)*bin + bin/2 + cutpoint

    # With Camacho and Conover (2011) data cell midpoints coincide with categorical data.
    if (binning_proced=="camacho_conover_2011"){cellmp <- seq(from=0, to=r, by=1)}


    if (is.null(bw)){
        # Default bandwidth calculation follows McCrary (2008)'s procedure.

        # Bin number just left of cutoff.
        leftofc <-  round((((floor((lc - cutpoint)/bin)*bin + bin/2 + cutpoint) - l)/bin) + 1)
        # Bin number just right of cutoff.
        rightofc <- round((((floor((rc - cutpoint)/bin)*bin + bin/2 + cutpoint) - l)/bin) + 1)


        # With already discretisized data above algorithm is not stable.
        if (binning_proced=="camacho_conover_2011"){
            # Increment by one due to indexing, i.e. first entry is zero.
            leftofc <- 47 + 1
            rightofc <- 48 + 1
        }

        if (rightofc - leftofc != 1) {
            stop("Error occurred in bandwidth calculation")
        }

        # Consider data on left and right of cutoff seperately.
        cellvalleft <- cellval[1:leftofc]
        cellvalright <- cellval[rightofc:j]
        cellmpleft <- cellmp[1:leftofc]
        cellmpright <- cellmp[rightofc:j]


        # Estimate 4th order polynomial to the left.
        P.lm <- lm(cellvalleft ~ poly(cellmpleft, degree=4, raw=T))

        mse4 <- summary(P.lm)$sigma^2
        lcoef <- coef(P.lm)
        fppleft <- 2*lcoef[3] + 6*lcoef[4]*cellmpleft +
                   12*lcoef[5]*cellmpleft*cellmpleft
        hleft <- 3.348*(mse4*(cutpoint - l) / sum(fppleft*fppleft))^(1/5)


        # Estimate 4th order polynomial to the right.
        P.lm <- lm(cellvalright ~ poly(cellmpright, degree=4, raw=T))

        mse4 <- summary(P.lm)$sigma^2
        rcoef <- coef(P.lm)
        fppright <- 2*rcoef[3] + 6*rcoef[4]*cellmpright +
                    12*rcoef[5]*cellmpright*cellmpright
        hright <- 3.348*(mse4*( r - cutpoint ) / sum(fppright*fppright))^(1/5)

        bw = .5*(hleft + hright)
    }

    if( sum(runvar>cutpoint-bw & runvar<cutpoint)==0 |
    sum(runvar<cutpoint+bw & runvar>=cutpoint)==0 )
    stop("Insufficient data within the bandwidth.")


    # Estimate size of discontinuity at cutoff and compute manipulation test.
    cmp <- cellmp
    cval <- cellval
    padzeros <- ceiling(bw/bin)
    jp <- j + 2 * padzeros

    if(padzeros >= 1) {
        # naive way:
        cval <- cellval
        cmp <- cellmp
    }

    # Distance from cell midpoint to cutoff depicts regressor.
    dist <- cmp - cutpoint


    ### Implementation follows McCrary (2008).
    # Estimate density left and right of cutoff separately.

    # Density estimate to the left.
    # Compute weights for weighted-least-squares.
    w <- 1 - abs(dist/bw)
    w <- ifelse(w > 0, w*(cmp <= cutpoint), 0)
    w <- (w/sum(w)) * jp

    fhatl <- predict(lm(cval ~ dist , weights=w), newdata=data.frame(dist=0))[[1]]

    # Density estimate to the right.
    # Compute weights for weighted-least-squares.
    w <- 1 - abs(dist/bw)
    w <- ifelse(w > 0, w*(cmp > cutpoint), 0)
    w <- (w/sum(w)) * jp

    fhatr <- predict(lm(cval ~ dist , weights=w), newdata=data.frame(dist=0))[[1]]

    if (fhatl <= 0){
        fhatl <- NA
        warning("Density function estimate left to cutoff is zero or negative.
                Log is not defined, NA value for estimate produced.")
    }
    if (fhatr <= 0){
        fhatr <- NA
        warning("Density function estimate right to cutoff is zero or negative.
                Log is not defined, NA value for estimate produced.")
    }

    # McCrary estimate.
    thetahat <- log(fhatr) - log(fhatl)
    sethetahat <- sqrt( (1/(rn*bw)) * (24/5) * ((1/fhatr) + (1/fhatl)) )
    # Wald test stat under H0: theta=0.
    z <- thetahat/sethetahat
    p <- 2*pnorm(abs(z), lower.tail=FALSE)

    # Estimate of actual gap. (Not McCrary estimate)
    gap_mcc <- fhatl - fhatr


    ### Alternative implementation follows Camacho and Conover (2011)
    # Estimate density for whole support and use dummy variable to indicate
    # cutoff value. Its coefficient is jump in density.

    w <- 1 - abs(dist/bw)
    w <- ifelse(w > 0, w, 0)
    w <- (w/sum(w)) * jp

    jump <- cmp
    jump <- ifelse(cellmp <= cutpoint, 1, 0)
    jump_dist <- jump * dist

    llmDens <- lm(cval ~ dist + jump + jump_dist , weights=w)
    gap_cc <- llmDens$coefficients["jump"]
    # Robust standard error (reproduces Stata's robust)
    gap_cc_se <- coeftest(llmDens, vcov=vcovHC(llmDens, "HC1"))[6]


    return(list(theta=thetahat,
                se=sethetahat,
                z=z,
                p=p,
                binsize=bin,
                bw=bw,
                data=data.frame(cellmp, cellval),
                gap_cc_se=gap_cc_se,
                gap_cc=gap_cc,
                gap_mcc=gap_mcc,
                nbins=j
              ))
}
