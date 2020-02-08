################################################################################
############## This file generates the results in Table 5.1 ####################
################################################################################

# Apply McCrary test to data on Poverty Index Score
# from first SISBEN dataset for two choices of bandwidth.

rm(list=ls())

library(dplyr)
library(xtable)
library(lmtest)
library(sandwich)
library(foreign)


# Housekeeping
data_path <- "../data/"


manipulation_test <-
    function(runvar,
             cutpoint,
             bw = NULL,
             bin = NULL,
             h = NULL,
             binning_proced = "mccrary_2008") {
        # Perform McCrary (2008) manipulation test
        # in regression discontinuity designs, accomodating
        # for different binning procedures and different
        # degrees of undersmoothing.

        # Args:
        #       runvar (vector of floats): Data on running variable.
        #       cutpoint (float): RDD cutoff.
        #       bw (float): Bandwidth. Default is NULL and McCrary's procedure
        #                              is performed.
        #       bin (float): Binsize. Default is McCrary's procedure.
        #       binning_proced (str): Specify binning procedure, either
        #                             McCrary (2008) (default) or specific for
        #                             Camacho and Conover (2011) ("camacho_conover_2011")
        #       h (float): 1 or NULL takes default McCrary bandwidth. A value
        #                  smaller than one denotes the degree of undersmoothing,
        #                  i.e. percentage of default McCrary bandwidth.)
        #       out (str): Select return value of function. Default is "theta".
        #                  Options are: "theta", "se", "z", "p", "data"

        # Returns:
        #       theta (float): Discontinuity estimate transformed by ln().
        #       se (float): Standard error of theta.
        #       z (float): Z-statistic of the test.
        #       p (float): P-value of the test.
        #       binsize (float): Binsize.
        #       bw (float): Bandwidth.
        #       cutpoint (float): RDD cutoff.
        #       data (dataframe): Dataframe for binning of the histogram.
        #       gap_cc (float): Gap size estimated via the procedure employed
        #                       by Camacho and Conover (2011)
        #       gap_cc_se (float): Standard error of the gap size estimated
        #                          via the procedure employed by Camacho and
        #                          Conover (2011)
        #       gap_mcc (float): Estimated gap via the McCrary procedure before
        #                        taking the log

        # Drop missing values.
        runvar <- runvar[complete.cases(runvar)]

        # Summary statistics frequently used.
        rn <- length(runvar)
        rsd <- sd(runvar)
        rmin <- min(runvar)
        rmax <- max(runvar)

        if (cutpoint <= rmin | cutpoint >= rmax) {
            stop("Cutpoint must lie within range of runvar")
        }

        if (is.null(bin)) {
            # Default binsize stems from McCrary (2008)'s
            # automatic bandwidth selection, step 1.
            bin <- 2 * rsd * rn ^ (-1 / 2)

            # Manual Binning technique for Camacho and Conover 2011.
            if (binning_proced == "camacho_conover_2011") {
                bin <- 1
            }
        }


        if (binning_proced == "mccrary_2008") {
            # Original McCrary Binning technique.
            l <-
                floor((rmin - cutpoint) / bin) * bin + bin / 2 + cutpoint # Midpoint of lowest bin
            r <-
                floor((rmax - cutpoint) / bin) * bin + bin / 2 + cutpoint # Midpoint of highest bin
            lc <-
                cutpoint - (bin / 2) # Midpoint of bin just left of cutoff
            rc <-
                cutpoint + (bin / 2) # Midpoint of bin just right of cutoff
            j <- floor((rmax - rmin) / bin) + 2 # Total number of bins.
        }

        if (binning_proced == "camacho_conover_2011") {
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
        binnum <-
            round((((
                floor((runvar - cutpoint) / bin) * bin + bin / 2 + cutpoint
            ) - l) / bin) + 1)

        # With Camacho and Conover (2011) data no need to assign bins as categorical data
        # already correspond to bin.
        if (binning_proced == "camacho_conover_2011") {
            binnum <- runvar
        }

        # Create histogram.
        cellval <- rep(0, j)
        for (i in seq(1, rn)) {
            # Count observations of running variable in each bin.
            cnum <- binnum[i]

            # Increment index cnum by 1 since R indexing starts
            # with 1. E.g. individual has score zero, hence cnum
            # will be zero but cellval for zero cannot be accessed
            # by index zero. In fact, it corresponds to index 1.
            if (binning_proced == "camacho_conover_2011") {
                cnum <- cnum + 1
            }

            cellval[cnum] <- cellval[cnum] + 1
        }
        # Normalize bin height.
        cellval <- (cellval / rn) / bin

        # Calculate bin midpoints.
        cellmp <- seq(from = 1, to = j, by = 1)
        cellmp <-
            floor(((l + (cellmp - 1) * bin) - cutpoint) / bin) * bin + bin / 2 + cutpoint

        # With Camacho and Conover (2011) data cell midpoints coincide with categorical data.
        if (binning_proced == "camacho_conover_2011") {
            cellmp <- seq(from = 0,
                          to = r,
                          by = 1)
        }


        if (is.null(bw)) {
            # Default bandwidth calculation follows McCrary (2008)'s procedure.

            # Bin number just left of cutoff.
            leftofc <-
                round((((
                    floor((lc - cutpoint) / bin) * bin + bin / 2 + cutpoint
                ) - l) / bin) + 1)
            # Bin number just right of cutoff.
            rightofc <-
                round((((
                    floor((rc - cutpoint) / bin) * bin + bin / 2 + cutpoint
                ) - l) / bin) + 1)


            # With already discretisized data above algorithm is not stable.
            if (binning_proced == "camacho_conover_2011") {
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
            P.lm <- lm(cellvalleft ~ poly(cellmpleft, degree = 4, raw = T))

            mse4 <- summary(P.lm)$sigma ^ 2
            lcoef <- coef(P.lm)
            fppleft <- 2 * lcoef[3] + 6 * lcoef[4] * cellmpleft +
                12 * lcoef[5] * cellmpleft * cellmpleft
            hleft <-
                3.348 * (mse4 * (cutpoint - l) / sum(fppleft * fppleft)) ^ (1 / 5)


            # Estimate 4th order polynomial to the right.
            P.lm <-
                lm(cellvalright ~ poly(cellmpright, degree = 4, raw = T))

            mse4 <- summary(P.lm)$sigma ^ 2
            rcoef <- coef(P.lm)
            fppright <- 2 * rcoef[3] + 6 * rcoef[4] * cellmpright +
                12 * rcoef[5] * cellmpright * cellmpright
            hright <-
                3.348 * (mse4 * (r - cutpoint) / sum(fppright * fppright)) ^ (1 / 5)

            bw = .5 * (hleft + hright)
        }

        # Specify degree of undersmoothing
        if (is.null(h) || h == 1) {
            bw = bw
        }

        else if (h < 1 && h > 0) {
            # Undersmooth the default McCrary bandwidth.
            bw <- h * bw
        }
        else{
            stop("Stop: specified degree of undersmoothing is incorrect.")
        }


        if (sum(runvar > cutpoint - bw & runvar < cutpoint) == 0 |
            sum(runvar < cutpoint + bw & runvar >= cutpoint) == 0)
            stop("Insufficient data within the bandwidth.")


        # Estimate size of discontinuity at cutoff and compute manipulation test.
        cmp <- cellmp
        cval <- cellval
        padzeros <- ceiling(bw / bin)
        jp <- j + 2 * padzeros

        if (padzeros >= 1) {
            # Manually add observations (here: histogram value pairs (Xi, Yi))
            # to left and right of support of density with zero probability mass.

            # Original implementation:
            #cval <- c(rep(0, padzeros), cellval, rep(0, padzeros))
            #cmp <- c(seq(l-padzeros*bin, l-bin, bin), cellmp,
            #         seq(r+bin, r+padzeros*bin, bin))

            # naive way:
            cval <- cellval
            cmp <- cellmp
        }

        # Distance from cell midpoint to cutoff depicts regressor.
        dist <- cmp - cutpoint


        ### 1. way of computing discontinuity estimate (McCrary, 2008)
        #   Estimate density left and right of cutoff, and take difference.

        # Density estimate to the left.
        # Compute weights for weighted-least-squares.
        w <- 1 - abs(dist / bw)
        w <- ifelse(w > 0, w * (cmp <= cutpoint), 0)
        w <- (w / sum(w)) * jp
        fhatl <-
            predict(lm(cval ~ dist , weights = w), newdata = data.frame(dist = 0))[[1]]

        # Density estimate to the right.
        # Compute weights for weighted-least-squares.
        w <- 1 - abs(dist / bw)
        w <- ifelse(w > 0, w * (cmp > cutpoint), 0)
        w <- (w / sum(w)) * jp
        fhatr <-
            predict(lm(cval ~ dist , weights = w), newdata = data.frame(dist = 0))[[1]]

        if (fhatl <= 0) {
            fhatl <- NA
            warning(
                "Density function estimate left to cutoff is zero or negative.
                Log is not defined, NA value for estimate produced."
            )
        }
        if (fhatr <= 0) {
            fhatr <- NA
            warning(
                "Density function estimate right to cutoff is zero or negative.
                Log is not defined, NA value for estimate produced."
            )
        }

        # Estimate of actual gap. (Not McCrary estimate)
        gap_mcc <- fhatl - fhatr


        ### 2. way of computing discontinuity estimate (Camacho and Conover, 2011)
        #   Estimate density for whole support and use dummy variable to indicate
        #   cutoff value. Its coefficient is jump in density.

        w <- 1 - abs(dist / bw)
        w <- ifelse(w > 0, w, 0)
        w <- (w / sum(w)) * jp

        jump <- cmp
        jump <- ifelse(cellmp <= cutpoint, 1, 0)
        jump_dist <- jump * dist
        cval <- cval * 100

        llmDens <- lm(cval ~ dist + jump + jump_dist , weights = w)
        gap_cc <- llmDens$coefficients["jump"]
        # Robust standard error (reproduces Stata's robust)
        gap_cc_se <- coeftest(llmDens, vcov = vcovHC(llmDens, "HC1"))[6]


        # McCrary estimate.
        thetahat <- log(fhatr) - log(fhatl)
        sethetahat <-
            sqrt((1 / (rn * bw)) * (24 / 5) * ((1 / fhatr) + (1 / fhatl)))
        # Wald test stat under H0: theta=0.
        z <- thetahat / sethetahat
        p <- 2 * pnorm(abs(z), lower.tail = FALSE)

        return(
            list(
                theta = thetahat,
                se = sethetahat,
                z = z,
                p = p,
                binsize = bin,
                bw = bw,
                data = data.frame(cellmp, cellval),
                gap_cc_se = gap_cc_se,
                gap_cc = gap_cc,
                gap_mcc = gap_mcc
            )
        )
    }


# Get data.
#
# Create subset of first SISBEN dataset
# (only observations of Poverty Index Score
# by year).
sisben <- read.dta(paste0(data_path, "sisben_aejep.dta"))
sisben_small <- with(sisben, data.frame(score=puntaje, date=fencuesta))
year <- format(as.Date(sisben_small$date, format="%d/%m/%Y"),"%Y")
sisben_small <- mutate(sisben_small, year)
sisben_small <- select(sisben_small, -matches("date"))
sisben_small <- filter(sisben_small, year >= 1994 & year <= 2003)
rm(sisben, year)

# Alternatively, load pre-processed data, "sisben_small.csv".
#sisben_small <- read.table(paste0(data_path, "sisben_small.csv"), header = TRUE)
#sisben_small <- data[, -3]
#names(sisben_small)[1] <- "score"
#names(sisben_small)[2] <- "year"

# Produce table to compare Camacho and Conover (2011) results
# with McCrary density test estimated gap and decisions.
# No undersmoothing.
output = data.frame()
output <- do.call(rbind.data.frame,
                  lapply(1994:2003,
                         function(i) {
                             sisben_year = filter(sisben_small, sisben_small$year == i)
                             out = manipulation_test(
                                 sisben_year$score,
                                 cutpoint = 47,
                                 bw = NULL,
                                 bin = NULL,
                                 h = NULL,
                                 binning_proced = "camacho_conover_2011"
                             )
                             result = cbind(out$theta, out$se, out$p)
                             return(result)
                         }))
rownames(output) = 1994:2003
colnames(output) = c("Estimator_McCrary", "SE_McCrary", "p-value")

# Uncomment the following lines to generate latex table
#table = xtable(output, digits=3)
#print.xtable(table, type = "latex")

# Produce table to show McCrary's density test decisions
# for 50 percent undersmoothing of the pilot bandwidth.
output_bw = data.frame()
output_bw <- do.call(rbind.data.frame,
                     lapply(1994:2003,
                            function(i) {
                                sisben_year = filter(sisben_small, sisben_small$year == i)
                                out = manipulation_test(
                                    sisben_year$score,
                                    cutpoint = 47,
                                    bw = NULL,
                                    bin = NULL,
                                    h = 0.5,
                                    binning_proced = "camacho_conover_2011"
                                )
                                result = cbind(out$theta, out$se, out$p)
                                return(result)
                            }))
rownames(output_bw) = 1994:2003
colnames(output_bw) = c("Estimator_McCrary_undersmoothed",
                        "SE_McCrary_undersmoothed",
                        "p-value_undersmoothed")

# Uncomment the following lines to generate latex table
#table_bw = xtable(output_bw, digits=3)
#print.xtable(table_bw, type = "latex")
