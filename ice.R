
## ======================================================================
## R code that accompanies the article
##
## P. F. Craigmile and D. Mondal
## 
## Estimation of long range dependence in gappy Gaussian time series
##
## NOTE: Requires the installation of the 'dwt' and 'spectral'
##       R libraries which are available from:
##
## https://github.com/petercraigmile/dwt
## https://github.com/petercraigmile/spectral
## ======================================================================

library(dwt)
library(Rcpp)
library(spectral)

sourceCpp("modwt_gappy_variance.cpp")

source("modwt_variance.R")
source("modwt_gappy_variance.R")
source("modwt_wls_LM_est.R")



## Read in the data, and set up the variables
df <- read.table("draft-profile.txt")

km <- seq(0, max(df[,1]))

profile <- rep(NA, length(km))
profile[df[,1]+1] <- df[,2]

rm(df)



## What percent of the data is missing?

cat(round(mean(is.na(profile)) * 100, 1), "% of the values are missing\n\n")


not.missing <- as.numeric(!is.na(profile))

js <- 1:7

ltaus <- log2(2^(js-1))

est.d <- modwt.wls.LM.est(profile, "D4", js, eta=not.missing)
pm <- 1.96*sqrt(diag(est.d$D))


print(round(c(est.d$abry,
              est.d$abry - 1.96 * est.d$se.abry,
              est.d$abry + 1.96 * est.d$se.abry), 2))

print(round(c(est.d$abry,
              est.d$abry - 1.96 * est.d$se.sand,
              est.d$abry + 1.96 * est.d$se.sand), 2))

print(round(c(est.d$full,
              est.d$full - 1.96 * est.d$se.full,
              est.d$full + 1.96 * est.d$se.full), 2))

est.d37 <- modwt.wls.LM.est(profile, "D4", 3:7, eta=not.missing)

model <- lm(est.d$z~ltaus, w=1/diag(est.d$D))



## Interpolate the data

ice <- profile
ice[is.na(ice)] <- -100

index <- which(ice >-100)
series <- ice[index]
y <- diff(series)/ diff(index)

z <- list()
for (i in 1: length(diff(index))) {
    z[[i]]  <-  rep(y[i], diff(index)[i])
}

v <- unlist(z)

interpolate.x <- diffinv(v) + ice[1]


## Refit the model

interp.est.d <- modwt.wls.LM.est(interpolate.x, "D4", js, eta=rep(TRUE,length(interpolate.x)))

interp.pm <- 1.96*sqrt(diag(interp.est.d$D))



## Produce Figure 1

dd <- 0.1

par(mfrow=c(1,2), cex=0.75, mar=c(3.5, 3.7, 1.5, 0.5), mgp=c(2.3,1,0), bty="L")

plot(km, profile, type="l",
     xlab="kilometers", ylab="draft profile (meters)")
mtext(side=3, "(a)", line=0.5, cex=0.8)

plot(ltaus+dd, est.d$z,
     xlab=expression(log(tau[j])), ylab=expression(log(nu[j]^2)),
     ylim=c(-5.5, -2), xlim=c(-0.2, 6.2))
mtext(side=3, "(b)", line=0.5, cex=0.8)

points(ltaus-dd, interp.est.d$z, pch=2, col="gray")

segments(ltaus+dd, est.d$z-pm, ltaus+dd, est.d$z+pm)

segments(ltaus-dd, interp.est.d$z-interp.pm,
         ltaus-dd, interp.est.d$z+interp.pm, col="gray")

lines(ltaus, fitted(model))
