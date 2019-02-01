
calc.hj.tilde <- function (filter, js) {

    hj <- dwt.all.filters(filter, max(js))

    sapply(js, function (j) hj[[j]]$wavelet / sqrt(2^j))    
}



modwt.gappy.var <- function (x, eta, wavelet, js,
                             hj=calc.hj.tilde(dwt.filter(wavelet), js)) {
    ## ======================================================================
    ## Purpose: Calculate the MODWT gappy variance for a time series
    ## 'x' with observation pattern 'eta' and wavelet filter of name
    ## 'wavelet' on the levels 'js'.  'hj' is a list of the wavelet
    ## filters at each level
    ##
    ## 'eta' is assumed to be the same length as 'x' -- in each
    ## element of 'eta': FALSE means that case is missing; TRUE that
    ## case is observed.
    ##
    ## Updated: pfc@stat.osu.edu, Jan 2019
    ## ======================================================================
    
    eta.binary <- as.numeric(eta)
    
    z <- x
    z[!eta] <- 0
    
    sapply(hj, function (h) mean(CalcYjtHat(z, eta.binary, h)))
}



sigma.hat.Slepian <- function (x, eta, wavelet, js, R=7,
                               hj=calc.hj.tilde(dwt.filter(wavelet), js),
                               CORES=1) {
    ## ======================================================================
    ## Purpose: Calculate the multitaper spectral estimate based
    ## estimate of Sigma for data 'x' with observation pattern 'eta'
    ## and wavelet filter of name 'wavelet' on the levels 'js'.  We
    ## use 'R' tapers and 'hj' is a list of the wavelet filters at
    ## each level.  If CORES>1 then calculate Sigma using parallel
    ## computing operations.
    ##
    ## 'eta' is assumed to be the same length as 'x' -- in each
    ## element of 'eta': FALSE means that case is missing; TRUE that
    ## case is observed.
    ##
    ## Updated: pfc@stat.osu.edu, Jan 2019
    ## ======================================================================


    eta.binary <- as.numeric(eta)

    z <- x
    z[!eta] <- 0

    K <- length(js)

    if (CORES==1) {
        
        ys <- lapply(1:K, function (k) {
            CalcYjtHat(z, eta, hj[[k]])
            })
        
    } else {
        
        library(parallel)
        ys <- mclapply(1:K, function (k) {
            CalcYjtHat(z, eta, hj[[k]])
        }, mc.cores=CORES)
    }

    
    Mj <- sapply(ys, length)
    
    lambda <- lapply(Mj, function (M) dpss.taper(M, R, (R+2)/2))
    
    res <- sapply(1:K, function (k) {

        Qr <- as.numeric(crossprod(lambda[[k]], ys[[k]]))
        lam.plus <- colSums(lambda[[k]])
        beta.hat <- sum(lam.plus * Qr) / sum(lam.plus^2)
        Qr - beta.hat * lam.plus
    })
    
    Sigma <- matrix(NA, K, K)
    
    for (k1 in 1:K)
        for (k2 in 1:k1)
            Sigma[k1, k2] <- Sigma[k2, k1] <- mean(res[,k1] * res[,k2])
    
    list(vars=sapply(ys, mean),
         Mj=Mj,
         js=js,
         Sigma=Sigma)
}
