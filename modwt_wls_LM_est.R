

## ======================================================================
## Code to estimate the LRD parameter using a least squares estimator
## based on the MODWT wavelet variance.

## If eta=NULL, the data is gap free and we estimate Sigma using sums
## of squares; otherwise a multitaper spectral estimation method is
## used.
##
## Note: Below Abry means the diagonal LRD estimator.
##
## Last updated by pfc@stat.osu.edu, Jan, 2019.
## ======================================================================

modwt.wls.LM.est.prep <- function (x, wavelet, js, eta=NULL, R=7, CORES=1) {
    ## ======================================================================
    ## The helper function (arguments as 'modwt.wls.LM.est').
    ## ======================================================================
    
    if (is.null(eta)) { ## Use the gap free estimator
        
        dx  <- modwt(x, wavelet, max(js))
        v   <- modwt.var(dx)[js]

        K  <- length(js)        
        Sigma <- matrix(NA, K, K)
        for (k1 in 1:K) {
            for (k2 in 1:K) {
                j1 <- js[k1]
                j2 <- js[k2]
                Sigma[k1, k2] <- Sigma[k2, k1] <-
                    SigmaHat(dx$W[[j1]], dx$W[[j2]], dx$Mj[j1], dx$Mj[j2], dx$Lj[j1], dx$Lj[j2])
            }
        }

        list(vars=v,
             Mj=dx$Mj[js],
             js=js,
             Sigma=Sigma)

        
    } else { ## Use the gappy estimator
        
        filt <- dwt.filter(wavelet)
        hj   <- calc.hj.tilde(filt, js)
        sigma.hat.Slepian(x, eta, wavelet, js, R, hj, CORES)
    }
}



modwt.wls.LM.est <- function (x, wavelet, js, eta=NULL, R=7, CORES=1, include.full=TRUE,
                              prep=modwt.wls.LM.est.prep(x, wavelet, js, eta, R, CORES)) {
    ## ======================================================================
    ## Purpose: Estimate the LRD parameter using data 'x', a wavelet
    ## filter of name 'wavelet', on level 'js'.  'eta' is the missing
    ## data pattern: NULL means no missing data, otherwise 'eta' is
    ## assumed to be the same length as 'x' -- in each element of
    ## 'eta': FALSE means that case is missing; TRUE that case is
    ## observed.
    ## 
    ## If eta=NULL, the data is gap free and we estimate Sigma using
    ## sums of squares; otherwise a multitaper spectral estimation
    ## method is used -- we use 'R' tapers and if CORES>1 then we
    ## calculate Sigma using parallel computing operations.
    ##
    ## If 'include.full' is TRUE include the diagonal and full LRD
    ## estimators, otherwise only include the diagonal LRD estimator.
    ## ======================================================================


    z  <- log(prep$v)
    w  <- 2 * log(2) * prep$js
    K  <- length(z)

    D.full <- matrix(NA, K, K)
    for (k1 in 1:K) {
        for (k2 in 1:k1) {
            D.full[k1, k2] <- D.full[k2, k1] <-
                prep$Sigma[k1, k2] / (prep$v[k1] * prep$v[k2] * sqrt(prep$Mj[k1] * prep$Mj[k2]))
        }
    }

    dd <- diag(D.full)
    
    ## Abry/diagonal LRD estimator
    
    ddi <- 1/dd
    B1 <- sum(ddi)
    B2 <- sum(w * ddi)
    B3 <- sum(w^2 * ddi)
    
    at.A <- (B1 * w - B2) * ddi / (B1 * B3 - B2^2)
    
    delta.hat.abry <- 0.5 + sum(at.A * z)

    if (include.full) {

        ## Full LRD estimator

        D.full.inv <- solve(D.full)
        C  <- crossprod(w, D.full.inv)
        B1 <- sum(D.full.inv)
        B2 <- sum(C)
        B3 <- drop(C %*% w)
        
        at.F <- drop( crossprod(B1 * w - B2, D.full.inv) / (B1 * B3 - B2^2))
        
        delta.hat.full <- 0.5 + sum(at.F * z)
        
        list(abry=delta.hat.abry, full=delta.hat.full,
             v=prep$v, js=prep$js,
             z=z, D=D.full, Sigma=prep$Sigma, at.A=at.A, at.F=at.F, Mj=prep$Mj,
             se.abry=sqrt(crossprod(at.A, diag(diag(D.full))) %*% at.A),
             se.sand=sqrt(crossprod(at.A, D.full %*% at.A)),
             se.full=sqrt(crossprod(at.F, D.full) %*% at.F))
    } else {

        list(abry=delta.hat.abry, full=NA,
             v=prep$v, js=prep$js,
             z=z, D=D.full, Sigma=prep$Sigma, at.A=at.A, at.F=NA, Mj=prep$Mj,
             se.abry=sqrt(crossprod(at.A, diag(diag(D.full))) %*% at.A),
             se.sand=sqrt(crossprod(at.A, D.full %*% at.A)),
             se.full=NA)
    }
}


