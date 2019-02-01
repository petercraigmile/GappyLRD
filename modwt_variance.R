
modwt.var <- function (dx, js=1:dx$nlevels, biased=FALSE) {
    ## ======================================================================
    ## Purpose: Calculates either the unbiased or biased sample
    ## wavelet variance based on the MODWT object 'dx' at levels 'js'.
    ## (equation (306b) or equation (306c) of Percival and Walden
    ## (2000)).
    ##
    ## History: pfc@stat.osu.edu, Apr 2002.
    ## Updated: pfc@stat.osu.edu, Apr 2010.
    ## ======================================================================
  
    if (biased)
        sapply(js, function(j) mean(dx$W[[j]]^2))
    else 
        sapply(js, function (j) mean(dx$W[[j]][dx$Lj[j]:dx$N]^2))
}



modwt.var.approx.ci <- function (dx, vars, js=1:dx$nlevels, alpha=0.05) {
  ## ======================================================================
  ## Purpose: Calculates an approximate 100(1-'alpha')% confidence
  ## interval for the variances 'vars' calculated from the MODWT
  ## object 'dx' at levels 'js', using Equations (313c) and (314c) of
  ## Percival and Walden (2000).  Notes: This is very approximate!
  ##          
  ## History: pfc@stat.osu.edu, May 2002.
  ## Updated: pfc@stat.osu.edu, Apr 2010.
  ## ======================================================================

  etas <- max(dx$Mj[js]/2^js,1)
  
  lower <- etas * vars / qchisq(1-alpha/2,1)
  upper <- etas * vars / qchisq(alpha/2,1)

  list(vars=vars, lower=lower, upper=upper)
}



