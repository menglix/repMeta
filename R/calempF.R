bootdat <- function(dat1,SEED,tau2){
  set.seed(SEED)
  dat1$s2 <- sample(dat1$s2,nrow(dat1),replace=TRUE)
  w_as <- 1/(dat1$s2+tau2)
  theta_r <- weighted.mean(dat1$y,w_as)
  dat1.boot <- data.frame(y=rnorm(nrow(dat1), theta_r, sqrt(dat1$s2 + tau2)), s2=dat1$s2)
  return(dat1.boot)
}
#' Obtain multiple \eqn{R_m} values from parametric bootstrap
#'
#' This function obtains numbers of \eqn{nb} \eqn{R_m} values from parametric bootstrap
#'
#' @param dat1 data frame, a meta-analysis dataset with after \code{\link{to.dat.repMeta}}.
#' @param m numeric, \eqn{m} value.
#' @param nb numeric, number of iterations in parametric bootstrap.
#' @return \eqn{R_m} values from parametric bootstrap, a vector of length \eqn{nb}.
#' @examples
#' # 1000 bootstrapped \eqn{R_1} values
#' calemF(dat1,1)
#' # 500 bootstrapped \eqn{R_1} values
#' calemF(dat1,1,nb=500)
#' @export
calempF <- function(dat1,m,nb=1000){
  set.seed(9876)
  tau2 <- meta.reml(dat1)$tau2
  ## 1000 bootstrap runs
  seed.ls <- sample(1:200000,nb)
  Rmls <- parallel::mclapply(seed.ls,function(SEED){
    dat1.boot <- bootdat(dat1,SEED,tau2)
    ## leave one study out
    maxR <- calR(dat1.boot,m)
    Rm <- maxR$maxR
    return(Rm)
  },mc.cores=6)
  ## use bootstrap to estimate the sample variance
  Rval2 <- data.frame(do.call(cbind,Rmls))
  Rval2 <- as.vector(Rval2)
  return(Rval2)
}
