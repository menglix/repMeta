
calI.cont <- function(dat1){
  weight <- 1/dat1$s2
  theta_f <- weighted.mean(dat1$y,weight)
  Q <- sum(weight*(dat1$y-theta_f)^2)
  H2 <- Q/(nrow(dat1)-1)
  I2 <- (Q>(nrow(dat1)-1))*(1-1/H2)
  return(I2)
}


meta.reml <- function(dat){
  # obtain reml estimator of tau
  return(
    tryCatch(
      {
        res3 <- metafor::rma(yi=y,vi=s2,data=dat,method='REML',control = list(stepadj=0.5,maxiter=1000))
        return(res3)
      },
      error=function(e) {
        dat$study.name <- c(1:nrow(dat))
        res3 <- metafor::rma.mv(y, s2, random = ~ 1 | study.name, data=dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
        res3$tau2 <- res3$sigma2
        return(res3)}
    ))
}
#' Quantify replicability in a meta-analysis by the observed \eqn{R_m} statistic
#'
#' This function obtains the observed \eqn{R_m} statistic to qunaitfy replicability of all studies in a meta-analysis.
#' The meta-analysis dataset needs be transformed by \code{\link{data.trans.bin}} or \code{\link{data.trans.cont}}.
#'
#' @param dat1 data.frame, a meta-analysis dataset with \eqn{n} studies after transformed by \code{\link{data.trans.bin}} or \code{\link{data.trans.cont}}.
#' @param m numeric, \eqn{m} value.
#' @return A list contains the observed Rm and \eqn{\mathcal{C}^n_m} values of \eqn{R_{\mathcal{A}_{m,k}}} from the meta-analysis.
#' \describe{
#'   \item{\code{maxR}}{the observed Rm.}
#'   \item{\code{Rvals}}{\eqn{R_{\mathcal{A}_{m,k}}} values, a vector of length \eqn{\mathcal{C}^n_m}.}
#' }
#' @examples
#' # Obtain the R1
#' data.case <- data.trans.bin(data=moller12,ai = r1, n1i = n1, ci = r2, n2i = n2,measure="OR")
#' calR(data.case,1)
#' @export
calR <- function(dat1,m){
  # obtain Rm statistic
  lev_set <- combn(nrow(dat1),m)
  ncombn <- ncol(lev_set)
  Rval_ls <- sapply(c(1:ncombn),function(c){
    dat2 <- dat1[-lev_set[,c],]
    weight <- 1/dat2$s2
    theta_f <- weighted.mean(dat2$y,weight)
    res <- meta.reml(dat2)
    tau2 <- res$tau2
    var.tr <- res$se^2
    theta_r <- res$beta
    dat3 <- dat1[lev_set[,c],]
    d_i <- dat3$y-theta_r
    var_di <- dat3$s2+tau2+var.tr
    ri <- d_i/sqrt(var_di)
    R_value <- sum(ri^2)/m
    return(R_value)
  })
  return(list(maxR=max(Rval_ls),Rvals=Rval_ls))
}
