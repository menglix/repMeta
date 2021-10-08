#' Identify non-replicable studies in a meta-analysis from \eqn{R_m} (the Gumbel approximation)
#'
#' This function identifies non-replicable studies from \eqn{R_m} and replicability test (using the Gumbel approximation).
#' The meta-analysis dataset needs be transformed by \code{\link{data.trans.bin}} or \code{\link{data.trans.cont}}.
#'
#' @param dat1 data.frame, a meta-analysis dataset with \eqn{n} studies after transformed by \code{\link{data.trans.bin}} or \code{\link{data.trans.cont}}.
#' @param m numeric, \eqn{m} value, the default is 1.
#' @return A list object containing following components:
#' \describe{
#'   \item{\code{stat}}{a vector of the \eqn{R_m} value and its corresponding \eqn{p}-value after non-replicable studies are removed.}
#'   \item{\code{P_1}}{\eqn{p}-value of replicability test for the dataset using the Gumbel approximation.}
#'   \item{\code{out_studies}{Row indices of the non-replicable studies in \code{dat1}}
#' }
#' @examples
#' # Identify the non-replicable study using \eqn{R_1}
#' data.case <- data.trans.bin(data=moller12,ai = r1, n1i = n1, ci = r2, n2i = n2,measure="OR")
#' iden.ls <- Rm.func.iterative(data.case,1)
#' # \eqn{p}-value of replicability test
#' iden.ls$P_1
#' # index of non-replicable study
#' iden.ls$out_studies
#' nonrep.id <- iden.ls$out_studies
#' # \eqn{R_1} and its \eqn{p}-value among replicable studies
#' iden.ls$stat
#' @export
Rm.func.iterative <- function(dat1,m=1){
  # function to detect and identify non-replicable study under Gumbel approximation
  Rls <- calR(dat1,m)
  lev_set <- combn(nrow(dat1),m)
  ncombn <- ncol(lev_set)
  cn <- 2
  dn=2*(log(ncombn)+(m/2-1)*log(log(ncombn))-log(gamma(m/2)))
  Pval <- evd::pgumbel((1/cn)*(m*Rls$maxR-dn),loc=0,lower.tail = F)
  Pval_leave1 <- Pval
  dat1.copy <- dat1
  out.n <- 1
  # while(Pval<0.05 & out.n < min(4,floor(nrow(dat1)/2))){
  while(Pval<0.05 & out.n < (nrow(dat1.copy))){
    top2.ind <- order(Rls$Rvals,decreasing = T)[1:2]
    R.top2 <- sapply(top2.ind,function(x){
      dat1.sub <- dat1[-x,]
      Rls <- calR(dat1,1)
      return(Rls$maxR)
    })
    sm.ind <- which.min(R.top2)
    dat1.sub <- dat1[-top2.ind[sm.ind],]
    Rls <- calR(dat1.sub,m)
    lev_set <- combn(nrow(dat1.sub),m)
    ncombn <- ncol(lev_set)
    cn <- 2
    dn=2*(log(ncombn)+(m/2-1)*log(log(ncombn))-log(gamma(m/2)))
    Pval <- evd::pgumbel((1/cn)*(m*Rls$maxR-dn),loc=0,lower.tail = F)
    out.n <- out.n+1*(Pval<0.05)
    dat1 <- dat1.sub
  }
  out.stdies <- which(!dat1.copy$y %in% dat1$y)
  return(list(P_1=Pval_leave1,out_studies=out.stdies,stat=c(Rls$maxR,Pval)))
}
