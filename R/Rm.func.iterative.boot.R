#' Identify non-replicable studies in a meta-analysis from \eqn{R_m} (the parametric bootstrap)
#'
#' This function identifies non-replicable studies from \eqn{R_m} and replicability test (using the parametric bootstrap).
#' The meta-analysis dataset needs be transformed by \code{\link{to.dat.repMeta}}.
#'
#' @param dat1 data.frame, a meta-analysis dataset with \eqn{n} studies after transformed by \code{\link{to.dat.repMeta}}.
#' @param m numeric, \eqn{m} value, the default is 1.
#' @return A list object containing following components:
#' \describe{
#'   \item{\code{stat}}{a vector of the \eqn{R_m} value and its corresponding \eqn{p}-value after non-replicable studies are removed, using parametric bootstrap.}
#'   \item{\code{P_1}}{\eqn{p}-value of replicability test for the dataset using parametric bootstrap.}
#'   \item{\code{out_studies}{Row indices of the non-replicable studies in \code{dat1}}}
#' }
#' @examples
#' # Identify the non-replicable study using \eqn{R_1}
#' data.case <- to.dat.repMeta(data=moller2012,ai = r1, n1i = n1, ci = r2, n2i = n2,measure="OR")
#' iden.ls.boot <- Rm.func.iterative.boot(data.case,1)
#' # \eqn{p}-value of replicability test
#' iden.ls.boot$P_1
#' # index of non-replicable study
#' iden.ls.boot$out_studies
#' nonrep.id <- iden.ls.boot$out_studies
#' # \eqn{R_1} and its \eqn{p}-value among replicable studies
#' iden.ls.boot$stat
#' @export
Rm.func.iterative.boot <- function(dat1,m=1){
  # function to detect and identify non-replicable studies for bootstrap method
  Rls <- calR(dat1,m)
  Rls$empdist <- calempF(dat1,m)
  Pval <- (sum(Rls$empdist>Rls$maxR)+1)/(length(Rls$empdist)+1)
  Pval_leave1 <- Pval
  dat1.copy <- dat1
  out.n <- 1
  while(Pval<0.05 & out.n < (nrow(dat1.copy))){
    top2.ind <- order(Rls$Rvals,decreasing = T)[1:2]
    R.top2 <- sapply(top2.ind,function(x){
      dat1.sub <- dat1[-x,]
      Rls <- calR(dat1,1)
      return(Rls$maxR)
    })
    sm.ind <- which.min(R.top2)
    dat1.sub <- dat1[-top2.ind[sm.ind],]
    Rls <- calR(dat1.sub,1)
    Rls$empdist <- calempF(dat1,m)
    Pval <- (sum(Rls$empdist>Rls$maxR)+1)/(length(Rls$empdist)+1)
    out.n <- out.n+1*(Pval<0.05)
    dat1 <- dat1.sub
  }
  out.stdies <- which(!dat1.copy$y %in% dat1$y)
  return(list(P_1=Pval_leave1,out_studies=out.stdies,stat=c(Rls$maxR,Pval)))
}
