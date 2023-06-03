#' Test the replicability of multiple studies in a meta-analysis from \eqn{R_m} (the Gumbel approximation)
#'
#' This function obtains the \eqn{p}-value of the replicability test (detect non-replicability), using the Gumbel approximation.
#' The meta-analysis dataset needs be transformed by \code{\link{to.dat.repMeta}}.
#'
#' @param dat1 data.frame, a meta-analysis dataset with \eqn{n} studies after transformed by \code{\link{to.dat.repMeta}}.
#' @param m numeric, \eqn{m} value.
#' @return The \eqn{p}-value of \eqn{R_m} from the Gumbel approximation.
#' @examples
#' # Obtain the R1
#' data.case <- to.dat.repMeta(data=moller12,ai = r1, n1i = n1, ci = r2, n2i = n2,measure="OR")
#' pRmasym(data.case,1)
#' @export
pRmasym <- function(dat1,m){
  Rls <- calR(dat1,m)
  lev_set <- combn(nrow(dat1),m)
  ncombn <- ncol(lev_set)
  cn <- 2
  dn=2*(log(ncombn)+(m/2-1)*log(log(ncombn))-log(gamma(m/2)))
  Pval <- evd::pgumbel((1/cn)*(m*Rls$maxR-dn),loc=0,lower.tail = F)
  return(Pval)
}
