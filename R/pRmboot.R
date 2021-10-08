#' Test the replicability of multiple studies in a meta-analysis from \eqn{R_m} (the parametric bootstrap)
#'
#' This function obtains the \eqn{p}-value of the replicability test (detect non-replicability), using the parametric bootstrap.
#' The meta-analysis dataset needs be transformed by \code{\link{data.trans.bin}} or \code{\link{data.trans.cont}}.
#'
#' @param dat1 data.frame, a meta-analysis dataset with \eqn{n} studies after transformed by \code{\link{data.trans.bin}} or \code{\link{data.trans.cont}}.
#' @param m numeric, \eqn{m} value.
#' @return The \eqn{p}-value of \eqn{R_m} from the parametric bootstrap.
#' @examples
#' # Obtain the R1
#' data.case <- data.trans.bin(data=moller12,ai = r1, n1i = n1, ci = r2, n2i = n2,measure="OR")
#' pRmboot(data.case,1)
#' @export
pRmboot <- function(dat1,m=1){
  Rls <- calR(dat1,m)
  Rls$empdist <- calempF(dat1,m)
  Pval <- (sum(Rls$empdist>Rls$maxR)+1)/(length(Rls$empdist)+1)
  return(Pval)
}
