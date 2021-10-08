#' Transform a meta-analysis dataset with binary outcomes.
#'
#' This function transforms a meta-analysis dataset with binary outcomes to the dataset compatible of the package.
#'
#' @param data data.frame, a meta-analysis dataset with \eqn{n} studies.
#' @param ai column name in \code{data} to specify the 2 × 2 table frequencies (upper left cell).
#' @param bi column name in \code{data} to specify the 2 × 2 table frequencies (upper right cell).
#' @param ci column name in \code{data} to specify the 2 × 2 table frequencies (lower left cell).
#' @param di column name in \code{data} to specify the 2 × 2 table frequencies (lower right cell).
#' @param n1i column name in \code{data}  to specify the group sizes or row totals (ﬁrst group/row).
#' @param n2i column name in \code{data}  to specify the group sizes or row totals (second group/row).
#' @param measure a character string indicating which effect size or outcome measure should be calculated. See \code{\link[metafor]{escalc}} for details.
#' @param ... other arguments in \code{\link[metafor]{escalc}}.
#' @return An object of class c("escalc","data.frame").
#' @examples
#' # Obtain the transformed data frame format with odds ratio as the summary measure
#' data.case <- data.trans.bin(data=moller12,ai = r1, n1i = n1, ci = r2, n2i = n2,measure="OR")
#' @export
data.trans.bin <- function(data=data,n1i,ai,bi,n2i,ci,di,measure,...){
  mf <- match.call()
  mf.ai <- mf[[match("ai", names(mf))]]
  mf.bi <- mf[[match("bi", names(mf))]]
  mf.ci <- mf[[match("ci", names(mf))]]
  mf.di <- mf[[match("di", names(mf))]]
  mf.n1i <- mf[[match("n1i", names(mf))]]
  mf.n2i <- mf[[match("n2i", names(mf))]]
  ai <- eval(mf.ai, data, enclos = sys.frame(sys.parent()))
  bi <- eval(mf.bi, data, enclos = sys.frame(sys.parent()))
  ci <- eval(mf.ci, data, enclos = sys.frame(sys.parent()))
  di <- eval(mf.di, data, enclos = sys.frame(sys.parent()))
  n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
  n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
  if (is.null(bi))
    bi <- n1i - ai
  if (is.null(di))
    di <- n2i - ci
  data <- escalc(measure=measure,ai = ai, n1i = n1i, ci = ci, n2i = n2i,data=data,...)
  data$y <- data$yi
  data$s2 <- data$vi
  return(data)
}
