#' Transform a meta-analysis dataset for replicability analysis.
#'
#' This function transforms a meta-analysis dataset to the dataset compatible of the package. More details are in \code{\link[metafor]{escalc}}.
#' @param measure a character string to specify which effect size or outcome measure should be calculated. See \code{Details} in \code{\link[metafor]{escalc}} for possible options and how the data needed to compute the selected effect size or outcome measure should then be specified (i.e., which of the following arguments need to be used).
#' @param ai vector with the \mjeqn{2 \times 2}{2x2} table frequencies (upper left cell).
#' @param bi vector with the \mjeqn{2 \times 2}{2x2} table frequencies (upper right cell).
#' @param ci vector with the \mjeqn{2 \times 2}{2x2} table frequencies (lower left cell).
#' @param di vector with the \mjeqn{2 \times 2}{2x2} table frequencies (lower right cell).
#' @param n1i vector with the group sizes or row totals (first group/row).
#' @param n2i vector with the group sizes or row totals (second group/row).
#' @param x1i vector with the number of events (first group).
#' @param x2i vector with the number of events (second group).
#' @param t1i vector with the total person-times (first group).
#' @param t2i vector with the total person-times (second group).
#' @param m1i vector with the means (first group or time point).
#' @param m2i vector with the means (second group or time point).
#' @param sd1i vector with the standard deviations (first group or time point).
#' @param sd2i vector with the standard deviations (second group or time point).
#' @param xi vector with the frequencies of the event of interest.
#' @param mi vector with the frequencies of the complement of the event of interest or the group means.
#' @param ri vector with the raw correlation coefficients.
#' @param ti vector with the total person-times or t-test statistics.
#' @param fi vector with the F-test statistics.
#' @param pi vector with the (signed) p-values.
#' @param sdi vector with the standard deviations.
#' @param r2i vector with the \mjseqn{R^2} values.
#' @param ni vector with the sample/group sizes.
#' @param yi vector with the observed effect sizes or outcomes.
#' @param vi vector with the corresponding sampling variances.
#' @param sei vector with the corresponding standard errors.
#' @param data optional data frame containing the variables given to the arguments above.
#' @param slab optional vector with labels for the studies.
#' @param subset optional (logical or numeric) vector to specify the subset of studies that will be included in the data frame returned by the function.
#' @param include optional (logical or numeric) vector to specify the subset of studies for which the measure should be calculated.
#' @param add a non-negative number to specify the amount to add to zero cells, counts, or frequencies. See \code{Details} in \code{\link[metafor]{escalc}}.}
#' @param to a character string to specify when the values under \code{add} should be added (either \code{"all"}, \code{"only0"}, \code{"if0all"}, or \code{"none"}). See \code{Details} in \code{\link[metafor]{escalc}}.}
#' @param drop00 logical to specify whether studies with no cases/events (or only cases) in both groups should be dropped when calculating the observed effect sizes or outcomes. See \code{Details} in \code{\link[metafor]{escalc}}.}
#' @param vtype a character string to specify the type of sampling variances to calculate. See \code{Details} in \code{\link[metafor]{escalc}}.}
#' @param var.names character vector with two elements to specify the name of the variable for the observed effect sizes or outcomes and the name of the variable for the corresponding sampling variances (the defaults are \code{"yi"} and \code{"vi"}).
#' @param add.measure logical to specify whether a variable should be added to the data frame (with default name \code{"measure"}) that indicates the type of outcome measure computed. When using this option, \code{var.names} can have a third element to change this variable name.
#' @param append logical to specify whether the data frame provided via the \code{data} argument should be returned together with the observed effect sizes or outcomes and corresponding sampling variances (the default is \code{TRUE}).
#' @param replace logical to specify whether existing values for \code{yi} and \code{vi} in the data frame should be replaced. Only relevant when \code{append=TRUE} and the data frame already contains the \code{yi} and \code{vi} variables. If \code{replace=TRUE} (the default), all of the existing values will be overwritten. If \code{replace=FALSE}, only \code{NA} values will be replaced.
#' @param digits optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is 4. Note that the values are stored without rounding in the returned object.
#' @param \dots other arguments in \code{\link[metafor]{escalc}}.
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
  data <- escalc(measure, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i,
                 t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, r2i, ni,
                 yi, vi, sei, data, slab, subset, include, add = 1/2, to = "only0",
                 drop00 = FALSE, vtype = "LS", var.names = c("yi", "vi"),
                 add.measure = FALSE, append = TRUE, replace = TRUE, digits,...)
  data$y <- data$yi
  data$s2 <- data$vi
  return(data)
}
