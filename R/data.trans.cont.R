#' Renmae a meta-analysis dataset with continuous outcomes.
#'
#' This function renames a meta-analysis dataset with continuous outcomes to the dataset compatible of the replicability analysis.
#'
#' @param data data.frame, a meta-analysis dataset with \eqn{n} studies.
#' @param y column name in \code{data} to specify the study-specific effect sizes.
#' @param s2 column name in \code{data} to specify the variance of the study-specific effect sizes.
#' @return A data.frame object.
#' @examples
#' # Obtain the transformed data frame format
#' data.case.a1 <- data.trans.cont(data=lahart18,y = y, s2 = s2)
#' @export
#' @export
data.trans.cont <- function(data=data,y,s2){
  data$y <- y
  data$s2 <- s2
  return(data)
}
