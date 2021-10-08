#' Individual format physical activity versus control for reducing body fat in women with breast cancer after adjuvant therapy
#'
#' A meta-analysis dataset compared individual format physical activity versus control for reducing body fat in women with breast cancer after adjuvant therapy
#' in 10 studies, and the effect measure is the mean difference.
#'
#' @format A data frame with 10 rows and 7 variables:
#' \describe{
#'   \item{\code{CDSR.id}}{character, CDSR ID}
#'   \item{\code{data.type}}{character, show whether the outcome is continous (CONT) or binary (DICH)}
#'   \item{\code{study.name}{character, author and year}}
#'   \item{\code{y}{numeric, effect size}}
#'   \item{\code{s2}{numeric, variance of the effect size}}
#'   \item{\code{n1}{numeric, number of subjects in the treatment group}}
#'   \item{\code{n2}{numeric, number of subjects in the control group}}
#' }
#' @source \url{https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.CD011292.pub2/full}
"lahart18"

#' off-pump with heparin dose less than 300 mg/k vs. on-pump coronary artery bypass grafting (CABG) for reducing postoperative atrial fibrillation (POAF)
#'
#' A meta-analysis dataset compared off-pump with heparin dose less than 300 mg/k vs. on-pump coronary artery bypass grafting (CABG)
#' for reducing postoperative atrial fibrillation (POAF) in 17 studies, the effect measure we used is the odds ratio (OR).
#'
#' @format A data frame with 17 rows and 7 variables:
#' \describe{
#'   \item{\code{CDSR.id}}{character, CDSR ID}
#'   \item{\code{data.type}}{character, show whether the outcome is continous (CONT) or binary (DICH)}
#'   \item{\code{study.name}{character, author and year}}
#'   \item{\code{r1}{numeric, event counts in the treatment group}}
#'   \item{\code{r2}{numeric, event counts in the control group}}
#'   \item{\code{n1}{numeric, number of subjects in the treatment group}}
#'   \item{\code{n2}{numeric, number of subjects in the control group}}
#' }
#' @source \url{https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.CD007224.pub2/full}
"moller12"

