#' @keywords internal
"_PACKAGE"

#' The 'FAIRR' package.
#'
#' @description Functional Approximation of Impulse Responses on R.
#'
#' @name FAIRR-package
#' @useDynLib FAIRR, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' \itemize{
#'   \item Barnichon, R., & Matthes, C. (2018). Functional Approximation of Impulse Responses. Journal of Monetary Economics, 99, 41-55.
#'   \item Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#' }
NULL
