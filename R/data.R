#' A data frame of simulated multivariate time series data
#'
#' A data frame of simulated multivariate time series data for 150 individuals over one year.
#'
#' @format A data frame with 1919 observations of 8 variables:
#' \describe{
#' \item{`SubjectID`}{Identifier for each subject}
#' \item{`Time`}{Time variable ranging from 0 to 365.}
#' \item{`Var1`, `Var2`, `Var3`}{Variables to be smoothed and clustered on.}
#' \item{`x1`, `x2`}{Covariates for modeling the outcome.}
#' \item{`outcome`}{Binary outcome variable.}
#' }
#' @name TS.sim
NULL

