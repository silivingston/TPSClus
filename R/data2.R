#' A data frame of simulated multivariate time series data for predictions
#'
#' A data frame of simulated multivariate time series data for 30 individuals over one year.
#'
#' @format A data frame with 325 observations of 7 variables:
#' \describe{
#' \item{`SubjectID`}{Identifier for each subject}
#' \item{`Time`}{Time variable ranging from 0 to 365.}
#' \item{`Var1`, `Var2`, `Var3`}{Variables to be smoothed and clustered on.}
#' \item{`x1`, `x2`}{Covariates for modeling the outcome.}
#' }
#' @name TS.sim.new
NULL
