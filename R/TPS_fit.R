#'Fit tensor product splines to longitudinal data
#'
#'`TPSfit()` is used to fit multidimensional tensor product splines to
#'longitudinal data with three or more variable of interest prior to
#'implementation of a clustering algorithm.
#'
#'`TPSfit()` employs package `mgcv` to fit a tensor product splines to each
#'individual using a generalized additive model. The fitted splines are
#'two-dimensional, with one dimension being the variable identifier and the
#'other being time. An adequate number of observed time points are required for
#'each individual, and the number of knots should be less than the smallest
#'number of time points. If splines are unable to be fit for an individual, an
#'error message will be shown, but splines will be fit for remaining
#'individuals. A vector of identifiers for individuals with errors is included
#'in the output as `error_subjects`, and these subjects are not included in the
#'output `GAMSsfitted` or `GAMscoef`.
#'
#'
#' @param data A longitudinal dataset in long form with multiple variables
#'   measured over time.
#' @param time Name of the time variable (e.g. "Time").
#' @param vars A character vector of at least 3 variables of interest.
#' @param ID Name of the subject ID variable.
#' @param knots_time A numeric vector of knots for spline-fitting the time
#'   variable. Must supply knots_time or kt.
#' @param kt Number of evenly spaced knots for spline-fitting the time variable
#'   if knots_time is not given.
#' @param fit_times Optional vector for times where fitted values will be
#'   calculated. If fit_times and n_fit_times are not given, fitted values are
#'   calculated at knots.
#' @param n_fit_times Number of evenly spaced times where fitted values will be
#'   calculated if fit_times are not given.
#' @param st Logical expression indicating whether each variable should be
#'   standardized.
#'
#' @returns An object of class '`TPSfit`' containing the following components: \cr
#'* `GAMsfitted` A data frame containing the fitted spline values. \cr
#'* `GAMscoef` A data frame containing the tensor product spline coefficients \cr
#'* `data_long` A data frame containing data in long format for both time and
#'   variable \cr
#'* `knots` A list of two vectors containing the variable and time knots \cr
#'* `indiv_means` A list containing a data frame of individual means for each of
#'   the variables of interest \cr
#'* `GAMs` A list containing the generalized additive models for fitting splines
#'   on each individual \cr
#'* `nsubject` The number of subjects in the dataset \cr
#'* `IDmatch` A data frame matching the original subject ID and new consecutive
#'   ID numbers \cr
#'* `error_subjects` A vector of individuals that encountered errors in the
#'   spline-fitting process \cr
#'
#' @export
#'
#' @seealso The `mgcv` R package: <https://cran.r-project.org/web/packages/mgcv/index.html>
#'
#' @examples
#' library(tidyr); library(dplyr); library(mgcv)
#' data(TS.sim)
#'
#' fitsplines <- TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
#'      ID="SubjectID", knots_time=c(0, 91, 182, 273, 365), n_fit_times=10)
#'
#' fitsplines2 <- TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"),
#' time="Time", ID="SubjectID", knots_time=c(0, 91, 182, 273, 365),
#'      fit_times=c(46, 91, 137, 182, 228, 273, 319))
#'
TPSfit <- function(data, time, vars, ID, knots_time, kt, fit_times,
                    n_fit_times, st=TRUE) {
  # Error messages
  if (!is.data.frame(data)) {stop('Data must be a data frame.')}
  if (missing(time)) {stop('Please supply name of the time variable in the form time="your variable name". ')}
  if (missing(vars)) {stop('Please supply names of variables to smooth in the form vars=c("var1","var2", etc). ')}
  if (missing(ID)) {stop('Please supply name of the ID variable in the form ID="your variable name". ')}
  if (missing(knots_time) & missing(kt)) {stop('Please supply knots_time, a vector for placement of knots for the time variable, or
kt, the number of knots to be used for the time variable.')}
  # Rename variables of interest with Var1, Var2, etc
  numvars <- length(vars)
  if (numvars < 3) {stop('Please include at least three variables.')}
  varnumbers <- 1:numvars
  varnames <- paste("Var", varnumbers, sep="")
  Vars <- data[,vars]
  colnames(Vars) <- varnames
  Time <- data[,time]
  Id <- data[,ID]
  Id2 <- as.numeric(factor(data[,ID]))
  uId <- unique(Id)
  uId2 <- unique(Id2)
  IDmatch <- data.frame(uId, uId2)
  colnames(IDmatch)[1] <- ID; colnames(IDmatch)[2] <- "Id2"
  dat <- data.frame(Id2, Time, Vars)

  # Pivot to long Variable, requires tidyr
  data_long <- dat %>% tidyr::pivot_longer(
    cols=all_of(varnames),
    names_to="Variable",
    names_prefix="Var",
    values_to="Value"
  )
  data_long <- as.data.frame(data_long)
  data_long$Variable <- as.numeric(data_long$Variable)

  # Standardize data
  if (st==TRUE) {
    stdata <- standardize_vars(data_long)
    data_long <- stdata$stdata
  }
  if (st==TRUE) {data_long$x <- data_long$ValueS}
  if (st==FALSE) {data_long$x <- data_long$Value}

  # Create knots for spline-fitting
  knots_var <- varnumbers
  if (missing(knots_time)) {
    min_time <- min(data_long$Time)
    max_time <- max(data_long$Time)
    knot_dist <- (max_time-min_time)/(kt-1)
    knots_time <- min_time
    for (k in 2:kt) {knots_time <- append(knots_time, knot_dist*(k-1))}
  }
  knots <- list(knots_time, knots_var)

  # Times to calculate fitted values
  if (missing(fit_times) & missing(n_fit_times)) {fit_times <- knots_time}
  if (missing(fit_times) & !missing(n_fit_times)) {
    min_time <- min(data_long$Time)
    max_time <- max(data_long$Time)
    fit_dist <- (max_time-min_time)/(n_fit_times-1)
    fit_times <- min_time
    for (k in 2:n_fit_times) {fit_times <- append(fit_times, fit_dist*(k-1))}
  }

  # Set up data for fitting values
  Time <- rep(fit_times, times=numvars)
  Variable <- rep(varnumbers, each=length(fit_times))
  newdat <- data.frame(Time, Variable)
  subjects <- unique(dat$Id2)
  nsubjects <- length(subjects)

  # Lists to store fitted splines
  GAMs <- list()
  GAMfit <- list()
  GAMcoefs <- list()
  error_subjects <- NULL

  for (i in 1:nsubjects) {
    datasub <- dplyr::filter(data_long, Id2==i)
    tryCatch(
      {
        GAMsub <- mgcv::gam(x ~ te(Time, Variable, k=c(length(knots_time), length(knots_var))),
                      knots=knots, data=datasub)
        Id2 <- rep(i, times=length(fit_times)*numvars)
        newdata <- cbind(Id2, newdat)
        newdata$x <- mgcv::predict.gam(GAMsub, newdata)
        GAMfit[[i]] <- newdata
        GAMbasis <- GAMsub$coefficients
        Id2 <- rep(i, times=length(GAMbasis))
        GAMcoef <- as.data.frame(cbind(Id2, GAMbasis))
        GAMcoefs[[i]] <- GAMcoef
        GAMs[[i]] <- GAMsub
      },
      error=function(e){
        message(paste("Subject removed due to error in spline-fitting:", i))
        print(e)
        #error_subjects <<- i
        ifelse(is.null(error_subjects), error_subjects <<- i, error_subjects <<- c(error_subjects, i))
      },
      warning=function(w){
        message(paste("warning occured for subject", i))
        print(w)
      }
    )
  }

  GAMsfitted <- dplyr::bind_rows(GAMfit)
  GAMscoef <- dplyr::bind_rows(GAMcoefs)

  GAMsfitted$t <- rep(1:length(fit_times), times=(numvars*(nsubjects-length(error_subjects))))
  GAMsfitted <- dplyr::arrange(GAMsfitted, Id2, Time, Variable)

  # Center of individual mean
  center <- center_indiv(GAMsfitted)
  GAMsfitted <- cbind(GAMsfitted, center$object2)
  names(GAMsfitted)[names(GAMsfitted) == 'Time'] <- 'FitTime'

  # Dataset with spline coefficients
  GAMscoef$coef <- rep(1:(length(GAMbasis)), times=(nsubjects-length(error_subjects)))

  GAMsfitted <- merge(IDmatch, GAMsfitted, by="Id2")
  GAMscoef <- merge(IDmatch, GAMscoef, by="Id2")
  output <- list(GAMsfitted=GAMsfitted, GAMscoef=GAMscoef, fit_times=fit_times, vars=vars,
                 data_long=data_long, knots=knots, indiv_means=center$indiv_means, GAMs=GAMs,
                 nsubjects=nsubjects, ID=ID, IDmatch=IDmatch, error_subjects=error_subjects, call=match.call())
  class(output) <- "TPSfit"
  return(output)
}

#' Internal function to standardize variables
#'
#' @param data data
#'
#'
#' @noRd
standardize_vars <- function(data) {
  varid <- as.numeric(unique(data$Variable))
  data_by_var <- list()
  for (i in varid) {data_by_var[[i]] <- dplyr::filter(data, .data$Variable==i)}
  var_means <- rep(NA, times=length(varid))
  for (i in varid) {var_means[i] <- mean(data_by_var[[i]]$Value, na.rm=TRUE)}
  var_sd <- rep(NA, times=length(varid))
  for (i in varid) {var_sd[i] <- sd(data_by_var[[i]]$Value, na.rm=TRUE)}
  for (i in varid) {data_by_var[[i]]$ValueS <- (data_by_var[[i]]$Value-var_means[i])/var_sd[i]}

  data2 <- dplyr::bind_rows(data_by_var)
  data2 <- dplyr::arrange(data2, .data$Id2, .data$Time)
  answer <- list(stdata=data2, var_means=var_means, var_sd=var_sd)
  return(answer)
}


#' Internal function to center on individual trajectories
#'
#' @param object object
#'
#' @noRd
#'
center_indiv <- function(object) {
  varid <- as.numeric(unique(object$Variable))
  data_by_var <- list()
  for (i in varid) {data_by_var[[i]] <- dplyr::filter(object, .data$Variable==i)}
  indiv_means <- list()
  for (i in varid) {
    indiv_means[[i]] <- aggregate(x~ Id2, data=data_by_var[[i]], mean)
    colnames(indiv_means[[i]]) <- c("Id2", "mean_x")
    data_by_var[[i]] <- merge(data_by_var[[i]], indiv_means[[i]], by="Id2")
    data_by_var[[i]]$centered_x <- data_by_var[[i]]$x - data_by_var[[i]]$mean_x
  }
  object2 <- dplyr::bind_rows(data_by_var)
  object2 <- dplyr::arrange(object2, .data$Id2, .data$Time, .data$Variable)
  answer <- list(object2=object2[6:7], indiv_means=indiv_means)
  return(answer)
}

#' Internal function
#'
#' @param x x
#' @param byrow logical
#'
#' @noRd
#'
unmatrixB <- function (x, byrow = FALSE)   #adapted from gdata
{
  rnames <- rownames(x)
  cnames <- colnames(x)
  if (is.null(rnames))
    rnames <- paste("S", 1:nrow(x), sep = "")
  if (is.null(cnames))
    cnames <- paste("V", 1:ncol(x), sep = "")
  nmat <- outer(rnames, cnames, paste, sep = ":")
  if (byrow) {
    vlist <- c(t(x))
    names(vlist) <- c(t(nmat))
  }
  else {
    vlist <- c(x)
    names(vlist) <- c(nmat)
  }
  return(vlist)
}


#' @export
#'
print.TPSfit <- function(x, ...) {
  if (class(x) != "TPSfit") {stop('Please supply an object of class "TPSfit"')}
  cat("Object of type 'TPSfit'\n")
  cat("\n")
  cat("Tensor-product splines fit for", (x$nsubjects-length(x$error_subjects)),
      "out of", x$nsubjects, "subjects\n", sep=" ")
  cat("Variables of interest:", x$vars, "\n", sep=" ")
  cat("Time knots:", x$knots[[1]], "\n",sep=" ")
  if(!is.null(x$error_subjects)) {
    (cat("Errors: Unable to fit splines for subjects (Id2): ", x$error_subjects, "\n", sep=" "))}
  cat("Output: GAMscoef contains model coefficients\n")
  cat("Output: GAMsfitted has fitted values at times:", x$fit_times, "\n", sep=" ")
  cat("\n")
  cat("Available components:\n")
  print(names(x))
}



#' @export
#'
summary.TPSfit <- function(object, ...) {
  if (class(object) != "TPSfit") {stop('Please supply an object of class "TPSfit"')}
  cat("Object of class 'TPSfit'\n")
  cat("\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")
  cat("Tensor-product splines fit for", (object$nsubjects-length(object$error_subjects)),
      "out of", object$nsubjects, "subjects\n", sep=" ")
  cat("Variables of interest:", object$vars, "\n", sep=" ")
  cat("Time knots:", object$knots[[1]], "\n",sep=" ")
  cat("\n")
  if(!is.null(object$error_subjects)) {
    (cat("Errors: Unable to fit splines for subjects (Id2): ", object$error_subjects, "\n", sep=" "))}
  cat("\n")
  cat("Output: GAMscoef contains model coefficients: (first 6 rows)\n")
  print(head(object$GAMscoef))
  cat("\n")
  cat("Output: GAMsfitted has fitted values at times:", object$fit_times, "\n", sep=" ")
  cat("First 6 rows of GAMsfitted:\n")
  print(head(object$GAMsfitted))
  cat("\n")
  cat("Available components:\n")
  print(names(object))
}
