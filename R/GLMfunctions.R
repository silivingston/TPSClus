#' Fit a generalized linear model (GLM) using clusters as predictors
#'
#' `FKM.glm()` fits a generalized linear model (GLM) using clusters output from
#' [cluster.fitted()] or [cluster.coefs()] as predictors, along with additional covariates.
#'
#' `FKM.glm()` applies the [glm()] function to fit a generalized linear model
#' using clusters as predictors. Clusters are obtained using [cluster.fitted()]
#' or [cluster.coefs()], and the output object of class `FKM.TPS` is input into
#' the `FKM.glm()` function, along with a dataset containing the output variable
#' and additional covariates of interest. Clusters are included using the
#' "partial assignment" method that employs the degree of cluster membership for
#' each individual to account for uncertainty in the cluster assignment.
#'
#' @param FKM_object An object of class '`FKM.TPS`' output from [cluster.fitted()]
#'   or [cluster.coefs()].
#' @param data A data frame with the same subjects used for spline-fitting and
#'   clustering that includes an outcome variable of interest and optional
#'   covariates.
#' @param y Name of the outcome variable (e.g. `y="Death"`)
#' @param covariates A vector of covariates of interest to be included in the
#'   model.
#' @param refclus Numeric identification of the cluster to be used as the
#'   reference cluster. Default is cluster 1 (`refclus=1`). Use `refclus=0` to
#'   identify the noise cluster as the reference cluster.
#' @param family A description of the error distribution and link function to be
#'   used in the model.
#' @param ... Additional arguments for the [glm()] function.
#'
#' @returns An object of class '`FKM.glm`' containing the following components:
#' * `FKM_object ` The inputted object of class '`FKM.TPS`'.
#' * `model_data ` A data frame containing the variables used in the model, including degree of cluster membership.
#' * `formula ` The formula used in the model.
#' * `family ` The family call used in the model.
#' * `covariates` The covariates that were included.
#' * `model_full ` The GLM model using clusters as predictors and any additional
#'    covariates of interest.
#' * `model_noclusters ` The GLM model using the covariates of interest but no
#'    clusters.
#' * `anova ` ANOVA comparing the models with and without clusters as predictors.
#' * `anova_pval ` P-value for the ANOVA comparing the models with and without
#'    clusters as predictors.
#'
#' @export
#'
#' @seealso
#' The glm function: [glm()]
#'
#' @examples
#' data(TS.sim)
#'
#' fitsplines <- TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
#'      ID="SubjectID", knots_time=c(0, 91, 182, 273, 365), n_fit_times=10)
#'
#' clusters1 <- cluster.fitted(fitsplines, k=3, m=1.3, seed=12345, RS=5, noise=TRUE)
#'
#' model <- FKM.glm(clusters1, TS.sim, y="outcome", covariates=c("x1", "x2"),
#' family="binomial")
#' summary(model)
#' model$anova_pval
#'
FKM.glm <- function(FKM_object, data, y, covariates, refclus=1, family="gaussian", ...) {
  if (class(FKM_object) != "FKM.TPS") {stop('Please supply an object of class "FKM.TPS"')}
  if (missing(data)) {stop('Please supply dataset')}
  # Check covariates
  if (!missing(covariates)) {
    nocov <- 0; ncov <- length(covariates)
    for (i in 1:ncov) {
      if (!(covariates[i] %in% names(data))) {stop(paste("Covariate '", covariates[i],"' not in dataset", sep=""))}
    }
  }
  if (missing(covariates)) {nocov <- 1; ncov <- 0; covariates <- NULL}

  # Get cluster variables
  clustervars <- colnames(FKM_object$FKM_TPS_U)[3:(length(colnames(FKM_object$FKM_TPS_U))-1)]
  if (refclus==0) {refclus <- length(clustervars)}
  modelvars <- clustervars[-refclus]

  # Get dataset with one entry per person and cluster info
  data2 <- merge(FKM_object$IDmatch, data, by=FKM_object$TPSdata$ID)
  data2 <- data2 %>% dplyr::group_by(.data$Id2) %>% dplyr::filter(dplyr::row_number()==1)
  if (nocov==0) {data3 <- subset(data2, select=c(FKM_object$TPSdata$ID, y, covariates))}
  if (nocov==1) {data3 <- subset(data2, select=c(FKM_object$TPSdata$ID, y))}
  data3 <- merge(data3, FKM_object$FKM_TPS_U, by=FKM_object$TPSdata$ID)

  # Get formula for full model
  if (nocov==0) {allvars <- c(modelvars, covariates)}
  if (nocov==1) {allvars <- modelvars}
  x <- rep("0", times=length(allvars)) #filler
  for (i in 1:(length(allvars)-1)) { x[i] <- paste(allvars[i], " + ", sep="") }
  x[length(allvars)] <- allvars[length(allvars)]
  x <- paste(x, sep="", collapse="")
  f1 <- paste(y, " ~ ", x, sep="")

  # Get formula for model without clusters
  if (nocov==0) {
    x2 <- rep("0", times=length(covariates)) #filler
    for (i in 1:(length(covariates)-1)) { x2[i] <- paste(covariates[i], " + ", sep="") }
    x2[length(covariates)] <- covariates[length(covariates)]
    x2 <- paste(x2, sep="", collapse="")
    f2 <- paste(y, " ~ ", x2, sep="")
  }

  if (nocov==1) {f2 <- paste(y, " ~ 1", sep="")}

  m1 <- glm(formula=f1, family=family, data=data3)
  m2 <- glm(formula=f2, family=family, data=data3)
  test <- anova(m1, m2, test="Chisq")
  pval <- test$`Pr(>Chi)`[2]

  output <- list(FKM_object=FKM_object, model_data=data3, formula=f1,
                 family=family,
                 covariates=covariates,
                 model_full= m1, model_noclusters=m2,
                 anova=test, anova_pval=pval)
  class(output) <- "FKM.glm"
  return(output)
}



#'
#' @export
print.FKM.glm <- function(x, ...) {
  if (class(x) != "FKM.glm") {stop('Please supply an object of class "FKM.glm"')}
  cat("Full model:\n")
  cat("Formula (f1): ", x$formula, "\n")
  cat("Family:", x$family, "\n")
  print(x$model_full)
  cat("\n")
  cat("ANOVA chi-square p-value for significance of clusters in model:\n")
  cat(x$anova_pval)
  cat("\n")
}


#'
#' @export
summary.FKM.glm <- function(object, ...) {
  if (class(object) != "FKM.glm") {stop('Please supply an object of class "FKM.glm"')}
  cat("Full model:\n")
  cat("Formula (f1): ", object$formula, "\n")
  cat("Family:", object$family, "\n")
  print(summary(object$model_full))
  cat("\n")
  cat("ANOVA chi-square p-value for significance of clusters in model:\n")
  cat(object$anova_pval)
  cat("\n")
}



#' Predict the model outcome for a new set of data
#'
#' Predict the model outcome for a new set of data by first fitting
#' tensor-product splines to the new dataset and then identifying degree of
#' cluster membership for the previously identified clusters.
#'
#' @param object An object of class '`FKM.glm`' found by modeling an
#'   outcome based on cluster membership degrees.
#' @param newdata A new data containing the same variables clustered on and used
#'   as covariates in the model.
#' @param ... Additional arguments.
#'
#' @returns A data frame containing the degree of cluster membership for each
#'   individual in the new dataset and predicted outcome.
#' @export
#'
#' @examples
#' # Fit initial model
#' data(TS.sim)
#'
#' fitsplines <- TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
#'      ID="SubjectID", knots_time=c(0, 91, 182, 273, 365), n_fit_times=10)
#'
#' clusters1 <- cluster.fitted(fitsplines, k=3, m=1.3, seed=12345, RS=5, noise=TRUE)
#'
#' model <- FKM.glm(clusters1, TS.sim, y="outcome", covariates=c("x1", "x2"),
#' family="binomial")
#'
#' # Get new dataset
#' data(TS.sim.new)
#'
#' # Predict outcome for new dataset
#' predicted <- predict(model, TS.sim.new)
#'
predict.FKM.glm <- function(object, newdata, ...) {
  if (class(object) != "FKM.glm") {stop('Please supply an object of class "FKM.glm"')}
  FKM_predict <- predict.FKM.TPS(object$FKM_object, newdata)

  # Get dataset with one entry per person and cluster info
  data2 <- merge(FKM_predict$IDmatch, newdata, by=object$FKM_object$TPSdata$ID)
  data2 <- data2 %>% dplyr::group_by(.data$Id2) %>% dplyr::filter(dplyr::row_number()==1)
  data3 <- subset(data2, select=c(object$FKM_object$TPSdata$ID, object$covariates))
  data3 <- merge(data3, FKM_predict$predicted_U, by=object$FKM_object$TPSdata$ID)

  # predicted GLM
  predicted <- predict.glm(object$model_full, newdata=data3, type="response")
  FKM_glm_predict <- cbind(data3, predicted)

  return(FKM_glm_predict)
}


