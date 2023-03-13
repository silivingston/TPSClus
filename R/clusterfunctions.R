
#' Fuzzy c-medoids clustering on fitted values of splines
#'
#' `cluster.fitted()` is used to apply a fuzzy c-medoids clustering algorithm to
#'   the fitted values that are output from the `TPSfit()` function.
#'
#' `cluster.fitted()` employs package \pkg{fclust} to apply the fuzzy c-medoids
#'   clustering algorithm to the fitted tensor-product spline values that are
#'   output from the [TPSfit()] function. The user has options of including a noise
#'   cluster and centering trajectories on the individual means. The `addslopes`
#'   option will cluster on the slopes between fitted time-points in addition to
#'   cross-sectional fitting at those points.
#'
#'
#' @param TPSdata Output from the [TPSfit()] function.
#' @param center A logical expression indicating that trajectories are to be
#'   centered on individual means.
#' @param addslopes A logical expression indicating that clustering will be
#'   conducted on the slopes between fitted time-points in addition to the
#'   cross-sectional values at each point.
#' @param k An integer value indicating the number of clusters.
#' @param m A numeric value called the 'fuzziness parameter.' Must be greater
#'   than 1, with a default of 1.2.
#' @param noise A logical expression indicating whether a noise cluster for
#'   outliers is to be included.
#' @param seed An optional numeric value for the seed, as the clustering
#'   algorithm uses random start values.
#' @param ... Additional optional values employed by \pkg{fclust}, including `RS`,
#'   the number of random starts, and `delta`, the distance for outliers.
#'
#' @returns An object of class '`FKM.TPS`' containing the following components:
#'* `FKM_TPS `  An object output by the `fclust` algorithm.
#'* `U `   A dataframe containing the degrees of cluster membership for each subject.
#'* `Umax `  A vector containing the modal cluster assignment for each subject.
#'* `FKM_TPS_U `  A data frame containing subjects with subject ID, degrees of cluster
#'   membership, and modal class assignment.
#'* `FKM_indices `  A vector containing the values for six cluster validity indices.
#'   See [fclust::Fclust.index()]
#'* `wide_data `  A data frame with the values that were clustered on in the
#'   clustering algorithm plus subject IDs.
#'* `TPSdata `  The '`TPSfit`' object that was input into the function.
#'* `center `  Logical value that was used to determine whether to center on individual
#'   trajectories.
#'* `noise `  Logical value that was used to determine whether a noise cluster is
#'   included.
#'* `k `  The number of clusters.
#'* `m `  The fuzziness parameter.
#'* `addslopes `  Logical value that was used to determine whether to include the
#'   slopes between fitted time-points.
#'* `IDmatch `  A data frame matching the original subject ID with new consecutive
#'   ID values.
#'
#'
#' @export
#'
#' @seealso
#' * [cluster.coefs()] is an alternative method that clusters on spline
#'   coefficients rather than fitted values.
#' * The \pkg{fclust} R package: <https://CRAN.R-project.org/package=fclust>
#' * Fuzzy k-medoids function: [fclust::FKM.med()]
#' * Fuzzy k-medoids function with noise cluster: [fclust::FKM.med.noise()]
#'
#'
#' @examples
#' data(TS.sim)
#'
#' fitsplines <- TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
#'      ID="SubjectID", knots_time=c(0, 91, 182, 273, 365), n_fit_times=10)
#'
#' clusters1 <- cluster.fitted(fitsplines, k=3, m=1.3, seed=12345, RS=5, noise=TRUE)
#' summary(clusters1)
#'
cluster.fitted <- function(TPSdata, center=TRUE, addslopes=TRUE, k, m=1.2,
                           noise=TRUE, seed, ...) {
  GAMsfitted <- TPSdata$GAMsfitted

  # Get correct values to cluster on (centered or not)
  if (center==TRUE) {GAMsfitted$X <- GAMsfitted$centered_x}
  if (center==FALSE) {GAMsfitted$X <- GAMsfitted$x}
  if (center!= TRUE & center!=FALSE) {stop('center must be TRUE/FALSE')}

  # Convert data from long to wide
  #long_data <- subset(GAMsfitted, select=c(Id2, Variable, t, X))
  long_data <- GAMsfitted[,c(1,4,6,9)]
  long_data <- dplyr::arrange(long_data, .data$Id2, .data$Variable, .data$t)
  wide_data <- as.data.frame(tidyr::pivot_wider(long_data, names_from=c(.data$Variable, .data$t), names_prefix="X", values_from=.data$X))

  # Add slopes:
  if (addslopes==TRUE) {
    npoints <- length(TPSdata$fit_times)
    nvars <- length(TPSdata$vars)
    wide_data_m <- as.matrix(wide_data)
    wide_slopes_m <- matrix(0, nrow=length(wide_data$Id2), ncol=nvars*(npoints-1))
    for (j in 1:length(wide_data$Id2)) {    #j in wide_data$Id2
      subjmat <- matrix(wide_data_m[j,-1], nrow=npoints, ncol=nvars)
      slopemat <- matrix(0, nrow=npoints-1, ncol=nvars)
      for (i in 2:npoints-1) {
        for (q in 1:nvars) {
          slopemat[i,q] <- subjmat[i+1,q]-subjmat[i,q]
        }
      }
      wide_slopes_m[j,] <- unmatrixB(slopemat)
    }
    wide_slopes <- as.data.frame(wide_slopes_m)
    colnames(wide_slopes) <- names(unmatrixB(slopemat))
    wide_data <- cbind(wide_data, wide_slopes)
  }

  # Remove ID
  cluster_data <- wide_data[,-1]

  # Clustering using fclust package
  if (!missing(seed)) {set.seed(seed)}
  if (noise==TRUE) {
    FKM_TPS <- fclust::FKM.med.noise(cluster_data, k=k, m=m, ...)
  }
  if (noise==FALSE) {
    FKM_TPS <- fclust::FKM.med(cluster_data, k=k, m=m, ...)
  }
  if (noise!=TRUE & noise!=FALSE) {stop('noise should be TRUE/FALSE')}

  # Dataset with degree of cluster membership
  U <- as.data.frame(FKM_TPS$U, row.names=F)
  if (noise==TRUE) {
    sum <- apply(U,1,sum)
    U$noise <- 1-sum
    Umax <- c()
    for (q in 1:length(U$noise)) {
      Umax[q] <- which(U[q,]==max(U[q,]))
      if (Umax[q]==dim(U)[2]) {Umax[q] <- 0}
    }
  }
  if (noise==FALSE) {
    Umax <- c()
    for (q in 1:length(U[,1])) {
      Umax[q] <- which(U[q,]==max(U[q,]))
    }
  }

  subjects <- wide_data$Id2
  FKM_TPS_U <- cbind(subjects,U,Umax)

  # Column names
  column_names <- c("Id2")
  for (i in 1:k) {
    nextval <- paste("Clus",i,sep="")
    column_names <- append(column_names, nextval)
  }
  if (noise==TRUE) {column_names <- append(column_names, "Noise")}
  column_names <- append(column_names,"ClusModal")
  colnames(FKM_TPS_U) <- column_names

  FKM_TPS_U <- merge(TPSdata$IDmatch, FKM_TPS_U, by="Id2")
  log <- capture.output({
    FKM_indices <- fclust::Fclust.index(FKM_TPS);
  })
  # Values to return
  output <- list(FKM_TPS=FKM_TPS, U=U, Umax=Umax, FKM_TPS_U=FKM_TPS_U, FKM_indices=FKM_indices, wide_data=wide_data, TPSdata=TPSdata, center=center, noise=noise, k=k, m=m, addslopes=addslopes, IDmatch=TPSdata$IDmatch, call=match.call())
  class(output) <- "FKM.TPS"
  return(output)
}


#' Fuzzy c-medoids clustering on tensor-product spline coefficients
#'
#' `cluster.coefs()` is used to apply the fuzzy c-medoids clustering algorithm to
#'   the tensor-product spline coefficients that are output from the [TPSfit()] function.
#'
#' `cluster.coefs()` employs package \pkg{fclust} to apply th fuzzy c-medoids
#' clustering algorithm to the tensor-product spline coefficients that are
#' output from the [TPSfit()] function. The user has the option to include a noise
#' cluster.
#'
#' @param TPSdata Output from the [TPSfit()] function.
#' @param k An integer value indicating the number of clusters.
#' @param m A numeric value called the 'fuzziness parameter.' Must be greater
#'   than 1, with a default of 1.2.
#' @param noise A logical expression indicating whether a noise cluster for
#'   outliers is to be included.
#' @param seed An optional numeric value for the seed, as the clustering
#'   algorithm uses random start values.
#' @param ... Additional optional values employed by \pkg{fclust}, including `RS`,
#'   the number of random starts, and `delta`, the distance for outliers.
#'
#' @returns An object of class '`FKM.TPS`' containing the following components:
#'* `FKM_TPS ` An object output by the \pkg{fclust} algorithm.
#'* `U `  A dataframe containing the degrees of cluster membership for each subject.
#'* `Umax ` A vector containing the modal cluster assignment for each subject.
#'* `FKM_TPS_U ` A data frame containing subjects with subject ID, degrees of cluster
#'   membership, and modal class assignment.
#'* `FKM_indices ` A vector containing the values for six cluster validity indices.
#'   See [fclust::Fclust.index()]
#'* `wide_data ` A data frame with the values that were clustered on in the
#'   clustering algorithm plus subject IDs.
#'* `TPSdata ` The '`TPSfit`' object that was input into the function.
#'* `noise ` Logical value that was used to determine whether a noise cluster is
#'   included.
#'* `k ` The number of clusters.
#'* `m ` The fuzziness parameter.
#'* `IDmatch ` A data frame matching the original subject ID with new consecutive
#'   ID values.
#'
#' @export
#'
#' @seealso
#' * [cluster.fitted()] is an alternative method that clusters fitted spline
#'   values rather than spline coefficients.
#' * The pkg{fclust} R package: <https://cran.r-project.org/web/packages/fclust/index.html>
#' * Fuzzy k-medoids function: [fclust::FKM.med()]
#' * Fuzzy k-medoids function with noise cluster: [fclust::FKM.med.noise()]
#'
#' @examples
#' data(TS.sim)
#'
#' fitsplines <- TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
#'      ID="SubjectID", knots_time=c(0, 91, 182, 273, 365), n_fit_times=10)
#'
#' clusters2 <- cluster.coefs(fitsplines, k=3, m=1.3, seed=12345, RS=5, noise=TRUE)
#' summary(clusters2)
#'
cluster.coefs <- function(TPSdata, k, m=1.2, noise=TRUE, seed, ...) {
  GAMscoef <- TPSdata$GAMscoef

  # Convert data from long to wide
  wide_data <- as.data.frame(tidyr::pivot_wider(GAMscoef, names_from=.data$coef, names_prefix="coef",
                                         values_from=.data$GAMbasis))

  # Remove ID
  cluster_data <- wide_data[, c(-1, -2)]

  # Clustering using fclust package
  if (!missing(seed)) {set.seed(seed)}
  if (noise==TRUE) {
    FKM_TPS <- fclust::FKM.med.noise(cluster_data, k=k, m=m, ...)
  }
  if (noise==FALSE) {
    FKM_TPS <- fclust::FKM.med(cluster_data, k=k, m=m, ...)
  }
  if (noise!=TRUE & noise!=FALSE) {stop('noise should be TRUE/FALSE')}

  # Dataset with degree of cluster membership
  U <- as.data.frame(FKM_TPS$U, row.names=F)
  if (noise==TRUE) {
    sum <- apply(U,1,sum)
    U$noise <- 1-sum
    Umax <- c()
    for (q in 1:length(U$noise)) {
      Umax[q] <- which(U[q,]==max(U[q,]))
      if (Umax[q]==dim(U)[2]) {Umax[q] <- 0}
    }
  }
  if (noise==FALSE) {
    Umax <- c()
    for (q in 1:length(U[,1])) {
      Umax[q] <- which(U[q,]==max(U[q,]))
    }
  }

  subjects <- wide_data$Id2
  FKM_TPS_U <- cbind(subjects,U,Umax)

  # Column names
  column_names <- c("Id2")
  for (i in 1:k) {
    nextval <- paste("Clus",i,sep="")
    column_names <- append(column_names, nextval)
  }
  if (noise==TRUE) {column_names <- append(column_names, "Noise")}
  column_names <- append(column_names,"ClusModal")
  colnames(FKM_TPS_U) <- column_names

  FKM_TPS_U <- merge(TPSdata$IDmatch, FKM_TPS_U, by="Id2")
  log <- capture.output({
    FKM_indices <- fclust::Fclust.index(FKM_TPS);
  })
  # Values to return
  output <- list(FKM_TPS=FKM_TPS, U=U, Umax=Umax, FKM_TPS_U=FKM_TPS_U, FKM_indices=FKM_indices, wide_data=wide_data,
                 TPSdata=TPSdata, noise=noise, k=k, m=m, IDmatch=TPSdata$IDmatch, call=match.call())
  class(output) <- "FKM.TPS"
  return(output)
}


#' @export
#'
print.FKM.TPS <- function(x, ...) {
  if (class(x) != "FKM.TPS") {stop('Please supply an object of class "FKM.TPS"')}
  checktype <- substr(x$call[[1]], 9, 11)
  cat("Object of type 'FKM.TPS'\n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  if (x$noise==FALSE) {
    cat(x$TPSdata$nsubjects, "subjects clustered into", x$k,
        "clusters using m =", x$m, "\n", sep=" ") }
  if (x$noise==TRUE) {
    cat(x$TPSdata$nsubjects, "subjects clustered into", x$k,
        "clusters + noise cluster using m =", x$m, "\n", sep=" ") }
  if (checktype=="coe") {cat("Clusters based on tensor-product splines coefficients\n")}
  if (checktype=="fit") {if (x$addslopes==TRUE) {
    cat("Clusters based on fitted values at times", x$TPSdata$fit_times, "\n", sep=" ")
    cat("and slopes between points\n")
  }
    if (x$addslopes==FALSE) {
      cat("Clusters based on fitted values at times", x$TPSdata$fit_times, "\n", sep=" ")
    }}
  cat("\n")
  tableU <- table(x$Umax)
  if (names(tableU[1])=="0") {names(tableU)[1] <- "noise" }
  cat("Frequency of modal cluster assignments:\n")
  print(tableU)

  cat("\n")
  meandegree <- round(colMeans(x$U), digits=4)
  clusters <- sort(unique(x$FKM_TPS_U$ClusModal))
  dimclus <- length(clusters)
  for (i in clusters) {
    if (i==0) {meandegree[dimclus] <- colMeans(x$FKM_TPS_U[which(x$FKM_TPS_U$ClusModal==0),])[dimclus+1]}
    if (i!=0) {meandegree[i] <- colMeans(x$FKM_TPS_U[which(x$FKM_TPS_U$ClusModal==i),])[i+1]}
  }
  cat("Mean degree of cluster membership for modal cluster assignments:\n")
  print(meandegree)
  cat("\n")
}



#' @export
#'
summary.FKM.TPS <- function(object, ...) {
  if (class(object) != "FKM.TPS") {stop('Please supply an object of class "FKM.TPS"')}
  checktype <- substr(object$call[[1]], 9, 11)
  print(object$call)
  cat("\n")
  if (object$noise==FALSE) {
    cat(object$TPSdata$nsubjects, "subjects clustered into", object$k,
        "clusters using m =", object$m, "\n", sep=" ") }
  if (object$noise==TRUE) {
    cat(object$TPSdata$nsubjects, "subjects clustered into", object$k,
        "clusters + noise cluster using m =", object$m, "\n", sep=" ") }
  if (checktype=="coe") {cat("Clusters based on tensor-product splines coefficients\n")}
  if (checktype=="fit") {if (object$addslopes==TRUE) {
    cat("Clusters based on fitted values at times", object$TPSdata$fit_times, "\n", sep=" ")
    cat("and slopes between points\n")
  }
    if (object$addslopes==FALSE) {
      cat("Clusters based on fitted values at times", object$TPSdata$fit_times, "\n", sep=" ")
    }}
  cat("\n")
  cat("Cluster summary:\n")
  clusters <- sort(unique(object$FKM_TPS_U$ClusModal))
  dimclus <- length(clusters)
  clusmat <- matrix(NA, nrow=dimclus, ncol=4)

  for (i in clusters) {
    clusi <- object$FKM_TPS_U[which(object$FKM_TPS_U$ClusModal==i),]
    if (i==0) {
      col <- dimclus+2
      clusmat[dimclus,1] <- length(clusi$Id2)
      clusmat[dimclus,2] <- round(min(clusi[,col]), digits=3)
      clusmat[dimclus,3] <- round(max(clusi[,col]), digits=3)
      clusmat[dimclus,4] <- round(mean(clusi[,col]), digits=3)
    }
    if (i!=0) {
      col <- i+2
      clusmat[i,1] <- length(clusi$Id2)
      clusmat[i,2] <- round(min(clusi[,col]), digits=3)
      clusmat[i,3] <- round(max(clusi[,col]), digits=3)
      clusmat[i,4] <- round(mean(clusi[,col]), digits=3)
    }
  }
  colnames(clusmat) <- c("Cl.size", "Min.degree", "Max.degree", "Mean.degree")
  rownames(clusmat) <- colnames(object$FKM_TPS_U)[3:(dimclus+2)]
  print(clusmat)
  cat("\n")
  cat("Component 'FKM_TPS' contains the fuzzy clustering details from package 'fclust'.\n")
  cat("Components of 'FKM_TPS':\n")
  print(names(object$FKM_TPS))
  cat("\n")
  cat("Cluster validity indices:\n")
  print(object$FKM_indices)
  cat("\n")
  cat("Output dataset 'FKM_TPS_U' contains degree of cluster membership and modal
      cluster assignment for each object.\n")
  cat("Head of dataset 'FKM_TPS_U':\n")
  print(head(object$FKM_TPS_U))
  cat("\n")
}



#' Predict the cluster assignment for new subjects
#'
#' After clustering a set of trajectory data, predict the cluster assignment for
#' new subjects.
#'
#' @param object An object of class '`FKM.TPS`' of previously identified
#'   clusters.
#' @param newdata A new data frame containing subjects with the same variables
#'   previously used for clustering.
#' @param ... Additional arguments.
#'
#' @returns An object of class '`FKM.predicted`' containing the following components:
#' * `predicted_U ` A data frame containing the subjects with their degree of cluster
#'    membership and modal cluster assignment.
#' * `U ` The degrees of cluster membership for each subject.
#' * `Umax ` The modal class assignment for each subject.
#' * `wide_data ` A data frame containing the values that were clustered on.
#' * `IDmatch ` A data frame matching the original subject ID variable to a new
#'    consecutive subject ID.
#' * `TPSdata ` An object of class '`TPSfit`' containing the fitted splines for the
#'    new subjects.
#' * `FKM_TPS ` The inputted object of class '`FKM.TPS`'.
#' * `noise ` Logical expression indicating whether a noise cluster is included for
#'    outliers.
#' * `k ` The number of clusters.
#' * `m ` The fuzziness parameter.
#'
#' @export
#'
#' @examples
#' data(TS.sim)
#'
#' fitsplines <- TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
#'      ID="SubjectID", knots_time=c(0, 91, 182, 273, 365), n_fit_times=10)
#'
#' clusters1 <- cluster.fitted(fitsplines, k=3, m=1.3, seed=12345, RS=5, noise=TRUE)
#'
#' predicted_clusters <- predict(clusters1, TS.sim.new)
#' summary(predicted_clusters)
#'
predict.FKM.TPS <- function(object, newdata, ...) {
  if (class(object) != "FKM.TPS") {stop('Please supply an object of class "FKM.TPS"')}
  if (missing(newdata)) {stop('Please supply dataset for prediction "newdata".')}
  TPSdata <- object$TPSdata
  newsplines <- update(TPSdata, data=newdata)
  IDmatch <- newsplines$IDmatch
  checktype <- substr(object$call[[1]], 9, 11)

  if (checktype=="fit") {
    GAMsfitted <- newsplines$GAMsfitted
    center <- object$center

    # Get correct values to cluster on (centered or not)
    if (center==TRUE) {GAMsfitted$X <- GAMsfitted$centered_x}
    if (center==FALSE) {GAMsfitted$X <- GAMsfitted$x}
    if (center!= TRUE & center!=FALSE) {stop('center must be TRUE/FALSE')}

    # Convert data from long to wide
    #long_data <- subset(GAMsfitted, select=c(Id2, Variable, t, X))
    long_data <- GAMsfitted[,c(1,4,6,9)]
    long_data <- dplyr::arrange(long_data, .data$Id2, .data$Variable, .data$t)
    wide_data <- as.data.frame(tidyr::pivot_wider(long_data, names_from=c(.data$Variable, .data$t), names_prefix="X", values_from=.data$X))

    # Add slopes:
    if (object$addslopes==TRUE) {
      npoints <- length(TPSdata$fit_times)
      nvars <- length(TPSdata$vars)
      wide_data_m <- as.matrix(wide_data)
      wide_slopes_m <- matrix(0, nrow=length(wide_data$Id2), ncol=nvars*(npoints-1))
      for (j in 1:length(wide_data$Id2)) {    #j in wide_data$Id2
        subjmat <- matrix(wide_data_m[j,-1], nrow=npoints, ncol=nvars)
        slopemat <- matrix(0, nrow=npoints-1, ncol=nvars)
        for (i in 2:npoints-1) {
          for (q in 1:nvars) {
            slopemat[i,q] <- subjmat[i+1,q]-subjmat[i,q]
          }
        }
        wide_slopes_m[j,] <- unmatrixB(slopemat)
      }
      wide_slopes <- as.data.frame(wide_slopes_m)
      colnames(wide_slopes) <- names(unmatrixB(slopemat))
      wide_data <- cbind(wide_data, wide_slopes)
    }
  }

  if (checktype=="coe") {
    GAMscoef <- object$TPSdata$GAMscoef

    # Convert data from long to wide
    wide_data <- as.data.frame(tidyr::pivot_wider(GAMscoef, names_from=.data$coef, names_prefix="coef",
                                           values_from=.data$GAMbasis))
  }

  # Remove ID
  cluster_data <- wide_data[,-1]
  noise <- object$noise

  # Bring in medoids from previous FKM object
  cluster_data_m <- rbind(object$wide_data[object$FKM_TPS$medoid,-1], cluster_data)

  # Find distances and functions for calculating U for new data
  distance_matrix <- as.matrix(dist(unname(as.matrix(cluster_data_m))))
  distance_matrix_sq <- distance_matrix^2
  dist3 <- function(a) {(1/a)^(1/(m-1))}
  delta <- object$FKM_TPS$delta
  m <- object$m
  k <- object$k
  delta_term <- (1/(delta^2))^(1/(m-1))

  if (noise==TRUE) {
    #calculate uic and uic+1
    nums <- dist3(distance_matrix_sq[1:k,])
    denoms <- colSums(nums)+delta_term
    denoms <- do.call("rbind", replicate(k, denoms, simplify=FALSE))
    U <- nums/denoms
    Unoise <- delta_term/(colSums(nums)+delta_term)
    U <- t(U)
    Unoise <- as.matrix(Unoise, ncol=1)
    U[is.nan(U)] <- 1
  }

  if (noise==FALSE) {
    #calculate uic
    nums <- dist3(distance_matrix_sq[1:k,])
    denoms <- colSums(nums)
    denoms <- do.call("rbind", replicate(k, denoms, simplify=FALSE))
    U <- nums/denoms
    U <- t(U)
    U[is.nan(U)] <- 1
  }

  U <- U[-(1:k),]  #Removes medoids
  U <- as.data.frame(U, row.names = FALSE)

  # Get modal cluster
  if (noise==TRUE) {
    sum <- apply(U,1,sum)
    U$noise <- 1-sum
    Umax <- c()
    for (q in 1:length(U$noise)) {
      Umax[q] <- which(U[q,]==max(U[q,]))
      if (Umax[q]==dim(U)[2]) {Umax[q] <- 0}
    }
  }
  if (noise==FALSE) {
    Umax <- c()
    for (q in 1:length(U[,1])) {
      Umax[q] <- which(U[q,]==max(U[q,]))
    }
  }

  # New data frame with subject IDs and U
  subjects <- wide_data$Id2
  predicted_U <- cbind(subjects,U,Umax)

  # Column names
  column_names <- c("Id2")
  for (i in 1:k) {
    nextval <- paste("Clus",i,sep="")
    column_names <- append(column_names, nextval)
  }
  if (noise==TRUE) {column_names <- append(column_names, "Noise")}
  column_names <- append(column_names,"ClusModal")
  colnames(predicted_U) <- column_names

  predicted_U <- merge(IDmatch, predicted_U, by="Id2")

  output <- list(predicted_U=predicted_U, U=U, Umax=Umax, wide_data=wide_data, IDmatch=IDmatch,
                 TPSdata=newsplines, FKM_TPS=object, noise=noise, k=k, m=m, call=match.call())
  class(output) <- "FKM.predicted"
  return(output)
}



#' @export
#'
print.FKM.predicted <- function(x, ...) {
  if (class(x) != "FKM.predicted") {stop('Please supply an object of class "FKM.predicted"')}
  checktype <- substr(x$FKM_TPS$call[[1]], 9, 11)
  nsub <- x$TPSdata$nsubjects-length(x$TPSdata$error_subjects)
  cat("Object of type 'FKM.predicted'\n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Tensor-product splines fit for", nsub,
      "out of", x$TPSdata$nsubjects, "subjects\n", sep=" ")
  if(!is.null(x$TPSdata$error_subjects)) {
    (cat("Unable to fit splines for subjects: ", x$TPSdata$error_subjects, "\n", sep=" "))}
  if (x$noise==FALSE) {
    cat(nsub, "subjects clustered into", x$k,
        "clusters\n", sep=" ") }
  if (x$noise==TRUE) {
    cat(nsub, "subjects clustered into", x$k,
        "clusters + noise cluster\n", sep=" ") }
  if (checktype=="coe") {cat("Clusters based on tensor-product splines coefficients\n")}
  if (checktype=="fit") {if (x$FKM_TPS$addslopes==TRUE) {
    cat("Clusters based on fitted values at times", x$TPSdata$fit_times, "\n", sep=" ")
    cat("and slopes between points\n")
  }
    if (x$FKM_TPS$addslopes==FALSE) {
      cat("Clusters based on fitted values at times", x$TPSdata$fit_times, "\n", sep=" ")
    }}
  cat("\n")
  tableU <- table(x$Umax)
  if (names(tableU[1])=="0") {names(tableU)[1] <- "noise" }
  cat("Frequency of modal cluster assignments:\n")
  print(tableU)

  cat("\n")
  meandegree <- round(colMeans(x$U), digits=4)
  clusters <- sort(unique(x$predicted_U$ClusModal))
  dimclus <- length(clusters)
  for (i in clusters) {
    if (i==0) {meandegree[dimclus] <- colMeans(x$predicted_U[which(x$predicted_U$ClusModal==0),])[dimclus+2]}
    if (i!=0) {meandegree[i] <- colMeans(x$predicted_U[which(x$predicted_U$ClusModal==i),])[i+2]}
  }
  cat("Mean degree of cluster membership for modal cluster assignments:\n")
  print(meandegree)
  cat("\n")
}


#' @export
#'
summary.FKM.predicted <- function(object, ...) {
  if (class(object) != "FKM.predicted") {stop('Please supply an object of class "FKM.predicted"')}
  checktype <- substr(object$FKM_TPS$call[[1]], 9, 11)
  nsub <- object$TPSdata$nsubjects-length(object$TPSdata$error_subjects)
  cat("Object of type 'FKM.predicted'\n")
  print(object$call)
  cat("\n")
  cat("Tensor-product splines fit for", nsub,
      "out of", object$TPSdata$nsubjects, "subjects\n", sep=" ")
  if(!is.null(object$TPSdata$error_subjects)) {
    (cat("Unable to fit splines for subjects: ", object$TPSdata$error_subjects, "\n", sep=" "))}
  cat("\n")
  cat("Degree of membership calculated based on clusters from input 'FKM.TPS' object.\n")
  if (object$noise==FALSE) {
    cat(nsub, "subjects clustered into", object$k,
        "clusters\n", sep=" ") }
  if (object$noise==TRUE) {
    cat(nsub, "subjects clustered into", object$k,
        "clusters + noise cluster\n", sep=" ") }
  cat("\n")
  if (checktype=="coe") {cat("Clusters based on tensor-product splines coefficients\n")}
  if (checktype=="fit") {if (object$FKM_TPS$addslopes==TRUE) {
    cat("Clusters based on fitted values at times", object$TPSdata$fit_times, "\n", sep=" ")
    cat("and slopes between points\n")
  }
    if (object$FKM_TPS$addslopes==FALSE) {
      cat("Clusters based on fitted values at times", object$TPSdata$fit_times, "\n", sep=" ")
    }}
  cat("\n")
  cat("Cluster summary for new data:\n")
  clusters <- sort(unique(object$predicted_U$ClusModal))
  dimclus <- length(clusters)
  clusmat <- matrix(NA, nrow=dimclus, ncol=4)

  for (i in clusters) {
    clusi <- object$predicted_U[which(object$predicted_U$ClusModal==i),]
    if (i==0) {
      col <- dimclus+2
      clusmat[dimclus,1] <- length(clusi$Id2)
      clusmat[dimclus,2] <- round(min(clusi[,col]), digits=3)
      clusmat[dimclus,3] <- round(max(clusi[,col]), digits=3)
      clusmat[dimclus,4] <- round(mean(clusi[,col]), digits=3)
    }
    if (i!=0) {
      col <- i+2
      clusmat[i,1] <- length(clusi$Id2)
      clusmat[i,2] <- round(min(clusi[,col]), digits=3)
      clusmat[i,3] <- round(max(clusi[,col]), digits=3)
      clusmat[i,4] <- round(mean(clusi[,col]), digits=3)
    }
  }
  colnames(clusmat) <- c("Cl.size", "Min.degree", "Max.degree", "Mean.degree")
  rownames(clusmat) <- colnames(object$predicted_U)[3:(dimclus+2)]
  print(clusmat)
  cat("\n")
  cat("Output dataset 'predicted_U' contains degree of cluster membership and modal
      cluster assignment for each object.\n")
  cat("Head of dataset 'predicted_U':\n")
  print(head(object$predicted_U))
  cat("\n")
}


#' Compare indices for fuzzy clusterings
#'
#' `comp.FKM.indices()` compares the clustering indices of two more more
#' clusterings that are output from [cluster.coefs()] or [cluster.fitted()].
#'
#' `comp.FKM.indices()` compares the clustering indices of two more more
#' clusterings that are output from [cluster.coefs()] or [cluster.fitted()].
#' Available indices are those found in the \pkg{fclust} package: `PC`
#' (partition coefficient), `PE` (partition entropy), `MPC` (modified partition
#' coefficient), `SIL` (silhouette), `SIL.F` (fuzzy silhouette), and `XB`
#' (Xie-Beni). The default `ALL` gives all indices. See [fclust::Fclust.index()]
#'
#' @param clusterings A list object of the clusterings to compare. Each should
#'   be of class `'FKM.TPS'`.
#' @param index Desired cluster validity index or indices. Default is `"ALL"`.
#'   Options are: `"PC"`, `"PE"`, `"MPC"`, `"SIL"`, `"SIL.F"`, `"XB"`, or
#'   `"ALL"`. Multiple options can be input as a character vector. See [fclust::Fclust.index()]
#' @param clusternames Optional character vector with names of the clusters being compared.
#'
#' @return A data frame with the desired validity indices for each of the clusterings.
#'
#' @export
#'
#' @seealso Values for cluster validity indices are calculated using the
#'   \pkg{fclust} package. See [fclust::Fclust.index()]
#'
#' @examples
#' data(TS.sim)
#'
#' fitsplines <- TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
#'      ID="SubjectID", knots_time=c(0, 91, 182, 273, 365), n_fit_times=10)
#'
#' ccoefs_2 <- cluster.coefs(fitsplines, k=2, seed=1234, RS=10)
#' ccoefs_3 <- cluster.coefs(fitsplines, k=2, seed=1234, RS=10)
#' ccoefs_4 <- cluster.coefs(fitsplines, k=2, seed=1234, RS=10)
#' ccoefs_5 <- cluster.coefs(fitsplines, k=2, seed=1234, RS=10)
#'
#' # Compare clusters using all indices and custom names
#' comp.FKM.indices(list(ccoefs_2, ccoefs_3, ccoefs_4, ccoefs_5),
#' clusternames=c("k=2", "k=3", "k=4", "k=5"))
#'
#' # Compare clusterings using a subset of the indices
#' comp.FKM.indices(list(ccoefs_2, ccoefs_3, ccoefs_4, ccoefs_5),
#' clusternames=c("k=2", "k=3", "k=4", "k=5"), index=c("SIL.F", "XB"))
#'
comp.FKM.indices <- function(clusterings, index="ALL", clusternames=NA) {
  if (class(clusterings[[1]]) != "FKM.TPS") {stop("Please supply two or more objects of class 'FKM.TPS' as a list")}
  if (!missing(clusternames) & length(clusternames) != length(clusterings)) {stop("Length of clusternames must be equal to the number of clusterings")}
  if (!missing(clusternames) & !is.vector(clusternames)) {stop("clusternames must be a character vector")}

  ncompare <- length(clusterings)
  comp_matrix <- matrix(NA, nrow=ncompare, ncol=6)
  colnames(comp_matrix) <- names(clusterings[[1]]$FKM_indices)
  for (i in 1:ncompare) {
    comp_matrix[i,] <- clusterings[[i]]$FKM_indices
  }

  #Get desired indices
  if (index[1] != "ALL" & length(index)==1) {
    comp_matrix2 <- as.data.frame(comp_matrix[,index])
    colnames(comp_matrix2) <- index
  }
  if (index[1]=="ALL") {
    comp_matrix2 <- comp_matrix
  }
  if (length(index)>1) {comp_matrix2 <- comp_matrix[,index]}
  if (missing(clusternames)) {row.names(comp_matrix2) <- names(clusterings)}
  if (!missing(clusternames)) {row.names(comp_matrix2) <- clusternames}
  return(comp_matrix2)
}

