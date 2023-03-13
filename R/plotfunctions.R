#' Plot an object of class '`FKM.TPS`'
#'
#' Produces a number of plots for an object of class `FKM.TPS` using [ggplot2::ggplot()]. Possible plots include mean trajectory by cluster, spaghetti plot of raw data by cluster overlaid or in grid format, and smoothed spaghetti plot by cluster overlaid or in grid format.
#'
#' @param x A object of class '`FKM.TPS`' output from [cluster.fitted()] or [cluster.coefs()]
#' @param center Logical expression indicating whether to center trajectories on individual means.
#' @param xmin,xmax Optional minimum and maximum values to show on x-axis.
#' @param ntime Optional number of times to calculate fitted values for smoothed plots.
#' @param lab_x,lab_y Optional labels for x- and y-axis.
#' @param bw Logical expression for black and white graphic.
#' @param title Optional title.
#' @param title_size Optional title size.
#' @param axis_label_size Optional size of axis labels.
#' @param axis_title_size Optional size for axis titles.
#' @param legend_label_size Optional size for legend.
#' @param strip_label_size Optional size for strip labels on graphics.
#' @param type Type of plot to produce. Options are `"mean"`, `"smooth"`, `"smooth_grid"`, `"raw"`, `"raw_grid"`, and `"all"`.
#' @param ... Additional arguments
#'
#' @returns A plot of the data by cluster.
#' @export
#'
#' @examples
#' data(TS.sim)
#'
#' fitsplines <- TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
#'      ID="SubjectID", knots_time=c(0, 91, 182, 273, 365), n_fit_times=10)
#'
#' clusters1 <- cluster.fitted(fitsplines, k=3, m=1.3, seed=12345, RS=5, noise=TRUE)
#' plot(clusters1, type="raw_grid", strip_label_size=10, axis_label_size=10)
#' plot(clusters1, type="mean", legend_label_size=15, lab_y="Outcome")

plot.FKM.TPS <- function(x, center=TRUE, xmin, xmax, ntime=100, lab_x, lab_y, bw=FALSE, title, title_size=15, axis_label_size=15, axis_title_size=15, legend_label_size=15, strip_label_size=15, type="mean", ...) {
  if (class(x) != "FKM.TPS") {stop('Please supply an object of class "FKM.TPS"')}
  checktype <- substr(x$call[[1]], 9, 11)
  noise <- x$noise

  GAMs <- x$TPSdata$GAMs
  if (missing(xmin)) {xmin <- min(x$TPSdata$fit_times)}
  if (missing(xmax)) {xmax <- max(x$TPSdata$fit_times)}

  #Get smooth trajectories
  timegraph <- seq(from=xmin, to=xmax, length.out=ntime)
  vars <- x$TPSdata$vars
  numvars <- length(vars)
  timegraph2 <- rep(timegraph, times=numvars)
  varnums <- 1:numvars
  vargraph <- rep(varnums, each=ntime)
  graphdat <- data.frame(timegraph2, vargraph)
  names(graphdat) <- c("Time", "Variable")

  graphvals <- list()
  nsubjects <- x$TPSdata$nsubjects
  for (i in 1:nsubjects) {
    Id2 <- rep(i, times=ntime*numvars)
    newgraph <- cbind(Id2, graphdat)
    if (!(i %in% x$TPSdata$error_subjects)) {newgraph$ValueS <- mgcv::predict.gam(GAMs[[i]], newgraph)}
    graphvals[[i]] <- newgraph
  }
  graphdata <- dplyr::bind_rows(graphvals)

  #Get centered value if needed
  if (center==TRUE) {
    indiv_means <- x$TPSdata$indiv_means
    for (c in 1:numvars) {
      indiv_means[[c]]$Variable <- c
    }
    indiv_means <- dplyr::bind_rows(indiv_means)
    graphdata <- merge(graphdata, indiv_means, by=c("Id2", "Variable"))
    graphdata$centered_ValueS <- graphdata$ValueS-graphdata$mean_x
    graphdata$X <- graphdata$centered_ValueS
  }
  if (center==FALSE) {graphdata$X <- graphdata$ValueS}

  graphdata <- merge(graphdata, x$FKM_TPS_U, by="Id2")
  graphdata$Clus <- as.factor(graphdata$ClusModal)
  graphdata$Var <- factor(graphdata$Variable, labels=vars)
  graphdata$Cluster <- "Noise"
  graphdata$Cluster[which(graphdata$ClusModal!=0)] <- as.character(graphdata$ClusModal[which(graphdata$ClusModal!=0)])

  graphdata2 <- graphdata[which(graphdata$ClusModal!=0),]

  rawdata <- merge(x$TPSdata$data_long, x$FKM_TPS_U, by="Id2")
  rawdata$Clus <- as.factor(rawdata$ClusModal)
  rawdata$Var <- factor(rawdata$Variable, labels=vars)
  rawdata$Cluster <- "Noise"
  rawdata$Cluster[which(rawdata$ClusModal!=0)] <- as.character(rawdata$ClusModal[which(rawdata$ClusModal!=0)])

  rawdata2 <- rawdata[which(rawdata$ClusModal!=0),]

  if (missing(lab_y)) {lab_y <- "Value"}
  if (missing(lab_x)) {lab_x <- "Time"}
  if (missing(title)) {title <- ""}

  #Get labels with number per cluster
  if (noise==TRUE) {tableU <- as.numeric(table(x$Umax)[-1])}
  if (noise==FALSE) {tableU <- as.numeric(table(x$Umax))}
  labels <- rep("Cluster ", times=length(tableU))
  clusnums <- 1:length(tableU)
  neq <- rep("n", times=length(tableU))
  for (i in 1:length(tableU)) {neq[i] <- paste(" (n=",tableU[i], ")  ", sep="")}
  labels <- paste(labels, clusnums, neq, sep="")

  #Get labels with number per cluster including noise
  if (noise==TRUE) {
    tableU2 <- as.numeric(table(x$Umax))
    labels2 <- rep("Cluster ", times=(length(tableU2)-1))
    clusnums2 <- 1:(length(tableU2)-1)
    #clusnums2 <- c("noise",clusnums2)
    neq2 <- rep("n", times=length(tableU2))
    for (i in 1:length(tableU2)) {neq2[i] <- paste(" (n=",tableU2[i], ")  ", sep="")}
    labels2 <- paste(labels2, clusnums2, sep="")
    labels2 <- c("Noise cluster", labels2)
    labels2 <- paste(labels2, neq2, sep="")
  }
  if (noise==FALSE) {labels2 <- labels}

  if (type=="raw" | type=="all") {
    if (bw==FALSE) {
      show(ggplot2::ggplot(rawdata, ggplot2::aes(x=.data$Time, y=.data$x, group=.data$Id2, colour=.data$Clus)) +
             ggplot2::geom_line() +
             ggplot2::facet_wrap(~Var) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::scale_color_discrete(name="", labels=labels2) +
             ggplot2::theme_bw() +
             ggplot2::theme(legend.position = "bottom", axis.text=ggplot2::element_text(size=axis_label_size),
                            legend.text = ggplot2::element_text(size=legend_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            legend.key.width = ggplot2::unit(1, 'cm'), strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
    if (bw==TRUE) {
      show(ggplot2::ggplot(rawdata, ggplot2::aes(x=.data$Time, y=.data$x, group=.data$Id2, linetype= .data$Clus)) +
             ggplot2::geom_line() +
             ggplot2::facet_wrap(~Var) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::scale_linetype_discrete(name="", labels=labels2) +
             ggplot2::theme_bw() +
             ggplot2::theme(legend.position = "bottom", axis.text=ggplot2::element_text(size=axis_label_size),
                            legend.text = ggplot2::element_text(size=legend_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            legend.key.width = ggplot2::unit(1, 'cm'), strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
  }

  if (type=="raw_grid" | type=="all") {
    if (bw==FALSE) {
      show(ggplot2::ggplot(rawdata, ggplot2::aes(x=.data$Time, y= .data$x, group= .data$Id2)) +
             ggplot2::geom_line(color="cornflowerblue") +
             ggplot2::facet_grid(Cluster~Var, labeller=ggplot2::labeller(Cluster=ggplot2::label_both)) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.text=ggplot2::element_text(size=axis_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
    if (bw==TRUE) {
      show(ggplot2::ggplot(rawdata, ggplot2::aes(x=.data$Time, y=.data$x, group=.data$Id2)) +
             ggplot2::geom_line(color="black") +
             ggplot2::facet_grid(Cluster~Var, labeller=ggplot2::labeller(Cluster=ggplot2::label_both)) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.text=ggplot2::element_text(size=axis_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
  }

  if (type=="smooth" | type=="all") {
    if (bw==FALSE) {
      show(ggplot2::ggplot(data=graphdata, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Id2, colour=.data$Clus)) +
             ggplot2::geom_line() + ggplot2::facet_wrap(~Var) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::scale_color_discrete(name="", labels=labels2) +
             ggplot2::theme(legend.position = "bottom", axis.text=ggplot2::element_text(size=axis_label_size),
                            legend.text = ggplot2::element_text(size=legend_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            legend.key.width = ggplot2::unit(1, 'cm'), strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
    if (bw==TRUE) {
      show(ggplot2::ggplot(data=graphdata, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Id2, colour=.data$Clus)) +
             ggplot2::geom_line() + ggplot2::facet_wrap(~Var) +
             ggplot2::theme_bw() +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::scale_colour_grey(start=0, end=0.6, labels=labels2, name="") +
             ggplot2::theme(legend.position = "bottom", axis.text=ggplot2::element_text(size=axis_label_size),
                            legend.text = ggplot2::element_text(size=legend_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            legend.key.width = ggplot2::unit(1, 'cm'), strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
  }

  if (type=="smooth_grid" | type=="all") {
    if (bw==FALSE) {
      show(ggplot2::ggplot(data=graphdata, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Id2)) +
             ggplot2::geom_line(color="cornflowerblue") +
             ggplot2::facet_grid(Cluster~Var, labeller=ggplot2::labeller(Cluster=ggplot2::label_both)) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.text=ggplot2::element_text(size=axis_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
    if (bw==TRUE) {
      show(ggplot2::ggplot(data=graphdata, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Id2)) +
             ggplot2::geom_line(color="black") +
             ggplot2::facet_grid(Cluster~Var, labeller=ggplot2::labeller(Cluster=ggplot2::label_both)) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.text=ggplot2::element_text(size=axis_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
  }

  if (type=="mean" | type=="all") {
    if (bw==FALSE) {
      show(ggplot2::ggplot(data=graphdata2, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Clus, colour=.data$Clus)) +
             ggplot2::facet_wrap(~Var) +
             ggplot2::geom_smooth(stat='smooth') +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::scale_color_discrete(name="", labels=labels) +
             ggplot2::theme_bw() +
             ggplot2::theme(legend.position = "bottom", axis.text=ggplot2::element_text(size=axis_label_size),
                            legend.text = ggplot2::element_text(size=legend_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            legend.key.width = ggplot2::unit(1.5, 'cm'), strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
    if (bw==TRUE) {
      show(ggplot2::ggplot(data=graphdata2, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Clus, linetype=.data$Clus)) +
             ggplot2::facet_wrap(~Var) +
             ggplot2::geom_smooth(stat='smooth', color='black') +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::scale_linetype_discrete(name="", labels=labels) +
             ggplot2::theme_bw() +
             ggplot2::theme(legend.position = "bottom", axis.text=ggplot2::element_text(size=axis_label_size),
                            legend.text = ggplot2::element_text(size=legend_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            legend.key.width = ggplot2::unit(1.5, 'cm'), strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
  }
}


#' Plot an object of class '`TPSfit`'.
#'
#' Plots an object of class '`TPSfit`', with option to plot either raw or smoothed data.
#'
#' @param x Object of class '`TPSfit`'.
#' @param center Logical expression indicating whether to center data on individual trajectories.
#' @param xmin,xmax Optional minimum and maximum values to show on x-axis.
#' @param ntime Optional number of times to calculate fitted values for smoothed plots.
#' @param lab_x,lab_y Optional labels for x- and y-axis.
#' @param axis_label_size Optional size of axis labels.
#' @param axis_title_size Optional size for axis titles.
#' @param legend_label_size Size of legend label, if applicable
#' @param strip_label_size Optional size for strip labels on graphics.
#' @param title Optional title.
#' @param title_size Optional title size.
#' @param bw Logical expression for black and white graphic.
#' @param type Type of plot to produce. Options are `"smooth"`, `"raw"`, and `"all"`.
#' @param ... Additional arguments
#'
#' @returns A plot of the raw or smoothed data.
#' @export
#'
#' @examples
#' data(TS.sim)
#'
#' fitsplines <- TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
#'      ID="SubjectID", knots_time=c(0, 91, 182, 273, 365), n_fit_times=10)
#' plot(fitsplines, type="smooth")
#'
plot.TPSfit <- function(x, center=TRUE, xmin, xmax, ntime=100, lab_x, lab_y, bw=FALSE, title, title_size=15, axis_label_size=15, axis_title_size=15, legend_label_size=15, strip_label_size=15, type="raw", ...) {
  if (class(x) != "TPSfit") {stop('Please supply an object of class "TPSfit"')}

  GAMs <- x$GAMs
  if (missing(xmin)) {xmin <- min(x$fit_times)}
  if (missing(xmax)) {xmax <- max(x$fit_times)}

  #Get smooth trajectories
  timegraph <- seq(from=xmin, to=xmax, length.out=ntime)
  vars <- x$vars
  numvars <- length(vars)
  timegraph2 <- rep(timegraph, times=numvars)
  varnums <- 1:numvars
  vargraph <- rep(varnums, each=ntime)
  graphdat <- data.frame(timegraph2, vargraph)
  names(graphdat) <- c("Time", "Variable")

  graphvals <- list()
  nsubjects <- x$nsubjects
  for (i in 1:nsubjects) {
    Id2 <- rep(i, times=ntime*numvars)
    newgraph <- cbind(Id2, graphdat)
    if (!(i %in% x$error_subjects)) {newgraph$ValueS <- mgcv::predict.gam(GAMs[[i]], newgraph)}
    graphvals[[i]] <- newgraph
  }
  graphdata <- dplyr::bind_rows(graphvals)

  #Get centered value if needed
  if (center==TRUE) {
    indiv_means <- x$indiv_means
    for (c in 1:numvars) {
      indiv_means[[c]]$Variable <- c
    }
    indiv_means <- dplyr::bind_rows(indiv_means)
    graphdata <- merge(graphdata, indiv_means, by=c("Id2", "Variable"))
    graphdata$centered_ValueS <- graphdata$ValueS-graphdata$mean_x
    graphdata$X <- graphdata$centered_ValueS
  }
  if (center==FALSE) {graphdata$X <- graphdata$ValueS}

  graphdata$Var <- factor(graphdata$Variable, labels=vars)

  rawdata <- x$data_long
  rawdata$Var <- factor(rawdata$Variable, labels=vars)

  if (missing(lab_y)) {lab_y <- "Value"}
  if (missing(lab_x)) {lab_x <- "Time"}
  if (missing(title)) {title <- ""}

  if (type=="raw" | type=="all") {
    if (bw==FALSE) {
      show(ggplot2::ggplot(rawdata, ggplot2::aes(x=.data$Time, y=.data$x, group=.data$Id2)) +
             ggplot2::geom_line(color="cornflowerblue") +
             ggplot2::facet_wrap(~Var) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme( axis.text=ggplot2::element_text(size=axis_label_size),
                    axis.title=ggplot2::element_text(size=axis_title_size),
                    strip.text=ggplot2::element_text(size=strip_label_size),
                    plot.title=ggplot2::element_text(size=title_size)) )
    }
    if (bw==TRUE) {
      show(ggplot2::ggplot(rawdata, ggplot2::aes(x=.data$Time, y=.data$x, group=.data$Id2)) +
             ggplot2::geom_line(color="black") +
             ggplot2::facet_wrap(~Var) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme( axis.text=ggplot2::element_text(size=axis_label_size),
                    axis.title=ggplot2::element_text(size=axis_title_size),
                    strip.text=ggplot2::element_text(size=strip_label_size),
                    plot.title=ggplot2::element_text(size=title_size)) )
    }
  }

  if (type=="smooth" | type=="all") {
    if (bw==FALSE) {
      show(ggplot2::ggplot(data=graphdata, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Id2)) +
             ggplot2::geom_line(color="cornflowerblue") +
             ggplot2::facet_wrap(~Var) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.text=ggplot2::element_text(size=axis_label_size),
                   axis.title=ggplot2::element_text(size=axis_title_size),
                   strip.text=ggplot2::element_text(size=strip_label_size),
                   plot.title=ggplot2::element_text(size=title_size)) )
    }
    if (bw==TRUE) {
      show(ggplot2::ggplot(data=graphdata, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Id2)) +
             ggplot2::geom_line(color="black") +
             ggplot2::facet_wrap(~Var) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.text=ggplot2::element_text(size=axis_label_size),
                   axis.title=ggplot2::element_text(size=axis_title_size),
                   strip.text=ggplot2::element_text(size=strip_label_size),
                   plot.title=ggplot2::element_text(size=title_size)) )
    }
  }
}


#' Plot an object of class '`FKM.predicted`'
#'
#' @param x An object of class '`FKM.predicted`'
#' @param center Logical expression indicating whether trajectories are centered
#'   on individual means.
#' @param xmin,xmax Optional minimum and maximum values to show on x-axis.
#' @param ntime Optional number of times to calculate fitted values for smoothed
#'   plots.
#' @param lab_x,lab_y Optional labels for x- and y-axis.
#' @param bw Logical expression for black and white graphic.
#' @param title Optional title.
#' @param title_size Optional title size.
#' @param axis_label_size Optional size of axis labels.
#' @param axis_title_size Optional size for axis titles.
#' @param legend_label_size Optional size for legend.
#' @param strip_label_size Optional size for strip labels on graphics.
#' @param type Type of plot to produce. Options are `"raw"`, `"raw_grid"`, `"smooth"`,
#'   and `"smooth_grid"`.
#' @param ... Additional arguments
#'
#' @returns A plot of raw or smoothed trajectories for new data by predicted cluster assignment.
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
#' plot(predicted_clusters, label_size=10, type="raw_grid")
plot.FKM.predicted <- function(x, center=TRUE, xmin, xmax, ntime=100, lab_x,
                               lab_y, bw=TRUE, title, title_size=15, axis_label_size=15, axis_title_size=15, legend_label_size=15, strip_label_size=15, type="raw_grid", ...) {
  if (class(x) != "FKM.predicted") {stop('Please supply an object of class "FKM.predicted"')}
  noise <- x$noise

  GAMs <- x$TPSdata$GAMs
  if (missing(xmin)) {xmin <- min(x$TPSdata$fit_times)}
  if (missing(xmax)) {xmax <- max(x$TPSdata$fit_times)}

  #Get smooth trajectories
  timegraph <- seq(from=xmin, to=xmax, length.out=ntime)
  vars <- x$TPSdata$vars
  numvars <- length(vars)
  timegraph2 <- rep(timegraph, times=numvars)
  varnums <- 1:numvars
  vargraph <- rep(varnums, each=ntime)
  graphdat <- data.frame(timegraph2, vargraph)
  names(graphdat) <- c("Time", "Variable")

  graphvals <- list()
  nsubjects <- x$TPSdata$nsubjects - length(x$TPSdata$error_subjects)
  for (i in x$wide_data$Id2) {
    Id2 <- rep(i, times=ntime*numvars)
    newgraph <- cbind(Id2, graphdat)
    if (!(i %in% x$FKM_TPS$TPSdata$error_subjects)) {newgraph$ValueS <- mgcv::predict.gam(GAMs[[i]], newgraph)}
    graphvals[[i]] <- newgraph
  }
  graphdata <- dplyr::bind_rows(graphvals)

  #Get centered value if needed
  if (center==TRUE) {
    indiv_means <- x$TPSdata$indiv_means
    for (c in 1:numvars) {
      indiv_means[[c]]$Variable <- c
    }
    indiv_means <- dplyr::bind_rows(indiv_means)
    graphdata <- merge(graphdata, indiv_means, by=c("Id2", "Variable"))
    graphdata$centered_ValueS <- graphdata$ValueS-graphdata$mean_x
    graphdata$X <- graphdata$centered_ValueS
  }
  if (center==FALSE) {graphdata$X <- graphdata$ValueS}

  graphdata <- merge(graphdata, x$predicted_U, by="Id2")
  graphdata$Clus <- as.factor(graphdata$ClusModal)
  graphdata$Var <- factor(graphdata$Variable, labels=vars)
  graphdata$Cluster <- "Noise"
  graphdata$Cluster[which(graphdata$ClusModal!=0)] <- as.character(graphdata$ClusModal[which(graphdata$ClusModal!=0)])

  graphdata2 <- graphdata[which(graphdata$ClusModal!=0),]

  rawdata <- merge(x$TPSdata$data_long, x$predicted_U, by="Id2")
  rawdata$Clus <- as.factor(rawdata$ClusModal)
  rawdata$Var <- factor(rawdata$Variable, labels=vars)
  rawdata$Cluster <- "Noise"
  rawdata$Cluster[which(rawdata$ClusModal!=0)] <- as.character(rawdata$ClusModal[which(rawdata$ClusModal!=0)])

  rawdata2 <- rawdata[which(rawdata$ClusModal!=0),]

  if (missing(lab_y)) {lab_y <- "Value"}
  if (missing(lab_x)) {lab_x <- "Time"}
  if (missing(title)) {title <- ""}

  #Get labels with number per cluster
  if (noise==TRUE) {tableU <- as.numeric(table(x$Umax)[-1])}
  if (noise==FALSE) {tableU <- as.numeric(table(x$Umax))}
  labels <- rep("Cluster ", times=length(tableU))
  clusnums <- 1:length(tableU)
  neq <- rep("n", times=length(tableU))
  for (i in 1:length(tableU)) {neq[i] <- paste(" (n=",tableU[i], ")  ", sep="")}
  labels <- paste(labels, clusnums, neq, sep="")

  #Get labels with number per cluster including noise
  if (noise==TRUE) {
    tableU2 <- as.numeric(table(x$Umax))
    labels2 <- rep("Cluster ", times=(length(tableU2)-1))
    clusnums2 <- 1:(length(tableU2)-1)
    #clusnums2 <- c("noise",clusnums2)
    neq2 <- rep("n", times=length(tableU2))
    for (i in 1:length(tableU2)) {neq2[i] <- paste(" (n=",tableU2[i], ")  ", sep="")}
    labels2 <- paste(labels2, clusnums2, sep="")
    labels2 <- c("Noise cluster", labels2)
    labels2 <- paste(labels2, neq2, sep="")
  }
  if (noise==FALSE) {labels2 <- labels}

  if (type=="raw" | type=="all") {
    if (bw==FALSE) {
      show(ggplot2::ggplot(rawdata, ggplot2::aes(x=.data$Time, y=.data$x, group=.data$Id2, colour=.data$Clus)) +
             ggplot2::geom_line() +
             ggplot2::facet_wrap(~Var) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::scale_color_discrete(name="", labels=labels2) +
             ggplot2::theme_bw() +
             ggplot2::theme(legend.position = "bottom", axis.text=ggplot2::element_text(size=axis_label_size),
                            legend.text = ggplot2::element_text(size=legend_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            legend.key.width = ggplot2::unit(1, 'cm'), strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
    if (bw==TRUE) {
      show(ggplot2::ggplot(rawdata, ggplot2::aes(x=.data$Time, y=.data$x, group=.data$Id2, linetype=.data$Clus)) +
             ggplot2::geom_line() +
             ggplot2::facet_wrap(~Var) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::scale_linetype_discrete(name="", labels=labels2) +
             ggplot2::theme_bw() +
             ggplot2::theme(legend.position = "bottom", axis.text=ggplot2::element_text(size=axis_label_size),
                            legend.text = ggplot2::element_text(size=legend_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            legend.key.width = ggplot2::unit(1, 'cm'), strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
  }

  if (type=="raw_grid" | type=="all") {
    if (bw==FALSE) {
      show(ggplot2::ggplot(rawdata, ggplot2::aes(x=.data$Time, y=.data$x, group=.data$Id2)) +
             ggplot2::geom_line(color="cornflowerblue") +
             ggplot2::facet_grid(Cluster~Var, labeller=ggplot2::labeller(Cluster=ggplot2::label_both)) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.text=ggplot2::element_text(size=axis_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
    if (bw==TRUE) {
      show(ggplot2::ggplot(rawdata, ggplot2::aes(x=.data$Time, y=.data$x, group=.data$Id2)) +
             ggplot2::geom_line(color="black") +
             ggplot2::facet_grid(Cluster~Var, labeller=ggplot2::labeller(Cluster=ggplot2::label_both)) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.text=ggplot2::element_text(size=axis_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
  }

  if (type=="smooth" | type=="all") {
    if (bw==FALSE) {
      show(ggplot2::ggplot(data=graphdata, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Id2, colour=.data$Clus)) +
             ggplot2::geom_line() + ggplot2::facet_wrap(~Var) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::scale_color_discrete(name="", labels=labels2) +
             ggplot2::theme(legend.position = "bottom", axis.text=ggplot2::element_text(size=axis_label_size),
                            legend.text = ggplot2::element_text(size=legend_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            legend.key.width = ggplot2::unit(1, 'cm'), strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
    if (bw==TRUE) {
      show(ggplot2::ggplot(data=graphdata, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Id2, colour=.data$Clus)) +
             ggplot2::geom_line() + ggplot2::facet_wrap(~Var) +
             ggplot2::theme_bw() +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::scale_colour_grey(start=0, end=0.6, labels=labels2, name="") +
             ggplot2::theme(legend.position = "bottom", axis.text=ggplot2::element_text(size=axis_label_size),
                            legend.text = ggplot2::element_text(size=legend_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            legend.key.width = ggplot2::unit(1, 'cm'), strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
  }

  if (type=="smooth_grid" | type=="all") {
    if (bw==FALSE) {
      show(ggplot2::ggplot(data=graphdata, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Id2)) +
             ggplot2::geom_line(color="cornflowerblue") +
             ggplot2::facet_grid(Cluster~Var, labeller=ggplot2::labeller(Cluster=ggplot2::label_both)) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.text=ggplot2::element_text(size=axis_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
    if (bw==TRUE) {
      show(ggplot2::ggplot(data=graphdata, ggplot2::aes(x=.data$Time, y=.data$X, group=.data$Id2)) +
             ggplot2::geom_line(color="black") +
             ggplot2::facet_grid(Cluster~Var, labeller=ggplot2::labeller(Cluster=ggplot2::label_both)) +
             ggplot2::labs(y=lab_y, x=lab_x, title=title) +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.text=ggplot2::element_text(size=axis_label_size), axis.title=ggplot2::element_text(size=axis_title_size),
                            strip.text=ggplot2::element_text(size=strip_label_size),
                            plot.title=ggplot2::element_text(size=title_size)))
    }
  }

}


#' Plot a GLM that uses clusters as predictors
#'
#' Employs the `plot.glm` function to plot the generalized linear model with clusters and additional covariates as predictors.
#'
#' @param x The output object of class '`FKM.glm`'.
#' @param ... Additional arguments.
#'
#' @return A `plot.glm` plot.
#' @export
#'
plot.FKM.glm <- function(x, ...) {
  if (class(x) != "FKM.glm") {stop('Please supply an object of class "FKM.glm"')}
  plot(x$model_full)
}
