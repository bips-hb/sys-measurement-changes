## Implementation of the statistical methods


# Help functions ---------------------------------------------------------------

compute_estimators <- function(y, syschange, cpts = NULL) {
  # Compute estimators to quantify the magnitude of the systematic change.
  # y         : observed measurements
  # syschange : systematic change values (true or estimated)
  # cpts      : change-point information; either NULL, a vector of indices,
  #             or a binary indicator vector of length nobs

  nobs <- length(y)

  ## compute measures of dispersion
  range <- max(syschange) - min(syschange)
  var <- var(syschange)
  median <- median(syschange)
  madm <- mean(abs(syschange - median(syschange)))

  ## extract change points
  if (is.null(cpts)) {
    # some methods are not able to estimate the change point
    cpts <- NA
    n_cpts <- NA
  } else {
    if (length(cpts) == nobs) {
      # binary indicator vector (1 = change point)
      cpts <- which(cpts == 1)
    }
    n_cpts <- length(cpts)
    # to return the output in a dataframe
    paste(cpts, collapse = ", ")
  }

  return(data.frame(
    range = range,
    var = var,
    madm = madm,
    n_cpts = n_cpts,
    cpts = paste(cpts, collapse = ", ")
  ))
}

find_change_points <- function(y_est, max_clusters = 4, lag = length(y_est) * 0.05) {
  # Detect change points in an estimated signal using k-means clustering.
  # y_est        : estimated signal or systematic change values
  # max_clusters : maximum number of clusters considered for k-means
  # lag          : minimum distance between detected change points
  library(cluster)

  ## limit clusters to number of unique values
  # with rounding to avoid negligible differences (floating point representation)
  max_clusters <- min(length(unique(round(y_est, 7))), max_clusters)

  if (max_clusters >= 2) {
    gskmn <- cluster::clusGap(as.matrix(y_est),
      FUN = kmeans,
      nstart = 10,
      K.max = max_clusters,
      B = 500,
      d.power = 2,
      spaceH0 = "scaledPCA"
    )

    n_clust <- cluster::maxSE(
      f = gskmn$Tab[, 3],
      SE.f = gskmn$Tab[, 4],
      method = "Tibs2001SEmax",
      SE.factor = 1
    )

    final_kmeans <- kmeans(y_est,
      centers = n_clust,
      nstart = 10
    )

    change_points <- final_kmeans$cluster
    names(change_points) <- NULL
    change_points <- cumsum(rle(change_points)$lengths)
    change_points <- change_points[1:(length(change_points) - 1)]

    y <- c(0, change_points)

    ## enforce minimum distance between change points
    d <- which(diff(y) > lag)
    change_points <- change_points[d]
  } else {
    change_points <- integer(0)
  }

  return(change_points)
}


gam_helper_basis_dim_check <- function(mod, k.sample = 5000, k.rep = 200) {
  # Helper function to check effective degrees of freedom in a mgcv GAM.
  # source:
  # https://stackoverflow.com/questions/42042822/how-to-save-edf-from-mgcvgam-check-and-skip-the-plots
  # mod      : fitted mgcv::gam model
  # k.sample : number of observations to subsample for the check (default: 5000)
  # k.rep    : number of repetitions for the check (default: 200)

  mgcv:::k.check(mod, subsample = k.sample, n.rep = k.rep)
}


# Apply methods ----------------------------------------------------------------
# compute fitted values and change points if possible

f_arima <- function(simdata) {
  # Estimate systematic change using an ARIMA model.
  # simdata : data frame containing the observed measurements in column y

  library(forecast)

  ## fit ARIMA model to observed measurements
  fit <- forecast::auto.arima(simdata$y)

  ## estimate systematic change
  y_est <- simdata$y - as.vector(fit$residuals)

  return(list(
    fitted = y_est,
    changepoints = NULL
  ))
}


f_flsa <- function(simdata, groups = 1:4) {
  # Estimate systematic change using fused lasso signal approximation (FLSA).
  # simdata : data frame containing the observed measurements in column y
  # groups  : candidate numbers of segments used for model selection

  library(flsa)
  
  fit <- flsa::flsaTopDown(simdata$y, groups = groups)

  ## select model via k-means clustering
  km <- list()
  for (i in groups) {
    km[[i]] <- kmeans(as.matrix(fit$Solution[i, ]), i)
  }
  best.model <- which.min(sapply(km, function(x) x$tot.withinss))

  ## detect change points
  jump_points <- km[[best.model]]$cluster
  jump_points <- cumsum(rle(jump_points)$lengths)
  # omit the last jump point at the number of observations
  jump_points <- jump_points[seq_len(length(jump_points) - 1)]

  return(list(
    fitted = fit$Solution[best.model, ],
    changepoints = jump_points
  ))
}


f_gam <- function(simdata) {
  # Estimate systematic change using a GAM with thin plate regression splines
  # (default in mgcv). The function automatically adjusts the basis dimension (k)
  # if estimated degrees of freedom are close to the upper limit.
  # simdata : data frame containing observed measurements (columns: index, y)

  library(mgcv)
  k_num <- 10

  ## fit initial gam
  res <- mgcv::gam(y ~ s(index, k = k_num), data = simdata, method = "REML")

  ## check basis dimension
  kcheck <- gam_helper_basis_dim_check(res)

  ## increase basis dimension if EDF close to upper bound
  while (all(c("edf", "k'") %in% colnames(kcheck)) &&
    (round(kcheck[1, "edf"]) >= 0.7 * kcheck[1, "k'"] &
      k_num < nrow(simdata) - 10)) {
    k_num <- k_num + 10
    res <- mgcv::gam(y ~ s(index, k = k_num), data = simdata, method = "REML")
    kcheck <- f_gam_helper_basis_dim_check(res)
  }

  return(list(
    fitted = res$fitted.values,
    changepoints = NULL
  ))
}


f_lowess <- function(simdata, degree = 1) {
  # Estimate systematic change using LOWESS.
  # simdata : data frame containing observed measurements (columns: index, y)
  # degree  : degree of the local polynomial

  library(fANCOVA)

  ## fit LOESS (span selected automatically via AICc)
  fit <- loess.as(x = simdata$index, y = simdata$y, degree = degree)
  fitted_vals <- predict(fit, simdata$index)

  return(list(
    fitted = fitted_vals,
    changepoints = NULL
  ))
}


f_moving_average <- function(simdata) {
  # Estimate systematic change using a centered moving average. The window size
  # (order) is automatically selected using AICc.
  # simdata : data frame containing observed measurements (column: y)

  library(smooth)

  ## fit moving average with automatic order selection
  cma <- smooth::cma(simdata$y, order = NULL, silent = TRUE)

  return(list(
    fitted = as.numeric(cma$fitted),
    changepoints = NULL
  ))
}


f_pelt <- function(simdata, penalty = "MBIC") {
  # Estimate systematic change using the PELT change-point detection.
  # simdata : data frame containing the observed measurements in column y
  # penalty : penalty type for PELT

  library(changepoint)
  nobs <- nrow(simdata)

  ## detect changes in mean and variance
  cpt_results <- cpt.meanvar(simdata$y, method = "PELT", penalty = penalty)
  cpts_est <- cpt_results@cpts
  means_est <- cpt_results@param.est[["mean"]]

  ## estimate systematic change using segment means
  if (length(cpts_est) == 1) {
    y_est <- rep(means_est, times = nobs)
  } else {
    cpts_full <- c(0, cpts_est)
    y_est <- as.vector(unlist(mapply(function(start, end, mean) rep(mean, end - start),
      start = cpts_full[-length(cpts_full)],
      end = cpts_full[-1],
      mean = means_est
    )))
  }

  ## omit last point (nobs is always returned as a change point)
  cpts_est <- cpts_est[-length(cpts_est)]

  return(list(
    fitted = y_est,
    changepoints = cpts_est
  ))
}


f_piecewise_reg <- function(simdata, n_cp = 1:3) {
  # Estimate systematic change using piecewise linear regression.
  # simdata : data frame containing observed measurements (columns: index, y)
  # n_cp    : candidate numbers of change points for model selection

  library(segmented)

  sm <- list()
  lm_model <- lm(y ~ index, data = simdata)

  ## fit segmented models with different numbers of breakpoints
  for (breakpoints in n_cp) {
    sm[[breakpoints]] <- try(segmented(lm_model, npsi = breakpoints), silent = TRUE)
  }

  ## keep only successful models and linear model (no breakpoints)
  successful <- sapply(sm, function(x) inherits(x, "segmented"))
  sm <- append(sm[successful], list(lm_model), after = 0)

  ## select model with lowest BIC
  bic <- sapply(sm, BIC)
  fit <- sm[[which.min(bic)]]

  ## extract change points
  if (is.null(fit$id.group)) {
    # no jump point identified resulting in linear fit
    jump_points <- integer(0)
  } else {
    groups <- fit$id.group
    jump_points <- cumsum(rle(groups)$lengths)
    # omit jump point at the end of the data
    jump_points <- jump_points[1:(length(jump_points) - 1)]
  }

  return(list(
    fitted = fit$fitted.values,
    changepoints = jump_points
  ))
}


# Wrapper batchtools -----------------------------------------------------------

bt_arima <- function(data, job, instance, max_clusters = 4, lag = length(y_est) * 0.05, ...) {
  # Wrapper to apply ARIMA in batchtools simulation framework.
  # data, job     : required by batchtools
  # instance      : single simulation dataset (data frame with y, syschange, cpts)
  # max_clusters  : maximum clusters for change-point detection
  # lag           : minimum distance between detected change points

  start <- Sys.time()

  ## estimate systematic change
  results <- f_arima(simdata = instance)
  y_est <- results$fitted

  ## detect change points
  cpts_est <- find_change_points(y_est,
    max_clusters = max_clusters,
    lag = lag
  )

  ## compute estimates
  estimators <- compute_estimators(
    y = instance$y,
    syschange = y_est,
    cpts = cpts_est
  )
  end <- Sys.time()

  runtime <- difftime(end, start, units = "secs")
  colnames(estimators) <- paste0(colnames(estimators), "_est")

  ## compute true values of the estimators
  true_values <- compute_estimators(
    y = instance$y,
    syschange = instance$syschange,
    cpts = instance$cpts
  )

  return(cbind(estimators, true_values, runtime))
}


bt_flsa <- function(data, job, instance, groups = 1:4, ...) {
  # Wrapper to apply FLSA in batchtools simulation framework.
  # data, job  : required by batchtools
  # instance   : single simulation dataset (data frame with y, syschange, cpts)
  # groups     : candidate numbers of segments used for model selection

  start <- Sys.time()

  ## estimate systematic change
  results <- f_flsa(simdata = instance, groups = groups)
  y_est <- results$fitted

  ## extract change points
  cpts_est <- results$changepoints

  ## compute estimates
  estimators <- compute_estimators(
    y = instance$y,
    syschange = y_est,
    cpts = cpts_est
  )
  end <- Sys.time()

  runtime <- difftime(end, start, units = "secs")
  colnames(estimators) <- paste0(colnames(estimators), "_est")

  ## compute true values of the estimators
  true_values <- compute_estimators(
    y = instance$y,
    syschange = instance$syschange,
    cpts = instance$cpts
  )

  return(cbind(estimators, true_values, runtime))
}


bt_gam <- function(data, job, instance, max_clusters = 4, lag = length(y_est) * 0.05, ...) {
  # Wrapper to apply GAM in batchtools simulation framework.
  # data, job     : required by batchtools
  # instance      : single simulation dataset (data frame with y, syschange, cpts)
  # max_clusters  : maximum clusters for change-point detection
  # lag           : minimum distance between detected change points

  start <- Sys.time()

  ## estimate systematic change
  results <- f_gam(simdata = instance)
  y_est <- results$fitted

  ## detect change points
  cpts_est <- find_change_points(y_est,
    max_clusters = max_clusters,
    lag = lag
  )
  ## compute estimates
  estimators <- compute_estimators(
    y = instance$y,
    syschange = y_est,
    cpts = cpts_est
  )
  end <- Sys.time()

  runtime <- difftime(end, start, units = "secs")
  colnames(estimators) <- paste0(colnames(estimators), "_est")

  ## compute true values of the estimators
  true_values <- compute_estimators(
    y = instance$y,
    syschange = instance$syschange,
    cpts = instance$cpts
  )

  return(cbind(estimators, true_values, runtime))
}


bt_lowess <- function(data, job, instance, degree = 1, max_clusters = 4, lag = length(y_est) * 0.05, ...) {
  # Wrapper to apply LOWESS in batchtools simulation framework.
  # data, job     : required by batchtools
  # instance      : single simulation dataset (data frame with y, syschange, cpts)
  # degree        : degree of the local polynomial
  # max_clusters  : maximum clusters for change-point detection
  # lag           : minimum distance between detected change points

  start <- Sys.time()

  ## estimate systematic change
  results <- f_lowess_aicc(simdata = instance, degree = degree)
  y_est <- results$fitted

  ## detect change points
  cpts_est <- find_change_points(y_est,
    max_clusters = max_clusters,
    lag = lag
  )

  ## compute estimates
  estimators <- compute_estimators(
    y = instance$y,
    syschange = y_est,
    cpts = cpts_est
  )
  end <- Sys.time()

  runtime <- difftime(end, start, units = "secs")
  colnames(estimators) <- paste0(colnames(estimators), "_est")

  ## compute true values of the estimators
  true_values <- compute_estimators(
    y = instance$y,
    syschange = instance$syschange,
    cpts = instance$cpts
  )

  return(cbind(estimators, true_values, runtime))
}


bt_moving_average <- function(data, job, instance, max_clusters = 4, lag = length(y_est) * 0.05, ...) {
  # Wrapper to apply moving average in batchtools simulation framework.
  # data, job     : required by batchtools
  # instance      : single simulation dataset (data frame with y, syschange, cpts)
  # max_clusters  : maximum clusters for change-point detection
  # lag           : minimum distance between detected change points

  start <- Sys.time()

  ## estimate systematic change
  results <- f_moving_average(simdata = instance)
  y_est <- results$fitted

  ## detect change points
  cpts_est <- find_change_points(y_est,
    max_clusters = max_clusters,
    lag = lag
  )

  ## compute estimates
  estimators <- compute_estimators(
    y = instance$y,
    syschange = y_est,
    cpts = cpts_est
  )
  end <- Sys.time()

  runtime <- difftime(end, start, units = "secs")
  colnames(estimators) <- paste0(colnames(estimators), "_est")

  ## compute true values of the estimators
  true_values <- compute_estimators(
    y = instance$y,
    syschange = instance$syschange,
    cpts = instance$cpts
  )

  return(cbind(estimators, true_values, runtime))
}


bt_pelt <- function(data, job, instance, penalty = "MBIC", ...) {
  # Wrapper to apply PELT in batchtools simulation framework.
  # data, job  : required by batchtools
  # instance   : single simulation dataset (data frame with y, syschange, cpts)
  # penalty    : penalty type for PELT

  start <- Sys.time()

  ## estimate systematic change
  results <- f_pelt(simdata = instance, penalty = penalty)
  y_est <- results$fitted

  ## extract change points
  cpts_est <- results$changepoints

  ## compute estimates
  estimators <- compute_estimators(
    y = instance$y,
    syschange = y_est,
    cpts = cpts_est
  )
  end <- Sys.time()

  runtime <- difftime(end, start, units = "secs")
  colnames(estimators) <- paste0(colnames(estimators), "_est")

  ## compute true values of the estimators
  true_values <- compute_estimators(
    y = instance$y,
    syschange = instance$syschange,
    cpts = instance$cpts
  )

  return(cbind(estimators, true_values, runtime))
}


bt_piecewise_reg <- function(data, job, instance, n_cp = 1:3, ...) {
  # Wrapper to apply piecewise regression in batchtools simulation framework.
  # data, job  : required by batchtools
  # instance   : single simulation dataset (data frame with y, syschange, cpts)
  # n_cp       : candidate numbers of change points for model selection

  start <- Sys.time()

  ## estimate systematic change
  results <- f_piecewise_reg(simdata = instance, n_cp = n_cp)
  y_est <- results$fitted

  ## extract change points
  cpts_est <- results$changepoints

  ## compute estimates
  estimators <- compute_estimators(
    y = instance$y,
    syschange = y_est,
    cpts = cpts_est
  )
  end <- Sys.time()

  runtime <- difftime(end, start, units = "secs")
  colnames(estimators) <- paste0(colnames(estimators), "_est")

  ## compute true values of the estimators
  true_values <- compute_estimators(
    y = instance$y,
    syschange = instance$syschange,
    cpts = instance$cpts
  )

  return(cbind(estimators, true_values, runtime))
}
