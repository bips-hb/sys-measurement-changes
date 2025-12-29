## Functions to create simulated data


# Simulation with normally distributed data ------------------------------------

f_simulate_norm_nochange <- function(data, job, nobs = 100, snr = 1) {
  # Simulate normally distributed measurements without systematic change.
  # data, job : placeholders for the simulation framework
  # nobs      : number of observations
  # snr       : signal-to-noise ratio controlling random error variance

  mean <- 0
  sd <- 1
  sd_unsys_error <- sd / sqrt(snr)

  ## measurements
  # index for the measurements
  index <- seq(from = 1, to = nobs, by = 1)
  # proportion for the measurements
  prop <- index / nobs
  # true values
  x <- rnorm(n = nobs, mean = mean, sd = sd)

  ## unsystematic error
  unsys_error <- rnorm(n = nobs, mean = 0, sd = sd_unsys_error)

  ## measurements with error
  return(data.frame(
    index = index,
    syschange = rep(0, times = nobs),
    cpts = rep(0, times = nobs),
    y = x + unsys_error
  ))
}


f_simulate_norm_jump <- function(data, job, nobs = 100, snr = 1, dmaxf = 0.8) {
  # Simulate normally distributed measurements with a jump.
  # data, job : placeholders for the simulation framework
  # nobs      : number of observations
  # snr       : signal-to-noise ratio controlling random error variance
  # dmaxf     : magnitude of the systematic change in relation to the standard
  #             deviation of the true outcome

  mean <- 0
  sd <- 1
  p1 <- 0.25
  sd_unsys_error <- sd / sqrt(snr)
  # Attention: for .5 R rounds to even. This is not what we want here.
  # Instead, we use the floor function to get the desired behavior.
  p1_index <- floor(p1 * nobs + 0.5)
  dmax <- dmaxf * sd

  ## measurements
  # index for the measurements
  index <- seq(from = 1, to = nobs, by = 1)
  # proportion for the measurements
  p <- index / nobs
  # true values
  x <- rnorm(n = nobs, mean = mean, sd = sd)

  ## unsystematic error
  unsys_error <- rnorm(n = nobs, mean = 0, sd = sd_unsys_error)

  ## systematic change
  sys_change <- rep(0, times = nobs)
  sys_change[c(p1_index:nobs)] <- dmax

  ## changepoints
  cpts <- rep(0, times = nobs)
  cpts[p1_index] <- 1

  ## measurements with error
  return(data.frame(
    index = index,
    syschange = sys_change,
    cpts = cpts,
    y = x + unsys_error + sys_change
  ))
}


f_simulate_norm_linear <- function(data, job, nobs = 100, snr = 1, dmaxf = 0.8) {
  # Simulate normally distributed measurements with a linear change.
  # data, job : placeholders for the simulation framework
  # nobs      : number of observations
  # snr       : signal-to-noise ratio controlling random error variance
  # dmaxf     : magnitude of the systematic change in relation to the standard
  #             deviation of the true outcome

  mean <- 0
  sd <- 1
  p1 <- 0.25
  p2 <- 0.5
  sd_unsys_error <- sd / sqrt(snr)
  # Attention: for .5 R rounds to even. This is not what we want here.
  # Instead, we use the floor function to get the desired behavior.
  p1_index <- floor(p1 * nobs + 0.5)
  p2_index <- floor(p2 * nobs + 0.5)
  dmax <- dmaxf * sd
  slope <- dmax / (nobs * (p2 - p1))
  intercept <- -(dmax * p1) / (p2 - p1)

  ## measurements
  # index for the measurements
  index <- seq(from = 1, to = nobs, by = 1)
  # proportion for the measurements
  p <- index / nobs
  # true values
  x <- rnorm(n = nobs, mean = mean, sd = sd)

  ## unsystematic error
  unsys_error <- rnorm(n = nobs, mean = 0, sd = sd_unsys_error)

  ## systematic change
  sys_change <- rep(0, times = nobs)
  sys_change[c(p1_index:p2_index - 1)] <- index[c(p1_index:p2_index)] * slope + intercept
  sys_change[c(p2_index:nobs)] <- dmax

  ## changepoints
  cpts <- rep(0, times = nobs)
  cpts[c(p1_index, p2_index)] <- 1

  ## measurements with error
  return(data.frame(
    index = index,
    syschange = sys_change,
    cpts = cpts,
    y = x + unsys_error + sys_change
  ))
}


f_simulate_norm_quadratic <- function(data, job, nobs = 100, snr = 1, dmaxf = 0.8) {
  # Simulate normally distributed measurements with a quadratic change.
  # data, job : placeholders for the simulation framework
  # nobs      : number of observations
  # snr       : signal-to-noise ratio controlling random error variance
  # dmaxf     : magnitude of the systematic change in relation to the standard
  #             deviation of the true outcome

  mean <- 0
  sd <- 1
  p1 <- 0.25
  p2 <- 0.5
  p3 <- p2 - p1 / sqrt(2) + p2 / sqrt(2)
  sd_unsys_error <- sd / sqrt(snr)
  # Attention: for .5 R rounds to even. This is not what we want here.
  # Instead, we use the floor function to get the desired behavior.
  p1_index <- floor(p1 * nobs + 0.5)
  p2_index <- floor(p2 * nobs + 0.5)
  p3_index <- floor(p3 * nobs + 0.5)
  dmax <- dmaxf * sd

  factor <- -dmax / ((p1 * nobs - p2 * nobs)^2)

  ## measurements
  # index for the measurements
  index <- seq(from = 1, to = nobs, by = 1)
  # proportion for the measurements
  p <- index / nobs
  # true values
  x <- rnorm(n = nobs, mean = mean, sd = sd)

  ## unsystematic error
  unsys_error <- rnorm(n = nobs, mean = 0, sd = sd_unsys_error)

  ## systematic change
  sys_change <- rep(0, times = nobs)
  sys_change[c(p1_index:p3_index - 1)] <- factor * (index[c(p1_index:p3_index)] - p2 * nobs)^2 + dmax
  sys_change[c(p3_index:nobs)] <- 0.5 * dmax

  ## changepoints
  cpts <- rep(0, times = nobs)
  cpts[c(p1_index, p2_index, p3_index)] <- 1

  ## measurements with error
  return(data.frame(
    index = index,
    syschange = sys_change,
    cpts = cpts,
    y = x + unsys_error + sys_change
  ))
}


f_simulate_norm_recurrent <- function(data, job, nobs = 100, snr = 1, dmaxf = 0.8) {
  # Simulate normally distributed measurements with recurrent jumps.
  # data, job : placeholders for the simulation framework
  # nobs      : number of observations
  # snr       : signal-to-noise ratio controlling random error variance
  # dmaxf     : magnitude of the systematic change in relation to the standard
  #             deviation of the true outcome

  mean <- 0
  sd <- 1
  p1 <- 0.25
  p2 <- 0.5
  p3 <- p2 - p1 / sqrt(2) + p2 / sqrt(2)
  sd_unsys_error <- sd / sqrt(snr)
  # Attention: for .5 R rounds to even. This is not what we want here.
  # Instead, we use the floor function to get the desired behavior.
  p1_index <- floor(p1 * nobs + 0.5)
  p2_index <- floor(p2 * nobs + 0.5)
  p3_index <- floor(p3 * nobs + 0.5)
  dmax <- dmaxf * sd

  ## measurements
  # index for the measurements
  index <- seq(from = 1, to = nobs, by = 1)
  # proportion for the measurements
  p <- index / nobs
  # true values
  x <- rnorm(n = nobs, mean = mean, sd = sd)

  ## unsystematic error
  unsys_error <- rnorm(n = nobs, mean = 0, sd = sd_unsys_error)

  ## systematic change
  sys_change <- rep(0, times = nobs)
  sys_change[c(p1_index:(p2_index - 1))] <- dmax
  sys_change[c(p3_index:nobs)] <- dmax

  ## changepoints
  cpts <- rep(0, times = nobs)
  cpts[c(p1_index, p2_index, p3_index)] <- 1

  ## measurements with error
  return(data.frame(
    index = index,
    syschange = sys_change,
    cpts = cpts,
    y = x + unsys_error + sys_change
  ))
}


# Simulation with log-normally distributed data --------------------------------

f_simulate_lognorm_nochange <- function(data, job, nobs = 100, snr = 1) {
  # Simulate log-normally distributed measurements without systematic change.
  # data, job : placeholders for the simulation framework
  # nobs      : number of observations
  # snr       : signal-to-noise ratio controlling random error variance

  mean <- 0
  sd <- 1
  sd_unsys_error <- sd / sqrt(snr)

  ## measurements
  # index for the measurements
  index <- seq(from = 1, to = nobs, by = 1)
  # proportion for the measurements
  prop <- index / nobs
  # true values
  x <- rlnorm(n = nobs, meanlog = mean, sdlog = sd)

  ## unsystematic error
  unsys_error <- rnorm(n = nobs, mean = 0, sd = sd_unsys_error)

  ## measurements with error
  return(data.frame(
    index = index,
    syschange = rep(0, times = nobs),
    cpts = rep(0, times = nobs),
    y = x + unsys_error
  ))
}


f_simulate_lognorm_jump <- function(data, job, nobs = 100, snr = 1, dmaxf = 0.8) {
  # Simulate log-normally distributed measurements with a jump.
  # data, job : placeholders for the simulation framework
  # nobs      : number of observations
  # snr       : signal-to-noise ratio controlling random error variance
  # dmaxf     : magnitude of the systematic change in relation to the standard
  #             deviation of the true outcome

  mean <- 0
  sd <- 1
  p1 <- 0.25
  sd_unsys_error <- sd / sqrt(snr)
  # Attention: for .5 R rounds to even. This is not what we want here.
  # Instead, we use the floor function to get the desired behavior.
  p1_index <- floor(p1 * nobs + 0.5)
  dmax <- dmaxf * sd

  ## measurements
  # index for the measurements
  index <- seq(from = 1, to = nobs, by = 1)
  # proportion for the measurements
  p <- index / nobs
  # true values
  x <- rlnorm(n = nobs, meanlog = mean, sdlog = sd)

  ## unsystematic error
  unsys_error <- rnorm(n = nobs, mean = 0, sd = sd_unsys_error)

  ## systematic change
  sys_change <- rep(0, times = nobs)
  sys_change[c(p1_index:nobs)] <- dmax

  ## changepoints
  cpts <- rep(0, times = nobs)
  cpts[p1_index] <- 1

  ## measurements with error
  return(data.frame(
    index = index,
    syschange = sys_change,
    cpts = cpts,
    y = x + unsys_error + sys_change
  ))
}


f_simulate_lognorm_linear <- function(data, job, nobs = 100, snr = 1, dmaxf = 0.8) {
  # Simulate log-normally distributed measurements with a linear change.
  # data, job : placeholders for the simulation framework
  # nobs      : number of observations
  # snr       : signal-to-noise ratio controlling random error variance
  # dmaxf     : magnitude of the systematic change in relation to the standard
  #             deviation of the true outcome

  mean <- 0
  sd <- 1
  p1 <- 0.25
  p2 <- 0.5
  sd_unsys_error <- sd / sqrt(snr)
  # Attention: for .5 R rounds to even. This is not what we want here.
  # Instead, we use the floor function to get the desired behavior.
  p1_index <- floor(p1 * nobs + 0.5)
  p2_index <- floor(p2 * nobs + 0.5)
  dmax <- dmaxf * sd
  slope <- dmax / (nobs * (p2 - p1))
  intercept <- -(dmax * p1) / (p2 - p1)

  ## measurements
  # index for the measurements
  index <- seq(from = 1, to = nobs, by = 1)
  # proportion for the measurements
  p <- index / nobs
  # true values
  x <- rlnorm(n = nobs, meanlog = mean, sdlog = sd)

  ## unsystematic error
  unsys_error <- rnorm(n = nobs, mean = 0, sd = sd_unsys_error)

  ## systematic change
  sys_change <- rep(0, times = nobs)
  sys_change[c(p1_index:p2_index - 1)] <- index[c(p1_index:p2_index)] * slope + intercept
  sys_change[c(p2_index:nobs)] <- dmax

  ## changepoints
  cpts <- rep(0, times = nobs)
  cpts[c(p1_index, p2_index)] <- 1

  ## measurements with error
  return(data.frame(
    index = index,
    syschange = sys_change,
    cpts = cpts,
    y = x + unsys_error + sys_change
  ))
}


f_simulate_lognorm_quadratic <- function(data, job, nobs = 100, snr = 1, dmaxf = 0.8) {
  # Simulate log-normally distributed measurements with a quadratic change.
  # data, job : placeholders for the simulation framework
  # nobs      : number of observations
  # snr       : signal-to-noise ratio controlling random error variance
  # dmaxf     : magnitude of the systematic change in relation to the standard
  #             deviation of the true outcome

  mean <- 0
  sd <- 1
  p1 <- 0.25
  p2 <- 0.5
  p3 <- p2 - p1 / sqrt(2) + p2 / sqrt(2)
  sd_unsys_error <- sd / sqrt(snr)
  # Attention: for .5 R rounds to even. This is not what we want here.
  # Instead, we use the floor function to get the desired behavior.
  p1_index <- floor(p1 * nobs + 0.5)
  p2_index <- floor(p2 * nobs + 0.5)
  p3_index <- floor(p3 * nobs + 0.5)
  dmax <- dmaxf * sd

  factor <- -dmax / ((p1 * nobs - p2 * nobs)^2)

  ## measurements
  # index for the measurements
  index <- seq(from = 1, to = nobs, by = 1)
  # proportion for the measurements
  p <- index / nobs
  # true values
  x <- rlnorm(n = nobs, meanlog = mean, sdlog = sd)

  ## unsystematic error
  unsys_error <- rnorm(n = nobs, mean = 0, sd = sd_unsys_error)

  ## systematic change
  sys_change <- rep(0, times = nobs)
  sys_change[c(p1_index:p3_index - 1)] <- factor * (index[c(p1_index:p3_index)] - p2 * nobs)^2 + dmax
  sys_change[c(p3_index:nobs)] <- 0.5 * dmax

  ## changepoints
  cpts <- rep(0, times = nobs)
  cpts[c(p1_index, p2_index, p3_index)] <- 1

  ## measurements with error
  return(data.frame(
    index = index,
    syschange = sys_change,
    cpts = cpts,
    y = x + unsys_error + sys_change
  ))
}


f_simulate_lognorm_recurrent <- function(data, job, nobs = 100, snr = 1, dmaxf = 0.8) {
  # Simulate log-normally distributed measurements with recurrent jumps.
  # data, job : placeholders for the simulation framework
  # nobs      : number of observations
  # snr       : signal-to-noise ratio controlling random error variance
  # dmaxf     : magnitude of the systematic change in relation to the standard
  #             deviation of the true outcome

  mean <- 0
  sd <- 1
  p1 <- 0.25
  p2 <- 0.5
  p3 <- p2 - p1 / sqrt(2) + p2 / sqrt(2)
  sd_unsys_error <- sd / sqrt(snr)
  # Attention: for .5 R rounds to even. This is not what we want here.
  # Instead, we use the floor function to get the desired behavior.
  p1_index <- floor(p1 * nobs + 0.5)
  p2_index <- floor(p2 * nobs + 0.5)
  p3_index <- floor(p3 * nobs + 0.5)
  dmax <- dmaxf * sd

  ## measurements
  # index for the measurements
  index <- seq(from = 1, to = nobs, by = 1)
  # proportion for the measurements
  p <- index / nobs
  # true values
  x <- rlnorm(n = nobs, meanlog = mean, sdlog = sd)

  ## unsystematic error
  unsys_error <- rnorm(n = nobs, mean = 0, sd = sd_unsys_error)

  ## systematic change
  sys_change <- rep(0, times = nobs)
  sys_change[c(p1_index:(p2_index - 1))] <- dmax
  sys_change[c(p3_index:nobs)] <- dmax

  ## changepoints
  cpts <- rep(0, times = nobs)
  cpts[c(p1_index, p2_index, p3_index)] <- 1

  ## measurements with error
  return(data.frame(
    index = index,
    syschange = sys_change,
    cpts = cpts,
    y = x + unsys_error + sys_change
  ))
}
