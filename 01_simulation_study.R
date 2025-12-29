## Implementation of the simulation study

# load packages & functions
library(batchtools)
library(changepoint)
library(cluster)
library(dplyr)
library(fANCOVA)
library(flsa)
library(forecast)
library(mgcv)
library(segmented)
library(smooth)

# set seed
set.seed(42)


# Registry ---------------------------------------------------------------------
reg <- makeExperimentRegistry(file.dir = "../simulation_reg", seed = 1)
# reg <- loadRegistry(file.dir="../simulation_reg", writeable=FALSE)


# Problems ---------------------------------------------------------------------
source("functions/problems.R")
addProblem(name = "prob_norm_nochange", fun = f_simulate_norm_nochange, seed = 43)
addProblem(name = "prob_norm_jump", fun = f_simulate_norm_jump, seed = 44)
addProblem(name = "prob_norm_linear", fun = f_simulate_norm_linear, seed = 45)
addProblem(name = "prob_norm_quadratic", fun = f_simulate_norm_quadratic, seed = 46)
addProblem(name = "prob_norm_recurrent", fun = f_simulate_norm_recurrent, seed = 47)
addProblem(name = "prob_lognorm_nochange", fun = f_simulate_lognorm_nochange, seed = 49)
addProblem(name = "prob_lognorm_jump", fun = f_simulate_lognorm_jump, seed = 50)
addProblem(name = "prob_lognorm_linear", fun = f_simulate_lognorm_linear, seed = 51)
addProblem(name = "prob_lognorm_quadratic", fun = f_simulate_lognorm_quadratic, seed = 52)
addProblem(name = "prob_lognorm_recurrent", fun = f_simulate_lognorm_recurrent, seed = 53)


# Algorithms -------------------------------------------------------------------
source("functions/algorithms.R")
addAlgorithm(name = "alg_arima", fun = bt_arima)
addAlgorithm(name = "alg_flsa", fun = bt_flsa)
addAlgorithm(name = "alg_gam", fun = bt_gam)
addAlgorithm(name = "alg_lowess", fun = bt_lowess)
addAlgorithm(name = "alg_moving_avg", fun = bt_moving_average)
addAlgorithm(name = "alg_pelt", fun = bt_pelt)
addAlgorithm(name = "alg_piecewise_reg", fun = bt_piecewise_reg)


# Experiments ------------------------------------------------------------------
algo_design <- list(
  alg_arima = expand.grid(stringsAsFactors = FALSE),
  alg_flsa = expand.grid(stringsAsFactors = FALSE),
  alg_gam = expand.grid(stringsAsFactors = FALSE),
  alg_lowess = expand.grid(stringsAsFactors = FALSE),
  alg_moving_avg = expand.grid(stringsAsFactors = FALSE),
  alg_pelt = expand.grid(stringsAsFactors = FALSE),
  alg_piecewise_reg = expand.grid(stringsAsFactors = FALSE)
)

## Testing
# prob_design <- list(prob_norm_notrend=expand.grid(nobs=500, snr=1),
#                     prob_norm_jump=expand.grid(nobs=500, snr=1, dmaxf=1),
#                     prob_norm_linear=expand.grid(nobs=500, snr=1, dmaxf=1),
#                     prob_norm_quadratic=expand.grid(nobs=500, snr=1, dmaxf=1),
#                     prob_norm_recurrent=expand.grid(nobs=500, snr=1, dmaxf=1),
#                     prob_lognorm_notrend=expand.grid(nobs=500, snr=1),
#                     prob_lognorm_jump=expand.grid(nobs=500, snr=1, dmaxf=1),
#                     prob_lognorm_linear=expand.grid(nobs=500, snr=1, dmaxf=1),
#                     prob_lognorm_quadratic=expand.grid(nobs=500, snr=1, dmaxf=1),
#                     prob_lognorm_recurrent=expand.grid(nobs=500, snr=1, dmaxf=1))
# addExperiments(prob.designs=prob_design, algo.designs=algo_design, repls=1)

## Interval 1
nobs <- seq(30, 50, by = 1)
nsim <- 10
prob_design <- list(
  prob_norm_notrend = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2)),
  prob_norm_jump = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_linear = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_quadratic = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_recurrent = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_notrend = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2)),
  prob_lognorm_jump = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_linear = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_quadratic = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_recurrent = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2))
)
addExperiments(prob.designs = prob_design, algo.designs = algo_design, repls = nsim)

## Interval 2
nobs <- seq(51, 100, by = 1)
nsim <- 2
prob_design <- list(
  prob_norm_notrend = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2)),
  prob_norm_jump = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_linear = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_quadratic = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_recurrent = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_notrend = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2)),
  prob_lognorm_jump = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_linear = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_quadratic = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_recurrent = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2))
)
addExperiments(prob.designs = prob_design, algo.designs = algo_design, repls = nsim)

## Interval 3
nobs <- seq(101, 200, by = 1)
nsim <- 1
prob_design <- list(
  prob_norm_notrend = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2)),
  prob_norm_jump = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_linear = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_quadratic = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_recurrent = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_notrend = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2)),
  prob_lognorm_jump = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_linear = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_quadratic = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_recurrent = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2))
)
addExperiments(prob.designs = prob_design, algo.designs = algo_design, repls = nsim)

## Interval 4
nobs <- seq(205, 500, by = 5)
nsim <- 1
prob_design <- list(
  prob_norm_notrend = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2)),
  prob_norm_jump = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_linear = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_quadratic = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_recurrent = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_notrend = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2)),
  prob_lognorm_jump = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_linear = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_quadratic = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_recurrent = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2))
)
addExperiments(prob.designs = prob_design, algo.designs = algo_design, repls = nsim)

## Interval 5
nobs <- seq(510, 1000, by = 10)
nsim <- 1
prob_design <- list(
  prob_norm_notrend = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2)),
  prob_norm_jump = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_linear = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_quadratic = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_norm_recurrent = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_notrend = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2)),
  prob_lognorm_jump = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_linear = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_quadratic = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2)),
  prob_lognorm_recurrent = expand.grid(nobs = nobs, snr = c(0.33, 0.66, 1, 2), dmaxf = c(0.25, 0.5, 1, 2))
)
addExperiments(prob.designs = prob_design, algo.designs = algo_design, repls = nsim)


# Tests ------------------------------------------------------------------------
# Test one job for each algorithm

## ARIMA
testJob(id = head(findExperiments(prob.name = "prob_norm_notrend", algo.name = "alg_arima"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_jump", algo.name = "alg_arima"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_linear", algo.name = "alg_arima"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_quadratic", algo.name = "alg_arima"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_recurrent", algo.name = "alg_arima"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_notrend", algo.name = "alg_arima"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_jump", algo.name = "alg_arima"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_linear", algo.name = "alg_arima"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_quadratic", algo.name = "alg_arima"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_recurrent", algo.name = "alg_arima"), 1))

## FLSA
testJob(id = head(findExperiments(prob.name = "prob_norm_notrend", algo.name = "alg_flsa"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_jump", algo.name = "alg_flsa"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_linear", algo.name = "alg_flsa"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_quadratic", algo.name = "alg_flsa"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_recurrent", algo.name = "alg_flsa"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_notrend", algo.name = "alg_flsa"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_jump", algo.name = "alg_flsa"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_linear", algo.name = "alg_flsa"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_quadratic", algo.name = "alg_flsa"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_recurrent", algo.name = "alg_flsa"), 1))

## GAM
testJob(id = head(findExperiments(prob.name = "prob_norm_notrend", algo.name = "alg_gam"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_jump", algo.name = "alg_gam"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_linear", algo.name = "alg_gam"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_quadratic", algo.name = "alg_gam"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_recurrent", algo.name = "alg_gam"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_notrend", algo.name = "alg_gam"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_jump", algo.name = "alg_gam"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_linear", algo.name = "alg_gam"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_quadratic", algo.name = "alg_gam"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_recurrent", algo.name = "alg_gam"), 1))

## LOWESS
testJob(id = head(findExperiments(prob.name = "prob_norm_notrend", algo.name = "alg_lowess"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_jump", algo.name = "alg_lowess"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_linear", algo.name = "alg_lowess"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_quadratic", algo.name = "alg_lowess"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_recurrent", algo.name = "alg_lowess"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_notrend", algo.name = "alg_lowess"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_jump", algo.name = "alg_lowess"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_linear", algo.name = "alg_lowess"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_quadratic", algo.name = "alg_lowess"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_recurrent", algo.name = "alg_lowess"), 1))

## Moving average
testJob(id = head(findExperiments(prob.name = "prob_norm_notrend", algo.name = "alg_moving_avg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_jump", algo.name = "alg_moving_avg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_linear", algo.name = "alg_moving_avg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_quadratic", algo.name = "alg_moving_avg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_recurrent", algo.name = "alg_moving_avg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_notrend", algo.name = "alg_moving_avg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_jump", algo.name = "alg_moving_avg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_linear", algo.name = "alg_moving_avg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_quadratic", algo.name = "alg_moving_avg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_recurrent", algo.name = "alg_moving_avg"), 1))

## PELT
testJob(id = head(findExperiments(prob.name = "prob_norm_notrend", algo.name = "alg_pelt"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_jump", algo.name = "alg_pelt"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_linear", algo.name = "alg_pelt"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_quadratic", algo.name = "alg_pelt"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_recurrent", algo.name = "alg_pelt"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_notrend", algo.name = "alg_pelt"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_jump", algo.name = "alg_pelt"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_linear", algo.name = "alg_pelt"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_quadratic", algo.name = "alg_pelt"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_recurrent", algo.name = "alg_pelt"), 1))

## Piecewise regression
testJob(id = head(findExperiments(prob.name = "prob_norm_notrend", algo.name = "alg_piecewise_reg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_jump", algo.name = "alg_piecewise_reg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_linear", algo.name = "alg_piecewise_reg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_quadratic", algo.name = "alg_piecewise_reg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_norm_recurrent", algo.name = "alg_piecewise_reg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_notrend", algo.name = "alg_piecewise_reg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_jump", algo.name = "alg_piecewise_reg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_linear", algo.name = "alg_piecewise_reg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_quadratic", algo.name = "alg_piecewise_reg"), 1))
testJob(id = head(findExperiments(prob.name = "prob_lognorm_recurrent", algo.name = "alg_piecewise_reg"), 1))


# Submit -----------------------------------------------------------------------
ids <- findNotDone()
submitJobs(ids)


# Get results ------------------------------------------------------------------
# the results are automatically saved in the registry directory
getStatus()

# check for errors
errs <- findErrors()
err_mes <- getErrorMessages(errs)
err_job_table <- getJobTable(errs)

## modify results for convenience
res <- flatten(ijoin(getJobPars(), reduceResultsDataTable()))
res$algorithm <- factor(res$algorithm, levels = sort(unique(res$algorithm)))

# extract information from problem name
res$dmaxf[grepl("nochange", res$problem)] <- 0
res$distribution <- sapply(strsplit(res$problem, "_"), function(x) x[2])
res$pattern <- sapply(strsplit(res$problem, "_"), function(x) x[3])

# create sample size categories
res$nobs_cat <- cut(
  res$nobs,
  breaks = c(30, 50, 100, 200, 500, 1000),
  labels = c("30–50", "51–100", "101–200", "201–500", "501–1000"),
  include.lowest = TRUE,
  right = TRUE
)

# compute difference between true and estimated values
res$range_dif <- res$range_est - res$range
res$var_dif <- res$var_est - res$var
res$madm_dif <- res$madm_est - res$madm
res$n_cpts_dif <- res$n_cpts_est - res$n_cpts

# add prob.id
job_table <- getJobTable()
job_table <- data.frame(
  job.id = job_table$job.id,
  repl = job_table$repl
)
res <- res %>%
  left_join(job_table, by = join_by(job.id)) %>%
  group_by(problem, nobs, snr, dmaxf, repl) %>%
  mutate(prob.id = cur_group_id()) %>%
  ungroup()

saveRDS(res, "/results/res.rds")
