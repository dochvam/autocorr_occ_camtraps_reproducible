#########################################################################
#' Author: Benjamin R. Goldstein
#' run_simulation2.R
#' 
#' This file executes simulation 2 from the manuscript and produces output
#' files giving results from those simulations. Simulation 2 looks at how
#' the impact of autocorrelation varies with detection window and gap.
#########################################################################

library(tidyverse)
library(nimble)
library(unmarked)
library(AICcmodavg)
library(nimbleEcology)
library(parabar)

# Arbitrary seed used for the manuscript
set.seed(5995235)

# Read in helper functions for simulation and model fitting
source("code_helper/helper_fn_clustered.R")

# Simulation parameters
nsim <- 1#000
par_method <- "parabar"
ncores <- 44
tag <- "" # optional suffix to append to output filenames
datasummary_only <- FALSE


#### Set up simulation conditions ####

# Fixed sampling and system conditions
parameter_grid <- expand.grid(
  logit_psi_mu = logit(c(0.4)),
  overall_det_prob = c(0.3),
  det_breakdown = 1,
  autocorr_str = c(1, 10)
) %>%
  as.data.frame() %>% 
  filter(overall_det_prob * det_breakdown < 1 &
           overall_det_prob / det_breakdown < 1)

# Fit varying windows and gaps
effort_scenarios <- data.frame(
  window = rep(1:35, each = 4),
  gap = 0:3
) %>% 
  cross_join(expand.grid(
    nsite = 80,
    n_true_event = c(21, 35, 70)
  )) %>% 
  mutate(total = window + gap) %>% 
  filter(window + gap <= n_true_event / 2)

scenario_grid <- cross_join(effort_scenarios, parameter_grid) %>% 
  mutate(scenario = paste0("W", row_number()))

# This data frame has one row per simulation condition:
all_iters_df <- left_join(scenario_grid,
                          as.data.frame(expand.grid(
                            scenario = unique(scenario_grid$scenario),
                            iter = 1:nsim
                          )))


#### Run the simulation in parallel using parabar ####

backend <- start_backend(ncores, cluster_type = "psock", backend_type = "async")

prep <- evaluate(backend, {
  library(tidyverse)
  library(nimble)
  library(unmarked)
  library(AICcmodavg)
  source("code_helper/helper_fn_clustered.R")
  source("code_helper/gof.R")
  NULL
})

prep <- parabar::export(backend, "all_iters_df")

# Execute the simulation
system.time(result <- par_lapply(backend, sample(1:nrow(all_iters_df),
                                                 size = nrow(all_iters_df)),
                                 fun = sim_and_evaluate_one_mod_clustered_par,
                                 datasummary_only = datasummary_only,
                                 design_df = all_iters_df, do_gof = FALSE, 
                                 fit_clustered_mod = FALSE, drop_intervals = TRUE))
stop_backend(backend)

# Process output
if (!datasummary_only) pred_perf_list <- bind_rows(lapply(result, function(x) x$predperf))
if (!datasummary_only) inference_list <- bind_rows(lapply(result, function(x) x$inference))
datasummary_list <- bind_rows(lapply(result, function(x) x$datasummary))
if (!datasummary_only) gof_list <- bind_rows(lapply(result, function(x) x$gof))
if (!datasummary_only) cdp_list <- bind_rows(lapply(result, function(x) x$CDP_df))


#### Save results ####

write_csv(bind_rows(datasummary_list), paste0("intermediate/windows_study_datasummary_results", tag, ".csv"))
write_csv(scenario_grid,               paste0("intermediate/windows_study_scenario_grid", tag, ".csv"))
if (!datasummary_only) write_csv(bind_rows(gof_list),       paste0("intermediate/windows_study_gof_results", tag, ".csv"))
if (!datasummary_only) write_csv(bind_rows(inference_list), paste0("intermediate/windows_study_inference_results", tag, ".csv"))
if (!datasummary_only) write_csv(bind_rows(pred_perf_list), paste0("intermediate/windows_study_predperf_results", tag, ".csv"))
if (!datasummary_only) write_csv(bind_rows(cdp_list),    paste0("intermediate/windows_study_CDP_results", tag, ".csv"))

