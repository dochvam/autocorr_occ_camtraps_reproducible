#########################################################################
#' run_simulation24.R
#' Author: Benjamin R. Goldstein
#' 
#' This file executes simulations 2 and 4 from the manuscript and produces output
#' files giving results from those simulations. Simulation 2 evaluates the 
#' performance of the join count test for identifying autocorrelated detections.
#' Simulation 4 estimates these datasets with the clustered occupancy model
#' to compare performance.
#########################################################################

library(tidyverse)
library(nimble)
library(unmarked)
library(AICcmodavg)
library(nimbleEcology)
library(parabar)

# Arbitrary seed used in the manuscript
set.seed(7233717)

# Read in functions for sim
source("code_helper/helper_fn_clustered.R")

nsim <- 1000
par_method <- "parabar"
ncores <- 44
tag <- "" # optional suffix to append to output filenames
datasummary_only <- FALSE


#### Simulation 1: Varying the strength of clustering by changing theta_2 ####

parameter_grid <- expand.grid(
  logit_psi_mu = logit(c(0.3, 0.4, 0.5, 0.6)),
  overall_det_prob = c(0.15, 0.2, 0.25, 0.3),
  det_breakdown = c(1/3, 1),
  autocorr_str = c(1, 2, 5, 10)
) %>%
  as.data.frame() %>% 
  filter(overall_det_prob * det_breakdown < 1 &
           overall_det_prob / det_breakdown < 1)


effort_scenarios <- data.frame(
  window = c(1),
  gap = c(0)
) %>% 
  cross_join(expand.grid(  # Effort variables
    nsite = c(40, 80, 120),
    n_true_event = c(21, 35, 70)
  )) %>% 
  mutate(total = window + gap) %>% 
  filter(total <= n_true_event / 2)


scenario_grid <- cross_join(effort_scenarios, parameter_grid) %>% 
  mutate(scenario = row_number())

all_iters_df <- left_join(scenario_grid,
                          as.data.frame(expand.grid(
                            scenario = 1:max(scenario_grid),
                            iter = 1:nsim
                          )))
#### Start parallel simulations ####

  
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

# Execute parallel sims
system.time(result <- par_lapply(backend, sample(1:nrow(all_iters_df),
                                                 # size = 1000),
                                                 size = nrow(all_iters_df)),
                                 fun = sim_and_evaluate_one_mod_clustered_par,
                                 datasummary_only = datasummary_only,
                                 design_df = all_iters_df, do_gof = TRUE, 
                                 fit_clustered_mod = TRUE))
stop_backend(backend)

# Process output
if (!datasummary_only) pred_perf_list <- bind_rows(lapply(result, function(x) x$predperf))
if (!datasummary_only) inference_list <- bind_rows(lapply(result, function(x) x$inference))
datasummary_list <- bind_rows(lapply(result, function(x) x$datasummary))
if (!datasummary_only) gof_list <- bind_rows(lapply(result, function(x) x$gof))
if (!datasummary_only) cdp_list <- bind_rows(lapply(result, function(x) x$CDP_df))

#### Save results ####
write_csv(bind_rows(datasummary_list), paste0("intermediate/cluster_study_datasummary_results_gof", tag, ".csv"))
write_csv(scenario_grid,               paste0("intermediate/cluster_study_scenario_grid_gof", tag, ".csv"))


if (!datasummary_only) write_csv(bind_rows(gof_list),       paste0("intermediate/cluster_study_gof_results_gof", tag, ".csv"))
if (!datasummary_only) write_csv(bind_rows(inference_list), paste0("intermediate/cluster_study_inference_results_gof", tag, ".csv"))
if (!datasummary_only) write_csv(bind_rows(pred_perf_list), paste0("intermediate/cluster_study_predperf_results_gof", tag, ".csv"))
if (!datasummary_only) write_csv(bind_rows(cdp_list),    paste0("intermediate/cluster_study_CDP_results_gof", tag, ".csv"))

