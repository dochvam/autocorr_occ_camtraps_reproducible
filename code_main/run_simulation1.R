#########################################################################
#' run_simulation1.R
#' Author: Benjamin R. Goldstein
#' 
#' This file executes simulation 1 from the manuscript and produces output
#' files giving results from those simulations. Simulation 1 looks at how
#' autocorrelation causes biased estimation and incorrect uncertainty 
#' quantification depending on sampling design (number of sites, duration of
#' surveys) and ecological conditions (mean occupancy, mean detection,
#' strength of autocorrelation)
#########################################################################


library(tidyverse)
library(nimble)
library(unmarked)
library(AICcmodavg)
library(nimbleEcology)
library(parabar)

# Arbitrary seed used for the manuscript
set.seed(0215296)

# Read in functions for simulating data and estimating models
source("code_helper/helper_fn_clustered.R")

# Simulation parameters
nsim <- 1#000 # Number of datasets per condition
par_method <- "parabar"
ncores <- 44
tag <- "" # optional suffix to append to output filenames
datasummary_only <- FALSE


#### Set up simulation conditions grid ####

parameter_grid <- expand.grid(
  logit_psi_mu = logit(c(0.3, 0.4, 0.5, 0.6)),
  overall_det_prob = c(0.15, 0.2, 0.25, 0.3),
  det_breakdown = c(1/3, 1),
  autocorr_str = c(1, 5, 10)
) %>%
  as.data.frame() %>% 
  filter(overall_det_prob * det_breakdown < 1 &
           overall_det_prob / det_breakdown < 1)


effort_scenarios <- data.frame(
  window = 1, gap = 0
  ) %>% 
  cross_join(expand.grid(  # Effort variables
    nsite = c(40, 80, 120),
    n_true_event = c(21, 35, 70)
  )) %>% 
  mutate(total = window + gap) %>% 
  filter(total <= n_true_event / 2)# %>% 
  # mutate(n_true_event = total * floor(n_true_event / total))

scenario_grid <- cross_join(effort_scenarios, parameter_grid) %>% 
  mutate(scenario = row_number())

# This data frame has one row per simulated dataset
all_iters_df <- left_join(scenario_grid,
                          as.data.frame(expand.grid(
                            scenario = 1:max(scenario_grid),
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

# Execute the simulations
system.time(result <- par_lapply(backend, sample(1:nrow(all_iters_df),
                                                 size = nrow(all_iters_df)),
                    fun = sim_and_evaluate_one_mod_clustered_par,
                    datasummary_only = datasummary_only,
                    design_df = all_iters_df, do_gof = FALSE, 
                    fit_clustered_mod = FALSE))
stop_backend(backend)

# Process output
if (!datasummary_only) pred_perf_list <- bind_rows(lapply(result, function(x) x$predperf))
if (!datasummary_only) inference_list <- bind_rows(lapply(result, function(x) x$inference))
datasummary_list <- bind_rows(lapply(result, function(x) x$datasummary))
if (!datasummary_only) gof_list <- bind_rows(lapply(result, function(x) x$gof))
if (!datasummary_only) cdp_list <- bind_rows(lapply(result, function(x) x$CDP_df))

#### Save results ####
  
write_csv(bind_rows(datasummary_list), paste0("intermediate/cluster_study_datasummary_results", tag, ".csv"))
write_csv(scenario_grid,               paste0("intermediate/cluster_study_scenario_grid", tag, ".csv"))

if (!datasummary_only) write_csv(bind_rows(gof_list),       paste0("intermediate/cluster_study_gof_results", tag, ".csv"))
if (!datasummary_only) write_csv(bind_rows(inference_list), paste0("intermediate/cluster_study_inference_results", tag, ".csv"))
if (!datasummary_only) write_csv(bind_rows(pred_perf_list), paste0("intermediate/cluster_study_predperf_results", tag, ".csv"))
if (!datasummary_only) write_csv(bind_rows(cdp_list),    paste0("intermediate/cluster_study_CDP_results", tag, ".csv"))

