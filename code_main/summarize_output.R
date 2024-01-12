#########################################################################
#' Author: Benjamin R. Goldstein
#' summarize_output.R
#' 
#' This file processes and summarizes the output data from simulations 
#' 1, 3, and 4 for use in visualizations. The goal is to end up with a 
#' single result for each observed metric for each simulation condition.
#' If you used a "tag" in your output, you may need to adjust the filenames.
#########################################################################


library(tidyverse)

#### Summarize output from simulation 1 ####

# Read in all the results
scenario_grid <- read_csv("intermediate/cluster_study_scenario_grid.csv")
inf_raw <- read_csv("intermediate/cluster_study_inference_results.csv")
datasummary_raw <- read_csv("intermediate/cluster_study_datasummary_results.csv")
cdp_raw <- read_csv("intermediate/cluster_study_CDP_results.csv")
pp_raw <- read_csv("intermediate/cluster_study_predperf_results.csv")
inds_to_drop <- inf_raw %>% 
  filter(is.na(abs_error) & modtype == "Discrete") %>% 
  distinct(iter, scenario)

# Summarize predictive performance results
pp_res <- pp_raw %>%  
  filter(RMSE_type == "Out-of-sample") %>% 
  select(-type) %>% 
  pivot_wider(names_from = modtype, values_from = RMSE) %>% 
  mutate(RMSE_diff = Discrete - Clustered) %>% 
  group_by(scenario) %>% 
  summarize(mean_RMSE_diff = mean(RMSE_diff),
            mean_RMSE_discrete = mean(Discrete),
            RMSE_sd = sd(Discrete))

# Summarize inference results related to bias
inf_res_bias <- inf_raw %>% 
  filter(modtype == "Discrete") %>% 
  select(-type, -true_val, -est_SE, -covered, -est) %>% 
  pivot_wider(names_from = par, values_from = abs_error) %>% 
  group_by(scenario) %>% 
  summarize(median_error_psi_mu   = median(logit_psi_mu),
            median_error_psi_beta = median(psi_beta1))

# Summarize inference results related to standard errors
inf_res_SE <- inf_raw %>% 
  filter(!paste0(iter, "_", scenario) %in%
           paste0(inds_to_drop$iter, "_", inds_to_drop$scenario)) %>% 
  filter(modtype == "Discrete") %>% 
  select(-type, -est, -true_val, -abs_error, -covered) %>% 
  pivot_wider(names_from = par, values_from = est_SE) %>% 
  group_by(scenario) %>% 
  summarize(median_SE_psi_mu   = median(logit_psi_mu, na.rm = TRUE),
            median_SE_psi_beta = median(psi_beta1, na.rm = TRUE))

# Summarize inference results related to the coverage of the 95% CI.
inf_res_coverage <- inf_raw %>% 
  filter(!paste0(iter, "_", scenario) %in%
           paste0(inds_to_drop$iter, "_", inds_to_drop$scenario)) %>% 
  filter(modtype == "Discrete") %>% 
  select(-type, -est, -true_val, -abs_error, -est_SE) %>% 
  pivot_wider(names_from = par, values_from = covered) %>% 
  group_by(scenario) %>% 
  summarize(coverage_95_psi_mu   = mean(logit_psi_mu, na.rm = TRUE),
            coverage_95_psi_beta = mean(psi_beta1, na.rm = TRUE))

# Summarize characteristics of the simulated data
data_res <- datasummary_raw %>% 
  group_by(scenario, z) %>% 
  filter(z == 1) %>% 
  summarize(nsite_occ = sum(n), nsite_occ_w_det = sum(n * as.numeric(ndet > 0))) %>% 
  mutate(realized_CDP = nsite_occ_w_det / nsite_occ)

# Summarize results related to estimating the cumulative detection prob.
cdp_res <- cdp_raw %>% 
  group_by(scenario) %>% 
  summarize(mean_true_p = mean(mean_true_p),
            mean_p_est = mean(mean_p_est),
            mean_CDP_est = mean(mean_CDP_est),
            nobs_per_site = mean(nobs_per_site),
            mean_CDP_error = mean(CDP_error)) %>% 
  left_join(data_res[, c("scenario", "realized_CDP")]) %>% 
  mutate(CDP_pct_error = mean_CDP_error/realized_CDP, by = "scenario") 


# Join it all together
all_res <- pp_res %>% 
  left_join(inf_res_bias, by = "scenario") %>% 
  left_join(inf_res_SE, by = "scenario") %>% 
  left_join(inf_res_coverage, by = "scenario") %>% 
  left_join(data_res, by = "scenario") %>% 
  left_join(cdp_res, by = "scenario") %>% 
  left_join(scenario_grid, by = "scenario") 

# Calculate real-scale error
all_res$median_error_psi_mu_realscale <-
  nimble::expit(all_res$logit_psi_mu + all_res$median_error_psi_mu) -
  nimble::expit(all_res$logit_psi_mu)

# Save the results.
# Note that the exact results from the manuscript as produced here are provided
# in the directory "results_tables"
write_csv(all_res, "intermediate/all_res_clustered.csv")


#### Summarize output from simulation 3 and 4 ####

# Read everything in
scenario_grid_gof <- read_csv("intermediate/cluster_study_scenario_grid_gof.csv")
inf_raw <- read_csv("intermediate/cluster_study_inference_results_gof.csv")
pp_raw <- read_csv("intermediate/cluster_study_predperf_results_gof.csv")
gof_raw <- read_csv("intermediate/cluster_study_gof_results_gof.csv")
inds_to_drop <- inf_raw %>% 
  filter(is.na(abs_error) & modtype == "Discrete") %>% 
  distinct(iter, scenario)

# Summarize predictive performance from clustered model
pp_res <- pp_raw %>%  
  filter(RMSE_type == "Out-of-sample") %>% 
  select(-type) %>% 
  pivot_wider(names_from = modtype, values_from = RMSE) %>% 
  mutate(RMSE_diff = Discrete - Clustered) %>% 
  group_by(scenario) %>% 
  summarize(mean_RMSE_diff = mean(RMSE_diff),
            mean_RMSE_discrete = mean(Discrete),
            mean_RMSE_clustered = mean(Clustered),
            RMSE_sd = sd(Discrete))

# Summarize inference from clustered model
inf_res_bias <- inf_raw %>% 
  filter(modtype == "Clustered") %>% 
  select(-type, -true_val, -est_SE, -covered, -est) %>% 
  pivot_wider(names_from = par, values_from = abs_error) %>% 
  group_by(scenario) %>% 
  summarize(median_error_psi_mu   = median(logit_psi_mu),
            median_error_psi_beta = median(psi_beta1))

inf_res_SE <- inf_raw %>% 
  filter(!paste0(iter, "_", scenario) %in%
           paste0(inds_to_drop$iter, "_", inds_to_drop$scenario)) %>% 
  filter(modtype == "Clustered") %>% 
  select(-type, -est, -true_val, -abs_error, -covered) %>% 
  pivot_wider(names_from = par, values_from = est_SE) %>% 
  group_by(scenario) %>% 
  summarize(median_SE_psi_mu   = median(logit_psi_mu, na.rm = TRUE),
            median_SE_psi_beta = median(psi_beta1, na.rm = TRUE))


inf_res_coverage <- inf_raw %>% 
  filter(!paste0(iter, "_", scenario) %in%
           paste0(inds_to_drop$iter, "_", inds_to_drop$scenario)) %>% 
  filter(modtype == "Clustered") %>% 
  select(-type, -est, -true_val, -abs_error, -est_SE) %>% 
  pivot_wider(names_from = par, values_from = covered) %>% 
  group_by(scenario) %>% 
  summarize(coverage_95_psi_mu   = mean(logit_psi_mu, na.rm = TRUE),
            coverage_95_psi_beta = mean(psi_beta1, na.rm = TRUE))

# Summarize goodness-of-fit results
gof_res_Adj <- gof_raw %>% 
  select(scenario, iter, test, pval) %>% 
  filter(test == "Wright Adj.") %>% 
  group_by(scenario) %>% 
  summarize(bad_fit_rate_adj = mean(pval < 0.05),
            median_pval_adj = median(pval),
            Q90_pval_adj = quantile(pval, 0.9))

gof_res_All <- gof_raw %>% 
  select(scenario, iter, test, pval) %>% 
  filter(test == "Wright All") %>% 
  group_by(scenario) %>% 
  summarize(bad_fit_rate_all = mean(pval < 0.05),
            median_pval_all = median(pval),
            Q90_pval_all = quantile(pval, 0.9))


gof_res <- pp_res %>% 
  left_join(inf_res_bias, by = "scenario") %>% 
  left_join(inf_res_SE, by = "scenario") %>% 
  left_join(inf_res_coverage, by = "scenario") %>% 
  left_join(gof_res_All, by = "scenario") %>% 
  left_join(gof_res_Adj, by = "scenario") %>% 
  left_join(scenario_grid_gof, by = "scenario") 


# Note that the exact results from the manuscript as produced here are provided
# in the directory "results_tables"
write_csv(gof_res, "intermediate/gof_res_clustered.csv")

