#########################################################################
#' Author: Benjamin R. Goldstein
#' vis.R
#' 
#' This file creates the plots in the manuscript.
#########################################################################

library(tidyverse)
library(gridExtra)

#### If necessary, reproduce summary files from raw output ####
source("code_helper/vis_fn.R")
source("code_helper/helper_fn_clustered.R")

if (!dir.exists("plots")) dir.create("plots")

#### Single-dataset example ####
# Set some parameters

logit_psi_mu <- logit(0.6)
overall_det <- 0.16
nsite <- 80
effort_days <- 21
det_breakdown <- 1

logit_p_mu <- logit(sqrt(overall_det / det_breakdown))

### Simulate data w no autocorrelation
autocorr_str <- 1
theta1 <- 1/(1/sqrt(overall_det*det_breakdown) + autocorr_str - 1)
theta2 <- theta1 * autocorr_str

set.seed(293528)

true_params_noautocorr <- c(logit_p_mu, -0.5, logit_psi_mu, 0.5, logit(theta1), logit(theta2))
dat_no_autocorr <- sim_occu_clustered(
  true_params = true_params_noautocorr, 
  nsite = nsite, 
  n_true_events = effort_days, 
  window = 1, 
  gap = 0
)


### Simulate data w STRONG autocorrelation
set.seed(293528)

autocorr_str <- 10
theta1 <- 1/(1/sqrt(overall_det*det_breakdown) + autocorr_str - 1)
theta2 <- theta1 * autocorr_str

true_params_autocorr <- c(logit_p_mu, -0.5, logit_psi_mu, 0.5, logit(theta1), logit(theta2))
dat_autocorr <- sim_occu_clustered(
  true_params = true_params_autocorr, 
  nsite = nsite, 
  n_true_events = effort_days, 
  window = 1, 
  gap = 0
)

#Plot the detection histories
results <- bind_rows(
  dat_no_autocorr$y_alt %>% dethist_mtx_to_long(tag = "No autocorr."),
  dat_autocorr$y_alt %>% dethist_mtx_to_long(tag = "Strong autocorr.")
)

dethist_panes <- results %>% 
  mutate(value = ifelse(value, "Detected", "Not detected")) %>% 
  bind_rows(data.frame(
    visit = -1, site = 1:nsite,
    value = ifelse(dat_no_autocorr$z, "Site occupied", "Site not occupied"),
    tag = "No autocorr."
  ),
  data.frame(
    visit = -1, site = 1:nsite,
    value = ifelse(dat_no_autocorr$z, "Site occupied", "Site not occupied"),
    tag = "Strong autocorr."
  )) %>% 
  ggplot() +
  geom_tile(aes(visit, site, fill = as.factor(value))) +
  facet_wrap(~tag) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual("", values = c("Not detected" = "#b2df8a",
                                   "Detected" = "#1f78b4",
                                   "Site occupied" = "#222222",
                                   "Site not occupied" = "#cccccc")) +
  ggtitle("A")


# Pct. occupied sites w/ detections
pct_sites_w_dets_plot <- data.frame(
  tag = c("No autocorr.", "Strong autocorr."),
  pct_sitewdet = c(sum(rowSums(dat_no_autocorr$y_alt) > 0) / sum(dat_no_autocorr$z),
                   sum(rowSums(dat_autocorr$y_alt) > 0) / sum(dat_autocorr$z))
) %>% 
  ggplot(aes(tag, pct_sitewdet)) +
  xlab("") + ylab("Realized CDP") +
  theme_minimal() +
  geom_col() +coord_flip() + ggtitle("C")
dets_per_site_w_dets_plot <- data.frame(
  tag = c("No autocorr.", "Strong autocorr."),
  pct_sitewdet = c(mean(rowSums(dat_no_autocorr$y_alt)[rowSums(dat_no_autocorr$y_alt) > 0]),
                   mean(rowSums(dat_autocorr$y_alt)[rowSums(dat_autocorr$y_alt) > 0]))
) %>% 
  ggplot(aes(tag, pct_sitewdet)) +
  xlab("") + ylab("Avg. detections per site with any") +
  theme_minimal() +
  geom_col() +coord_flip() + ggtitle("B")

layout_mtx <- matrix(c(
  1,1,
  1,1,
  1,1,
  2,3
), byrow = T, ncol = 2)

singleset_plot <- gridExtra::arrangeGrob(grobs = list(
  dethist_panes,dets_per_site_w_dets_plot, pct_sites_w_dets_plot
), layout_matrix = layout_mtx)
ggsave("plots/fig2_datexample.jpg", plot = singleset_plot, 
       width = 8, height = 6)

#### Lit review results ####

res <- read_csv("input_data/litreview_results.csv")

table(res$detection_window_days)

res$dur <- as.numeric(res$survey_duration_days)
res$window <- as.numeric(res$detection_window_days)

res$window[which(res$detection_window_days == "1, 5, 10")] <- 10
res$window[which(res$detection_window_days == "5, 11")] <- 11

p1 <- res %>% 
  filter(Qualified == "Yes") %>% 
  mutate(window = ifelse(window > 14, 14, window)) %>% 
  count(window) %>% 
  ggplot(aes(window, n)) +
  geom_col() +
  theme_minimal() +
  xlab("Detection window") + ylab("Num. papers") +
  scale_x_continuous(breaks = c(0:7*2, 14)) +
  ggtitle("A. Detection window") +
  theme(panel.grid.minor = element_blank())


p3 <- res %>% 
  filter(Qualified == "Yes") %>% 
  mutate(dur = ifelse(dur > 300, 300, dur)) %>% 
  ggplot(aes(dur)) +
  geom_histogram() +
  theme_minimal() +
  xlab("Avg. deployment duration") + ylab("Num. papers") +
  scale_x_continuous(breaks = 50*0:6) +
  ggtitle("C. Deployment duration") +
  theme(panel.grid.minor = element_blank())


p2 <- res %>% 
  filter(Qualified == "Yes") %>% 
  mutate(window_frac = window/dur) %>% 
  ggplot(aes(window_frac)) +
  geom_histogram() +
  theme_minimal() +
  xlab("Detection window / avg. deployment length") +
  ylab("Num. papers") +
  scale_x_continuous(breaks = 0:4/4) +
  ggtitle("B. Window vs. deployment length") +
  theme(panel.grid.minor = element_blank())


arrplot <- gridExtra::arrangeGrob(grobs = list(p1, p2, p3), nrow = 1)
ggsave("plots/fig3_litreview.jpg", arrplot, width = 10, height = 2.5)


#### Change in bias and RMSE ####
det_breakdown_ref <- 1
threshold <- 0.2#NULL

scenario_grid <- read_csv("intermediate/cluster_study_scenario_grid.csv")
all_res <- read_csv("intermediate/all_res_clustered.csv") %>% 
  mutate(effort_compare = ifelse(n_true_event <= 21, 21,
                                 ifelse(n_true_event <= 35, 35, 70)))

pp_res_comp <- all_res %>% 
  select(all_of(colnames(scenario_grid)), effort_compare, mean_RMSE_discrete) %>% 
  select(-scenario) %>% 
  filter(autocorr_str %in% c(1, 10)) %>% 
  group_by(nsite, effort_compare, logit_psi_mu, overall_det_prob, det_breakdown) %>% 
  mutate(autocorr_str = paste0("Autocorr", autocorr_str)) %>% 
  pivot_wider(names_from =  "autocorr_str", values_from = "mean_RMSE_discrete") %>% 
  mutate(autocorr_change_in_pp = Autocorr10 - Autocorr1,
         autocorr_str = "Autocorr 10 - 1")

drmse_plot <- pp_res_comp %>% 
  filter(window == 1, gap == 0) %>% 
  make_big_grid_plot("autocorr_change_in_pp", 
                     target_name = "dRMSE", 
                     min_max = c(min(pp_res_comp$autocorr_change_in_pp),
                                 max(pp_res_comp$autocorr_change_in_pp)), 
                     divergent = TRUE, ncol = 3,
                     title = "C. Change in RMSE due to autocorrelation",
                     discrete = FALSE)

drmse_plot_good <- drmse_plot

inf_res_comp_int <- all_res %>% 
  select(all_of(colnames(scenario_grid)), effort_compare, median_error_psi_mu) %>% 
  select(-scenario) %>% 
  filter(autocorr_str %in% c(1, 10)) %>%
  group_by(nsite, effort_compare, logit_psi_mu, overall_det_prob, det_breakdown) %>% 
  mutate(autocorr_str = paste0("Autocorr", autocorr_str)) %>% 
  pivot_wider(names_from =  "autocorr_str", values_from = "median_error_psi_mu") %>% 
  mutate(autocorr_change_in_error = Autocorr10 - Autocorr1,
         autocorr_error = Autocorr10,
         autocorr_str = "Autocorr 10 - 1")
inf_res_comp_int2 <- inf_res_comp_int %>% filter(window == 1, gap == 0, det_breakdown == det_breakdown_ref) 
dbias_plot_int <- inf_res_comp_int2 %>%
  make_big_grid_plot("autocorr_change_in_error",
                     target_name = "d error",
                     min_max = c(min(inf_res_comp_int2$autocorr_change_in_error),
                                 max(inf_res_comp_int2$autocorr_change_in_error)),                                  ncol = 3,
                     divergent = TRUE,
                     discrete = FALSE, use_threshold = FALSE, 
                     title = "A. Bias in logit-scale occupancy intercept")
dbias_plot_int <- inf_res_comp_int2 %>%
  make_big_grid_plot("autocorr_error",
                     target_name = "Error",
                     min_max = c(min(inf_res_comp_int2$autocorr_error),
                                 max(inf_res_comp_int2$autocorr_error)),                                  ncol = 3,
                     divergent = TRUE,
                     discrete = FALSE, use_threshold = FALSE, 
                     title = "A. Bias in logit-scale occupancy intercept")

inf_res_comp_beta <- all_res %>% 
  select(all_of(colnames(scenario_grid)), effort_compare, median_error_psi_beta) %>% 
  select(-scenario) %>% 
  filter(autocorr_str %in% c(1, 10)) %>%
  group_by(nsite, effort_compare, logit_psi_mu, overall_det_prob, det_breakdown) %>% 
  mutate(autocorr_str = paste0("Autocorr", autocorr_str)) %>% 
  pivot_wider(names_from =  "autocorr_str", values_from = "median_error_psi_beta") %>% 
  mutate(autocorr_change_in_error = Autocorr10 - Autocorr1,
         autocorr_error = Autocorr10,
         autocorr_str = "Autocorr 10 - 1")
inf_res_comp_beta2 <- inf_res_comp_beta %>% 
  filter(window == 1, gap == 0, det_breakdown == det_breakdown_ref)
dbias_plot_beta <- inf_res_comp_beta2 %>%
  make_big_grid_plot("autocorr_error",
                     target_name = "Error",
                     min_max = c(min(inf_res_comp_beta2$autocorr_error),
                                 max(inf_res_comp_beta2$autocorr_error)), 
                     ncol = 3,
                     divergent = TRUE,
                     discrete = FALSE, use_threshold = FALSE, 
                     title = "B. Bias in occupancy covariate effect")

bias_plot <- gridExtra::arrangeGrob(grobs = list(dbias_plot_int, dbias_plot_beta,
                                                 drmse_plot_good), ncol= 1)
ggsave("plots/fig4_autocorr_bias.jpg", bias_plot, width = 8, height = 12)

#### CDP vs. bias ####
all_res <- read_csv("intermediate/all_res_clustered.csv")
cdp_line_plot <- all_res %>% 
  filter(window == 1, gap == 0) %>% 
  # filter(det_breakdown == 1) %>%
  ggplot(aes(mean_CDP_error, median_error_psi_mu)) + 
  theme_minimal() + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(col = "darkblue") +
  geom_smooth(method = lm) +
  xlab("Error in estimated CDP (Window = 1)") +
  ylab("Bias in occu. intercept") + ggtitle("A. CDP vs. bias") +
  coord_fixed(1/3)

ggsave("plots/fig5_cdp.jpg", cdp_line_plot + theme(plot.title = element_blank()),
       width = 5, height = 5)


#### Windows plot ####
scenario_grid <- read_csv("intermediate/windows_study_scenario_grid.csv")
cdp_res       <- read_csv("intermediate/windows_study_CDP_results.csv")
inf_res       <- read_csv("intermediate/windows_study_inference_results.csv")
prd_res       <- read_csv("intermediate/windows_study_predperf_results.csv")

cdp_plots <- make_plots_bywindow(df = cdp_res, 
                                 outcome_column = "CDP_error", 
                                 outcome_column_name = "Error in est. CDP",
                                 title = "A. CDP error")

inf_int_plots <- make_plots_bywindow(inf_res[inf_res$par == "logit_psi_mu" &
                                               inf_res$modtype == "Discrete",],
                                     outcome_column = "abs_error", 
                                     outcome_column_name = "Error in est. of logit(psi)",
                                     title = "B. Occu. est. error")


SE_int_plots <- make_plots_bywindow(inf_res[inf_res$par == "logit_psi_mu" &
                                              inf_res$modtype == "Discrete",],
                                    outcome_column = "est_SE", hline = "none",
                                    outcome_column_name = "Standard error of est. of logit(psi)",
                                    title = "C. SE of occu. estimate")

coverage_int_plots <- make_plots_bywindow(inf_res[inf_res$par == "logit_psi_mu" &
                                                    inf_res$modtype == "Discrete",],
                                          summary_fn = mean,
                                          outcome_column = "covered", hline = 0.95,
                                          outcome_column_name = "95% CI coverage",
                                          title = "D. 95%CI coverage")

plotlist <- list(cdp_plots$dotplot, inf_int_plots$dotplot,
                 SE_int_plots$dotplot, coverage_int_plots$dotplot)
legend <- get_legend(plotlist[[1]])

plotlist <- lapply(plotlist, function(x) x + theme(legend.position = "none"))

four_plots <- gridExtra::arrangeGrob(
  grobs = plotlist,
  nrow = 2, right = legend
)

ggsave("plots/fig6_windowgap.jpg", four_plots, width = 9, height = 9)

#### GOF success rate ####

scenario_grid <- read_csv("intermediate/cluster_study_scenario_grid_gof.csv")
gof_res <- read_csv("intermediate/gof_res_clustered.csv") %>% 
  mutate(effort_compare = ifelse(n_true_event <= 21, 21,
                                 ifelse(n_true_event <= 35, 35, 70)))
gof_fig <- gof_res %>% 
  group_by(nsite, psi_mu = nimble::expit(logit_psi_mu), autocorr_str) %>% 
  summarize(reject_rate = mean(bad_fit_rate_adj)) %>% 
  mutate(autocorr_name = factor(paste0("Autocorr. str. = ", autocorr_str),
                                levels = paste0("Autocorr. str. = ", c(1,5,10)))) %>% 
  ggplot() +
  geom_col(aes(as.factor(nsite), reject_rate, 
               fill = as.factor(psi_mu), 
               group = psi_mu), position = "dodge",
           col = "black") +
  facet_wrap(~autocorr_name) +
  # theme_bw() +
  theme_minimal() +
  theme(strip.background = element_rect(color = "black", fill = "white")) +
  scale_fill_viridis_d("Mean occu.", end = 0.8) +
  xlab("Num. deployments") +
  ylab("Join count test detection rate (p < 0.05)")

ggsave("plots/fig7_gof.jpg", width = 7.5, height = 3)

#### Case study results ####
source("code_helper/case_study_fn.R")
source("code_helper/gof.R")
source("code_helper/helper_fn_clustered.R")
source("code_helper/vis_fn.R")

species_df <- read_csv("input_data/casestudy_target_specs.csv") %>% 
  left_join(hr_sizes, by = c("sciname" = "Species")) %>% 
  mutate(median_depl = NA)%>% 
  filter(!common_name %in% c("Desert Cottontail"))
gof_df <-         read_csv("intermediate/casestudy_gof_df.csv")
summary_df <-     read_csv("intermediate/casestudy_summary.csv")
param_space_df <- read_csv("intermediate/casestudy_param_space.csv")
alt_summary_df <- read_csv("intermediate/casestudy_alt_summary_df.csv")
alt_cdp_df <-     read_csv("intermediate/casestudy_alt_cdp_df.csv")
clustered_df <-   read_csv("intermediate/snapshot_clustered_results.csv")
lm_gof_dat <- gof_df %>% 
  left_join(species_df, by = c("species" = "common_name")) %>% 
  filter(grepl("Wright", test)) %>% 
  filter(test == "Wright Adj") %>% 
  mutate(pval_adj = p.adjust(pval, method = "fdr")) %>% 
  mutate(has_autocorr = pval < 0.05,
         has_autocorr_adj = pval_adj < 0.05)
panelA <- lm_gof_dat %>%
  mutate(species = gsub("_", " ", species)) %>% 
  ggplot(aes(species, pval_adj, shape = pval_adj < 0.05, 
             ymax = pval_adj, ymin = 0)) +
  geom_point(show.legend = F, size = 2.5) +
  ylim(c(0, 1)) +
  # facet_wrap(~test) + 
  coord_flip() +
  theme_bw() +
  geom_hline(yintercept = 0.05) +
  xlab("") + ylab("Join count test p-value") +
  ggtitle("A.") +
  scale_shape_manual(values = c(1, 19))

clustered_dat_toplot <- clustered_df %>% 
  filter(submodel == "occ", !is.na(submodel)) %>% 
  mutate(ymin = Estimate - 1.96*SE, ymax = Estimate + 1.96*SE) %>% 
  mutate(ymax = ifelse(ymax > 5, 5, ymax)) %>% 
  mutate(param = factor(gsub("_", " ", parname), 
                        levels = c("Intercept", "Forest", "Temp max")  )) %>% 
  mutate(window = "Clustered") %>% 
  filter(!species %in% c("California Ground Squirrel", "Snowshoe Hare"),
         !is.na(SE))


estimate_df_wider <- alt_summary_df %>%
  filter(submodel == "Occupancy", window %in% c(1, 10)) %>%
  mutate(window = paste0("Win", window, "_Est")) %>% 
  select(Estimate, param, submodel, species, window) %>% 
  pivot_wider(names_from = window, values_from = Estimate) %>% 
  mutate(param = gsub("_", " ", gsub("_scaled", "", param))) %>% 
  left_join(clustered_dat_toplot[clustered_dat_toplot$submodel == "occ",
                                 c("species", "param", "Estimate")]) %>% 
  rename(WinClustered_Est = Estimate)
res_df_wider <- alt_summary_df %>%
  filter(submodel == "Occupancy", window %in% c(1, 10)) %>%
  mutate(window = paste0("Win", window, "_SE")) %>% 
  select(SE, param, submodel, species, window) %>% 
  pivot_wider(names_from = window, values_from = SE) %>% 
  mutate(param = gsub("_", " ", gsub("_scaled", "", param))) %>% 
  left_join(clustered_dat_toplot[clustered_dat_toplot$submodel == "occ",
                                 c("species", "param", "SE")]) %>% 
  rename(WinClustered_SE = SE) %>% 
  left_join(estimate_df_wider) %>% 
  mutate(
    Win10_min = Win10_Est - 1.96*Win10_SE,
    Win10_max = Win10_Est + 1.96*Win10_SE,
    Win1_min = Win1_Est - 1.96*Win1_SE,
    Win1_max = Win1_Est + 1.96*Win1_SE
  )


compareplots <- list()
compareplots[[1]] <- res_df_wider %>% 
  filter(Win10_SE < 5) %>% 
  pivot_longer(cols = c("Win10_Est", "WinClustered_Est")) %>% 
  mutate(name = ifelse(name == "Win10_Est", "10-day window", "Clustered")) %>% 
  ggplot() +
  geom_boxplot(outlier.alpha = 0, aes(
    fill = name,
    factor(param, levels = c("Intercept", "Forest", "Temp max")),
    value - Win1_Est
  )) +
  theme_minimal() +
  coord_flip() +
  ylim(c(-0.25, 1)) +
  xlab("") +
  ylab("Change in param. estimate (vs. 1-day window)") +
  ggtitle("B.") +
  scale_fill_manual("", values = c("#1f78b4", "#b2df8a"))

all_legend <- get_legend(compareplots[[1]])
compareplots[[1]] <- compareplots[[1]] + theme(legend.position = "none")

compareplots[[2]] <- res_df_wider %>% 
  filter(Win10_SE < 5) %>% 
  pivot_longer(cols = c("Win10_SE", "WinClustered_SE")) %>% 
  mutate(name = ifelse(name == "Win10_SE", "10-day window", "Clustered")) %>% 
  ggplot() +
  geom_boxplot(outlier.alpha = 0, aes(
    fill = name, 
    factor(param, levels = c("Intercept", "Forest", "Temp max")),
    value - Win1_SE
  )) +
  theme_minimal() +
  coord_flip() +
  ylim(c(-0.025, 0.1)) +
  xlab("") +
  ylab("Change in param. SE (vs. 1-day window)") + 
  ggtitle("C.") +
  scale_fill_manual("", values = c("#1f78b4", "#b2df8a")) +
  theme(legend.position = "none")

plots <- arrangeGrob(grobs = compareplots, nrow = 2, right = all_legend)


casestudy_plots <- list(panelA, compareplots[[1]], compareplots[[2]])
layout_mtx <- matrix(c(
  1,2,2,
  1,3,3), byrow = TRUE, nrow = 2)
plots <- arrangeGrob(grobs = casestudy_plots, 
                     layout_matrix = layout_mtx,
                     right = all_legend, widths = c(1.8, 1, 1))

ggsave("plots/fig8_casestudy.jpg", plot = plots, 
       width = 9, height = 5)

