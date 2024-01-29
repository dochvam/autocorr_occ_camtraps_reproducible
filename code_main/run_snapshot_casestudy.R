#########################################################################
#' Author: Benjamin R. Goldstein
#' run_snapshot_casestudy.R
#' 
#' This file executes the analyses for the case study estimating occupancy
#' from Snapshot USA data.
#########################################################################

library(unmarked)
library(tidyverse)

# Load helper functions
source("code_helper/case_study_fn.R")
source("code_helper/gof.R")
source("code_helper/helper_fn_clustered.R")

# Load in the list of species that we want to iterate over
species_df <- read_csv("input_data/casestudy_target_specs.csv") %>% 
  mutate(median_depl = NA, num_cams = NA, num_cams_wdet = NA,
         num_dets = NA) %>% 
  filter(!common_name %in% c("Desert Cottontail"))


#### Read in the Snapshot detection histories and covariate data for all species ####
# These data are  formatted for use with unmarked. For more info see
# help("unmarkedFrameOccu")

umf_list <- list()
for (i in 1:nrow(species_df)) {
  
  common_name_clean <- gsub("'", "",
                            gsub("[ -]", "_", species_df$common_name[i]))
  
  umf_list[[i]] <- readRDS(paste0(
    "input_data/snapshot2020/datlist_", common_name_clean, "_1_0_2020.RDS"
  ))
  
  thisdat <- umf_list[[i]]@y
  
  species_df$median_depl[i]   <- median(rowSums(!is.na(thisdat), na.rm = T))
  species_df$num_cams[i]      <- nrow(thisdat)
  species_df$num_cams_wdet[i] <- sum(rowSums(thisdat, na.rm = T) > 0)
  species_df$num_dets[i]      <- sum(thisdat, na.rm = T)
}
write_csv(species_df, "intermediate/casestudy_species_summary.csv")


#### Fit the basic models ####
fit_list <- list()
summary_list <- list()
for (i in 1:nrow(species_df)) {

  fit_list[[i]] <- fit_om(umf = umf_list[[i]])
  
  summary_list[[i]] <- get_unmarked_summary(fit_list[[i]]) %>% 
    mutate(species = species_df$common_name[i],
           model = "Baseline",
           window = 1, gap = 0, duration = "Full")
}
summary_df <- bind_rows(summary_list)
convergence <- data.frame(
  conv = unlist(lapply(fit_list, function(x) x@opt$convergence)),
  species = species_df$common_name
)

write_csv(summary_df, "intermediate/casestudy_summary.csv")
write_csv(convergence, "intermediate/casestudy_convergence.csv")


#### Goodness-of-fit results ####

gof_df <- data.frame(
  test = "Wright Adj",
  species = species_df$common_name,
  pval = NA
)

for (i in 1:nrow(species_df)) {
  gof_df$pval[i] <- as.numeric(run_jc(fit_list[[i]], nsim = 500)$pval)
}

write_csv(gof_df, "intermediate/casestudy_gof_df.csv")



#### Alternative models: change detection window, gaps ####

# Choose what windows and gaps to iterate over
mods_to_fit <- data.frame(
  window = 1:30,
  gap = rep(0, 30)
)

alt_fit_list <- list()
alt_summary_list <- list()
alt_cdp_list <- list()

# Loop over species, fit each model
ct <- 0
pb <- progress::progress_bar$new(total = nrow(species_df) * nrow(mods_to_fit))
for (i in 1:nrow(species_df)) {
  chs <- umf_list[[i]]@obsCovs$Canopy_height_scaled[1 + (0:(nrow(umf_list[[i]]@y)-1))*ncol(umf_list[[i]]@y)]
  lrs <- umf_list[[i]]@obsCovs$log_roaddist_scaled[1 + (0:(nrow(umf_list[[i]]@y)-1))*ncol(umf_list[[i]]@y)]
  
  for (k in 1:nrow(mods_to_fit)) {
    ct <- ct + 1
    pb$tick()
    
    this_y_alt <- umf_list[[i]]@y %>% 
      y_to_y_alt(mods_to_fit$window[k], mods_to_fit$gap[k])
  
    this_umf_alt <-  unmarkedFrameOccu(this_y_alt, 
                                       siteCovs = cbind(
                                         umf_list[[i]]@siteCovs,
                                         Canopy_height_scaled = chs,
                                         log_roaddist_scaled = lrs
                                       ))
    
    alt_fit_list[[ct]] <- fit_om(this_umf_alt)
    
    alt_summary_list[[ct]] <- get_unmarked_summary(alt_fit_list[[ct]]) %>% 
      mutate(species = species_df$common_name[i],
             model = "Alt",
             window = mods_to_fit$window[k], 
             gap =  mods_to_fit$gap[k], duration = "Full")
  
    alt_cdp_list[[ct]] <- data.frame(
      window = mods_to_fit$window[k], 
      gap =  mods_to_fit$gap[k],
      species = species_df$common_name[i],
      nobs = sum(!is.na(this_umf_alt@y)),
      frac_sites_wreps = mean(rowSums(this_umf_alt@y, na.rm = T) > 1),
      ndays = sum(!is.na(this_umf_alt@y)) * mods_to_fit$window[k],
      mean_cdp = as.numeric(get_avg_cdp_unmarked(alt_fit_list[[ct]])[1])
    )
  }
}

alt_summary_df <- bind_rows(alt_summary_list)
alt_cdp_df <- bind_rows(alt_cdp_list)

write_csv(alt_summary_df, "intermediate/casestudy_alt_summary_df.csv")
write_csv(alt_cdp_df,     "intermediate/casestudy_alt_cdp_df.csv")

#### Run clustered model ####

umf_list <- list()
result_list <- list()
for (i in 1:nrow(species_df)) {
  writeLines(paste0("Modeling ", species_df$common_name[i]))
  umf_list[[i]] <- readRDS(paste0(
    "input_data/snapshot2020/datlist_", 
    gsub("[ -]", "_", gsub("'", "", species_df$common_name[i])), 
    "_1_0_2020.RDS"
  ))
  
  result_list[[i]] <- fit_hines_to_snapshot(umf = umf_list[[i]])
  result_list[[i]]$summary$species <- species_df$common_name[i]
}

result_df <- result_list %>% 
  lapply(function(x) x$summary) %>% 
  bind_rows()

write_csv(result_df, "intermediate/snapshot_clustered_results.csv")

#### Create summary table for paper ####

alt_df <- read_csv("intermediate/casestudy_alt_summary_df.csv") %>% 
  filter(window %in% c(1, 10)) %>% 
  mutate(modtype = paste0(window, "-day SOM")) %>% 
  select(common_name = species, modtype, submodel, param, Estimate, SE)
cluster_df <- read_csv("intermediate/snapshot_clustered_results.csv") %>% 
  mutate(modtype = "Clustered") %>% 
  select(common_name = species, modtype, submodel, param = parname, Estimate, SE)


output_table <- bind_rows(alt_df, cluster_df) %>% 
  mutate(submodel = recode(submodel, "occ" = "Occupancy", "det" = "Detection")) %>% 
  mutate(submodel = ifelse(grepl("Theta", param), "Site use", submodel)) %>% 
  mutate(param = recode(gsub("_scaled", "", param),
                           "Theta1" = "Theta",
                           "Theta2" = "Theta Prime",
                           "Temp_max" = "Max temp.", 
                           "Forest" = "Pct. forest cover", 
                           "log_roaddist" = "Log dist. to road", 
                           "Canopy_height" = "Canopy height")) %>% 
  left_join(species_df[, c("sciname", "common_name")]) %>% 
  select(Species = sciname, ModType = modtype, Submodel = submodel, 
         Param = param, Estimate, SE) %>% 
  arrange(Species, ModType)





