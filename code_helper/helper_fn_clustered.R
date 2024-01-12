#########################################################################
#' Author: Benjamin R. Goldstein
#' helper_fn_clustered.R
#' 
#' This file holds helper functions for use in simulating data and 
#' estimating models. These functions are required throughout the study
#' pipeline.
#########################################################################

library(nimble)
library(unmarked)
library(AICcmodavg)
library(nimbleEcology)

source('code_helper/gof.R')

# Define the order of parameters in the clustered model
true_parnames <- c("logit_p_mu", "p_beta2", 
                   "logit_psi_mu", "psi_beta1",
                   "logit_theta", "logit_theta_prime")

# Default values
true_params_baseline <- c(logit(0.8), -0.5, logit(0.5), 
                          0.5, logit(0.1), logit(0.5))

#' @name y_to_y_alt
#' @param y A dataset with observations according to a 1-day detection window
#' @param interval The detection window length PLUS gap length
#' @param gap The detection gap length
#' @description Aggregates a detection history according to the specified
#'              detection window and gap.
y_to_y_alt <- function(y, interval, gap) {
  
  y_alt <- matrix(0, nrow = nrow(y), ncol = ncol(y) / interval)
  
  if (interval == 1) {
    y_alt <- y
  } else {
    for (i in 1:nrow(y_alt)) {
      for (j in 1:ncol(y_alt)) {
        inds <- (1 + (j-1) * interval):(interval * j - gap)
      
        if (all(!is.na(y[i, inds]))) {
          y_alt[i,j] <- as.numeric(any(y[i, inds] == 1))
        } else {
          y_alt[i,j] <- NA
        }
      }
    }
  }
  y_alt
}

#' @name sim_occu_clustered
#' @param true_params The length-6 vector of parameter values. (See the
#'                    definition for "true_param_names" above)
#' @param nsite Number of closed camera deployment sites
#' @param n_true_events Number of days over which each camera is sampling
#' @param window Detection window to aggregate to
#' @param gap Detection gap to aggregate to
#' @param drop_intervals If (window + gap) % n_true_events > 0, should leftover
#'                       samples be discarded (TRUE), or should deployment lengths 
#'                       be adjusted such that n_true_events is the average
#'                       deployment length (FALSE)
#' @description Simulate data according to the clustered occupancy model.
sim_occu_clustered <- function(true_params,
                               nsite, n_true_events,
                               window, gap, drop_intervals) {
  logit_p_mu        <- true_params[1]
  p_beta2           <- true_params[2]
  logit_psi_mu      <- true_params[3]
  psi_beta1         <- true_params[4]
  theta             <- expit(true_params[5])
  theta_prime       <- expit(true_params[6])
  
  if (n_true_events <= 1) {
    stop("n_true_events must be more than 1")
  }
  
  interval <- window + gap
  
  if (n_true_events %% interval == 0) {
    nrep_vec <- rep(n_true_events, nsite)
    nobs_vec <- nrep_vec / interval
    
  } else if (!drop_intervals) {
    
    n_true_event_lower <- floor(n_true_events/interval) * interval
    n_true_event_upper <- ceiling(n_true_events/interval) * interval
    
    ewlo <- (n_true_events-n_true_event_upper)/(n_true_event_lower-n_true_event_upper)
    ewhi <- 1 - ewlo
    
    n_lo <- ceiling(ewlo * nsite - 0.0001) # Correct rounding errors
    n_hi <- floor(ewhi * nsite + 0.0001)
      
    nrep_vec <- c(rep(n_true_event_lower, n_lo),
                  rep(n_true_event_upper, n_hi))
    nobs_vec <- nrep_vec / interval
  } else {
    n_true_event_lower <- floor(n_true_events/interval) * interval
    
    nrep_vec <- rep(n_true_event_lower, nsite)
    nobs_vec <- nrep_vec / interval
  }

  # Simulate covariate data
  cov1 <- rnorm(nsite, 0, 1)
  cov2 <- rnorm(nsite, 0, 1)
  
  psi  <- expit(logit_psi_mu +  cov1 * psi_beta1)
  z    <- rbinom(nsite, size = 1, prob = psi)
  
  p    <- expit(logit_p_mu +  cov2 * p_beta2)
  
  theta_equil <- theta / (theta + 1 - theta_prime)
  
  occurred <- y <- matrix(NA, nrow = nsite, ncol = max(nrep_vec))
  
  for (i in 1:nsite) {
    if (z[i]) {
      occurred[i, 1] <- rbinom(1, 1, theta_equil) # Follow up: I think this should be the equilibrium of thetas, not theta
      if (nrep_vec[i] > 1) {
        for (j in 2:nrep_vec[i]) {
          occurred[i, j] <- rbinom(1, 1, theta * (1-occurred[i, j-1]) + theta_prime * occurred[i, j-1])
        }
      }
      
      y[i, 1:nrep_vec[i]] <- rbinom(nrep_vec[i], 1, p[i] * occurred[i, ])
    } else {
      occurred[i, 1:nrep_vec[i]] <- 0
      y[i, 1:nrep_vec[i]] <- 0
    }
  }
  
  y_alt <- y_to_y_alt(y, interval, gap)
  
  return(list(
    y_alt = y_alt,
    y_base = y, 
    occurred = occurred,
    true_params = true_params,
    true_psi = psi, true_p = p, 
    cov1 = cov1, cov2 = cov2, z = z,
    nrep_vec = nrep_vec, 
    nobs_vec = nobs_vec
  ))
}

#' @name sim_and_evaluate_one_mod_clustered
#' @param iter The iteration of the survey condition
#' @param design_df The data frame that defines the simulation study design
#' @param row The row of "design_df" that we're on
#' @param do_gof If TRUE, the join count test is conducted on the SOM; otherwise
#'              it is not
#' @param fit_clustered_mod If TRUE, a clustered model is fit to the data;
#'              otherwise, the clustered model is skipped
#' @param datasummary_only If TRUE, data are simulated and summarized, and then
#'              the function ends without estimating any models
#' @param drop_intervals Passed to "sim_occu_clustered"
#' @description Simulate data according to the clustered occupancy model.
sim_and_evaluate_one_mod_clustered <- function(
    iter, design_df, row, do_gof = TRUE, fit_clustered_mod = TRUE,
    datasummary_only = FALSE, drop_intervals = FALSE
  ) {
  
  overall_det   <- design_df$overall_det_prob[row]
  det_breakdown <- design_df$det_breakdown[row]
  autocorr_str  <- design_df$autocorr_str[row]
  
  theta1 <- 1/(1/sqrt(overall_det*det_breakdown) + autocorr_str - 1)
  theta2 <- theta1 * autocorr_str
  
  true_params <- true_params_baseline
  true_params[3] <- design_df$logit_psi_mu[row]
  
  true_params[1] <- logit(sqrt(overall_det / det_breakdown))
  
  true_params[5] <- logit(theta1)
  true_params[6] <- logit(theta2)
  
  window <- design_df$window[row]
  gap    <- design_df$gap[row]
  
  
  # Simulate data
  sim_dat <- sim_occu_clustered(
    true_params, 
    nsite = design_df$nsite[row], 
    n_true_events = design_df$n_true_event[row],
    window = window,
    gap = gap,
    drop_intervals = drop_intervals
  )
  
  
  dat_summary <- data.frame(z = sim_dat$z, 
                            ndet = rowSums(sim_dat$y_alt, na.rm = TRUE),
                            nobs = rowSums(!is.na(sim_dat$y_alt))) %>% 
    count(z, ndet) %>% 
    mutate(iter = iter, scenario = design_df$scenario[row])
  
  
  if (datasummary_only) return(list(
    datasummary = dat_summary
  ))
  
  
  sim_dat_oos <- sim_occu_clustered(
    true_params, 
    nsite = design_df$nsite[row], 
    n_true_events = design_df$n_true_event[row],
    window = window,
    gap = gap, 
    drop_intervals = drop_intervals
  )
  
  # Fit each of the discrete occupancy models
  discrete_inference_list <- list()
  discrete_predperf_list <- list()
  discrete_gof_list <- list()
  
  for (f in 2:2) {
    if (f == 1) {
      umf <- unmarked::unmarkedFrameOccu(
        y = sim_dat$y_base,
        siteCovs = data.frame(
          cov1 = sim_dat$cov1, 
          cov2 = sim_dat$cov2
        )
      )
    } else {
      umf <- unmarked::unmarkedFrameOccu(
        y = sim_dat$y_alt,
        siteCovs = data.frame(
          cov1 = sim_dat$cov1, 
          cov2 = sim_dat$cov2
        )
      )
    }
    

    this_fit_discrete <- occu(~cov2 ~cov1, dat = umf, 
                              control = list(maxit = 1000))
    captured <- capture.output(s <- summary(this_fit_discrete))

    est_vals <- c(s$det$Estimate, s$state$Estimate)
    est_SEs <- c(s$det$SE, s$state$SE)
    
    discrete_inference_list[[f]] <- data.frame(
      par = c("logit_p_mu", "p_beta2", "logit_psi_mu", "psi_beta1"),
      true_val = true_params[1:4],
      est = est_vals,
      abs_error = est_vals - true_params[1:4],
      est_SE = est_SEs,
      covered = 
        est_vals - est_SEs*1.96 < true_params[3:4] &
        est_vals + est_SEs*1.96 > true_params[3:4],
      modtype = "Discrete", type = ifelse(f == 1, "Baseline", "Window/Gap"),
      iter = iter, scenario = design_df$scenario[row]
    )
    
    if (do_gof) {
      # captured <- capture.output(
      #   gof_mb <- mb.gof.test(this_fit_discrete, nsim = 100,
      #                         parallel = FALSE)
      # )
      gof_result <- as.numeric(run_jc(this_fit_discrete, nsim = 500)$pval)
      
      discrete_gof_list[[f]] <- data.frame(
        test = "Wright Adj.",
        pval = gof_result,
        modtype = "Discrete",
        type = ifelse(f == 1, "Baseline", "Window/Gap"),
        iter = iter, scenario = design_df$scenario[row]
      )
      # captured <- capture.output(
      #   gof_parboot <- parboot(this_fit_discrete, fitstats, 
      #                          nsim = 100, parallel = FALSE)
      # )
      # 
      # discrete_gof_list[[f]] <- data.frame(
      #   test = c("SSE", "Wright Adj.", "Wright All" #"MB"
      #            ),
      #   pval = c(mean(gof_parboot@t.star[,1] >= gof_parboot@t0[1]),
      #            mean(gof_parboot@t.star[,2] >= gof_parboot@t0[2]),
      #            mean(gof_parboot@t.star[,3] >= gof_parboot@t0[3])#,
      #            #gof_mb$p.value
      #            ),
      #   modtype = "Discrete",
      #   type = ifelse(f == 1, "Baseline", "Window/Gap"), 
      #   iter = iter, scenario = design_df$scenario[row]
      # )
    } else {
      discrete_gof_list[[f]] <- data.frame(
        dummy=NA
      )
    }
    
    rmse_in <- get_rmse_unmarked(this_fit = this_fit_discrete, dat = sim_dat, which_y = f)
    rmse_out <- get_rmse_unmarked(this_fit = this_fit_discrete, dat = sim_dat_oos, which_y = f)
    discrete_predperf_list[[f]] <- data.frame(
      RMSE_type = c("In-sample", "Out-of-sample"),
      RMSE = c(rmse_in, rmse_out),
      modtype = "Discrete", 
      type = ifelse(f == 1, "Baseline", "Window/Gap"), 
      iter = iter, scenario = design_df$scenario[row]
    )
  }
  
  if (fit_clustered_mod) {
    fit_clustered <- tryCatch(
      fit_cluster_occ_mod(dat = sim_dat),
      error = function(e) NA
    )
  } else {
    fit_clustered <- NA
  }
  if (is.list(fit_clustered)) {
    
    discrete_inference_list[[3]] <- data.frame(
      par = true_parnames,
      true_val = true_params,
      est = fit_clustered$summary$Estimate,
      abs_error = fit_clustered$summary$Estimate - true_params,
      est_SE = fit_clustered$summary$SE
    ) %>% 
      mutate(
        covered = 
          fit_clustered$summary$Estimate - est_SE*1.96 < true_val &
          fit_clustered$summary$Estimate + est_SE*1.96 > true_val,
        modtype = "Clustered", type = "Clustered",
        iter = iter, scenario = design_df$scenario[row]
      )
    
    
    rmse_in <- get_rmse_clustered(this_fit = fit_clustered, dat = sim_dat)
    rmse_out <- get_rmse_clustered(this_fit = fit_clustered, dat = sim_dat_oos)
    
    
    
    discrete_predperf_list[[3]] <- data.frame(
      RMSE_type = c("In-sample", "Out-of-sample"),
      RMSE = c(rmse_in, rmse_out),
      modtype = "Clustered", 
      type = "Clustered", 
      iter = iter, scenario = design_df$scenario[row]
    ) 
  } else {
    discrete_inference_list[[3]] <- data.frame(
      par = true_parnames,
      true_val = true_params,
      abs_error = NA,
      est_SE = NA,
      covered = NA,
      modtype = "Clustered", type = "Clustered",
      iter = iter, scenario = design_df$scenario[row]
    )
    
    discrete_predperf_list[[3]] <- data.frame(
      RMSE_type = c("In-sample", "Out-of-sample"),
      RMSE = c(NA, NA),
      modtype = "Clustered", 
      type = "Clustered", 
      iter = iter, scenario = design_df$scenario[row]
    ) 
  }
  
  est_logit_p_mu <- discrete_inference_list[[2]]$est[1] 
  est_p_beta2 <- discrete_inference_list[[2]]$est[2]
  
  CDP_df_bysite <- data.frame(
    iter = iter, scenario = design_df$scenario[row],
    z = sim_dat$z, 
    nobs = rowSums(!is.na(sim_dat$y_alt)),
    true_p = sim_dat$true_p,
    p_est = expit(est_logit_p_mu + sim_dat$cov2 * est_p_beta2)
  ) %>% 
    mutate(CDP_est = 1 - (1 - p_est)^nobs)
  
  CDP_df = data.frame(
    iter = iter, scenario = design_df$scenario[row],
    mean_true_p = mean(CDP_df_bysite$true_p),
    mean_p_est = mean(CDP_df_bysite$p_est),
    mean_CDP_est = mean(CDP_df_bysite$CDP_est),
    true_CDP = sum(rowSums(sim_dat$y_alt, na.rm = TRUE) > 0) / sum(sim_dat$z),
    nobs_per_site = mean(CDP_df_bysite$nobs)
  ) %>% 
    mutate(CDP_error = mean_CDP_est - true_CDP)
  
  
  return(list(
    inference = bind_rows(discrete_inference_list),
    predperf = bind_rows(discrete_predperf_list),
    datasummary = dat_summary,
    gof = bind_rows(discrete_gof_list),
    CDP_df = CDP_df
  ))
}

nimbleOptions(enableDerivs = FALSE)
dHMM_compiled <- compileNimble(dHMM)


#' @name cluster_occ_ll
#' @param par A vector of parameters for the clustered occupancy model, of 
#'            length 6. Its elements are, in order: the logit-scale mean 
#'            detection probability; the effect of the detection covariate;
#'            the logit-scale mean occupancy probability; the effect of the
#'            occupancy covariate; logit(theta); and logit(theta').
#' @param dat A named list of covariate and detection data 
#' @param site_has_no_obs A vector of 1s and 0s. If element i=1, it means that 
#'          site i had no observations of the species. Used for faster 
#'          computation.
#' @description This function takes a parameter vector and data and calculates
#'              the log-likelihood of the data given the parameters and the
#'              clustered occupancy model. This is the objective function that 
#'              is optimized when estimating the clustered model for the simulation 
#'              study. It uses dHMM from nimbleEcology.
cluster_occ_ll <- function(par, dat, site_has_no_obs) {
  cov1 <- dat$cov1
  cov2 <- dat$cov2
  y <- dat$y_base
  
  logit_p_mu        <- par[1]
  p_beta2           <- par[2]
  logit_psi_mu      <- par[3]
  psi_beta1         <- par[4]
  theta             <- expit(par[5])
  theta_prime       <- expit(par[6])
  
  psi  <- expit(logit_psi_mu +  cov1 * psi_beta1)
  p    <- expit(logit_p_mu +  cov2 * p_beta2)  
  theta_equil <- theta / (theta + 1 - theta_prime)
  
  ll_per_site <- numeric(length(dat$event_times))
  nobs_vec <- dat$nobs_vec
  
  for (i in 1:nrow(y)) {
    obs_mtx <- matrix(c(1, 0, 1-p[i], p[i]), nrow = 2, byrow = TRUE)
    trans_mtx <- matrix(c(1-theta, theta, 1-theta_prime, theta_prime), nrow = 2, byrow = TRUE)
    
    lprob_obs_given_occu <- dHMM_compiled(y[i, 1:nobs_vec[i]] + 1, 
                                          init = c(1-theta_equil, theta_equil), 
                                          probObs = obs_mtx, 
                                          probTrans = trans_mtx,
                                          len = nobs_vec[i],
                                          checkRowSums = TRUE,
                                          log = TRUE)
    
    if (site_has_no_obs[i]) {
      # Site prob. is the (prob of the data|occupied) * prob occupied, plus prob it is not occupied
      ll_per_site[i] <- log((1-psi[i]) + exp(lprob_obs_given_occu) * psi[i])
    } else {
      # Site prob. is (prob of the data|occupied) * prob occupied
      ll_per_site[i] <- lprob_obs_given_occu + log(psi[i])
    }
  }
  
  return(sum(ll_per_site))
}



#' @name fit_cluster_occ_mod
#' @param dat A named list of covariate and detection data 
#' @param maxit Maximum number of optimization iterations (an optim control
#'              parameter)
#' @description This function estimates the clustered occupancy model for the
#'              simulation study. It calls R's optim and the log-likelihood
#'              objective function cluster_occ_ll.
fit_cluster_occ_mod <- function(dat, maxit = 10000) {
  
  site_has_no_obs <- as.numeric(rowSums(dat$y_base, na.rm=T) == 0)
  
  optimized_soln <- optim(par = c(0, 0, 0, 0, 0, 0),
                          fn = cluster_occ_ll, 
                          gr = NULL,
                          dat = dat,
                          site_has_no_obs = site_has_no_obs,
                          control = list(fnscale = -1,
                                         maxit = maxit))
  
  re_optimized_soln <- optim(par = optimized_soln$par,
                             fn = cluster_occ_ll, gr = NULL,
                             dat = dat,
                             hessian = TRUE,
                             site_has_no_obs = site_has_no_obs,
                             control = list(fnscale = -1,
                                            maxit = maxit),
                             method = "BFGS")
  
  std_error <- rep(NA, length(re_optimized_soln$par))
  tryCatch({
    std_error <- sqrt(diag(solve(-re_optimized_soln$hess)))
  }, error = function(e) { })
  
  
  return(list(summary = data.frame(
    param = c("logit_p_mu", "p_beta2", 
              "logit_psi_mu", "psi_beta1", "logit_theta", "logit_theta_prime"),
    submodel = c(
      "det", "det", "occ", "occ", NA, NA
    ),
    Estimate = re_optimized_soln$par,
    SE = std_error
  ), 
  convergence = optimized_soln$convergence,
  ll = optimized_soln$value
  ))
}


#' @name sim_and_evaluate_one_mod_clustered_par
#' @param x The row of "design_df" that we're on
#' @param design_df The data frame that defines the simulation study design
#' @param row The row of "design_df" that we're on
#' @param do_gof If TRUE, the join count test is conducted on the SOM; otherwise
#'              it is not
#' @param fit_clustered_mod If TRUE, a clustered model is fit to the data;
#'              otherwise, the clustered model is skipped
#' @param datasummary_only If TRUE, data are simulated and summarized, and then
#'              the function ends without estimating any models
#' @param drop_intervals Passed to "sim_occu_clustered"
#' @description A wrapper around sim_and_evaluate_one_mod_clustered that makes
#'               parallelization easier.
sim_and_evaluate_one_mod_clustered_par <- function(x, design_df = NULL, do_gof, 
                                                   datasummary_only = FALSE,
                                                   fit_clustered_mod = TRUE,
                                                   drop_intervals = FALSE) {
  
  
  
  if (is.null(design_df)) {
    design_df <- all_iters_df
  }
  
  tryCatch({
    sim_and_evaluate_one_mod_clustered(
      iter = design_df$iter[x], row = x, design_df = design_df,
      datasummary_only = datasummary_only,
      do_gof = do_gof, fit_clustered_mod = fit_clustered_mod, 
      drop_intervals = drop_intervals
    )
  }, error = function(e) {
    return(list(datasummary = data.frame(
      avg_dets_per_occ_site = "Errored",
      avg_occ_sites_w_dets = "Errored",
      total_dets = "Errored",
      total_occ_sites = "Errored",
      iter = iter, scenario = design_df$scenario[x]
    )
    ))
  })
}



