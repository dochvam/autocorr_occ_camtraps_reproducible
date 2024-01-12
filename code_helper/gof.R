#########################################################################
#' Author: Benjamin R. Goldstein
#' gof.R
#' 
#' This file holds helper functions that support estimating goodness-of-fit
#' of the occupancy model, specifically functions that conduct the join count
#' goodness-of-fit test.
#########################################################################

#' @name get_jc_adjacent_onesite
#' @param x A vector of detection-nondetection data
#' @description Calculates the join count statistic for one site using the
#'              adjacency rule.
get_jc_adjacent_onesite <- function(x) {
  obs <- x[!is.na(x)]
  sum(obs[1:(length(obs)-1)] == 1 &
        obs[2:(length(obs))] == 1)
}

#' @name get_jc_all_onesite
#' @param x A vector of detection-nondetection data
#' @description Calculates the join count statistic for one site using the
#'              all site rule. We do not discuss the all site rule in this
#'              manuscript, but see Wright et al. 2016.
get_jc_all_onesite <- function(x) {
  obs <- x[!is.na(x)]

  choose(sum(obs), 2)
}

#' @name get_join_count_stat_from_y_mtx
#' @param y_mtx A detection history matrix. Each row is a site, and each
#'              column represents replicate visits. Does not support gaps in 
#'              sampling (i.e. if 10 sampling occasions occur at a site, they 
#'              should be represented in columns 1-10.)
#' @param defn The join count rule to use; either "adjacent" or "all"
#' @description Calculates the join count statistic for an entire detection
#'              history, formatted as a matrix
get_join_count_stat_from_y_mtx <- function(y_mtx, defn = "adjacent") {
  
  y_mtx_raw <- y_mtx
  y_mtx <- y_mtx[rowSums(!is.na(y_mtx)) > 1,]
  
  if (defn == "adjacent") {
    jc_neigh <- apply(y_mtx, 1, get_jc_adjacent_onesite) %>% 
      sum()
  } else if (defn == "all") {
    jc_neigh <- apply(y_mtx, 1, get_jc_all_onesite) %>% 
      sum()
  } else {
    stop(paste0("Invalid defn: ", defn))
  }
  
  jc_neigh
}


#' @name run_jc
#' @param fm A fitted unmarked model as returned by `occu`.
#' @param nsim Number of simulated datasets for the parametric bootstrap
#' @description Runs the join count statistic parametric bootstrap using the
#'              adjacent scoring rule, following Wright et al. 2016.
run_jc <- function(fm, nsim) {
  observed_jc <- get_join_count_stat_from_y_mtx(getY(fm@data), "adjacent")
  
  newdat <- unmarked::simulate(fm, nsim = nsim)
  
  jc_adj_boot <- numeric(nsim)
  for (i in 1:nsim) {
    jc_adj_boot[i] <- get_join_count_stat_from_y_mtx(newdat[[i]], "adjacent")
  }
  
  return(list(
    observed = observed_jc,
    jc_bootstrapped = jc_adj_boot,
    pval = mean(jc_adj_boot > observed_jc)
  ))
}



#' @name get_rmse_unmarked
#' @param this_fit A fitted unmarked occu model
#' @param dat The simulated data (which include true psi)
#' @param which_y Are we using the simulated data ("Base") or the alternative
#'                window/gap transformed data
#' @description Calculates the RMSE of a fitted unmarked occu model
get_rmse_unmarked <- function(this_fit, dat, which_y) {
  if (which_y == "Base" || which_y == 1) {
    umf <- unmarked::unmarkedFrameOccu(
      y = dat$y_base,
      siteCovs = data.frame(
        cov1 = dat$cov1, 
        cov2 = dat$cov2
      )
    )
  } else {
    umf <- unmarked::unmarkedFrameOccu(
      y = dat$y_alt,
      siteCovs = data.frame(
        cov1 = dat$cov1, 
        cov2 = dat$cov2
      )
    )
  }
  
  prediction <- predict(this_fit, newdata = umf, type = "state")
  
  sqrt(mean((prediction$Predicted - dat$true_psi)^2))
}


#' @name get_rmse_unmarked
#' @param this_fit A fitted clustered occu model as returned by 
#'                 fit_cluster_occ_mod
#' @param dat The simulated data (which include true psi)
#' @description Calculates the RMSE of a fitted clustered occu model
get_rmse_clustered <- function(this_fit, dat) {
  prediction <- expit(
    this_fit$summary$Estimate[this_fit$summary$param == "logit_psi_mu"] +
      this_fit$summary$Estimate[this_fit$summary$param == "psi_beta1"] *
      dat$cov1
  )
  
  sqrt(mean((prediction - dat$true_psi)^2))
}

