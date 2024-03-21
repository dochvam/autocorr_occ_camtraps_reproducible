#########################################################################
#' vis_fn.R
#' 
#' This file holds helper functions for use in visualizations
#########################################################################

#' @name get_legend
#' @description Extract the legend from a ggplot object
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#' @name make_one_pane
#' @description Make a single pane of a multifaceted heatmap grid plot for vis.R
make_one_pane <- function(df, nsite_level) {
  pane <- df %>% 
    filter(nsite == nsite_level, det_breakdown == det_breakdown_ref) %>% 
    ggplot() +
    geom_tile(aes(as.factor(overall_det_prob), 
                  as.factor(nimble::expit(logit_psi_mu)),
                  fill = target_col)) +
    facet_grid(effort_compare_name~autocorr_str_named) +
    theme_minimal() +
    xlab("Detection rate") + ylab("Mean occupancy prob.") +
    ggtitle(paste0(nsite_level, " sites"))
  
  pane
}

#' @name make_big_grid_plot
#' @description Make multifaceted heatmap grid plots for vis.R
make_big_grid_plot <- function(df, target_col, target_name = "RMSE", min_max = NULL,
                               divergent = FALSE, discrete = FALSE, title = NULL,
                               use_threshold = FALSE, ncol = 1, midpoint = 0) {
  df$target_col <- unlist(df[, target_col])
  thresholded <- FALSE
  if (!discrete && !is.null(threshold) && use_threshold) {
    thresholded <- TRUE
    indices <- which(abs(df$target_col) > threshold)
    
    if (length(indices) > 0) {
      df$target_col[indices] <- sign(df$target_col[indices]) * threshold
    }
    
    if (!is.null(min_max) && !is.na(min_max[1]) && min_max[1] < -threshold) min_max[1] <- -threshold
    if (!is.null(min_max) && !is.na(min_max[2]) && min_max[2] >  threshold) min_max[2] <-  threshold
  }
  
  if (length(unique(df$autocorr_str)) > 1) {
    df$autocorr_str_named <- factor(
      paste0("Autocorr = ", df$autocorr_str),
      levels = paste0("Autocorr = ", sort(unique(scenario_grid$autocorr_str)))
    )
  } else {
    df$autocorr_str_named <- ""
  }
  
  if ("effort_compare" %in% colnames(df)) {
    df$effort_compare_name <- factor(
      paste0(df$effort_compare, " days"),
      levels = paste0(sort(unique(scenario_grid$n_true_event)), " days")
      
    )
  } else {
    df$effort_compare_name <- factor(
      paste0(df$n_true_event, " days"),
      levels = paste0(sort(unique(scenario_grid$n_true_event)), " days")
    )
  }
  
  if (discrete) {
    df$target_col <- factor(df$target_col, levels = names(modtype_colors))
  }
  
  pane1 <- make_one_pane(df, nsite_level = 40)
  pane2 <- make_one_pane(df, nsite_level = 80)
  pane3 <- make_one_pane(df, nsite_level = 120)
  
  if (discrete) {
    pane1 <- pane1 + scale_fill_manual(target_name, values = modtype_colors) 
    legend <- get_legend(pane1)
    pane1 <- pane1 + theme(legend.position = "none")
    pane2 <- pane2 + scale_fill_manual(target_name, values = modtype_colors) +
      theme(legend.position = "none")
    pane3 <- pane3 + scale_fill_manual(target_name, values = modtype_colors) +
      theme(legend.position = "none")
    
    complete_plot <- gridExtra::arrangeGrob(
      pane1, pane2, pane3, ncol = 1, top = title, right = legend
    )
    return(complete_plot)
  } else if (divergent) {
    pane1 <- pane1 +
      scale_fill_gradient2(target_name, high = scales::muted("red"), mid = "gray",
                           low = scales::muted("blue"), midpoint = midpoint,
                           limits = min_max)
    pane2 <- pane2 +
      scale_fill_gradient2(target_name, high = scales::muted("red"), mid = "gray",
                           low = scales::muted("blue"), midpoint = midpoint,
                           limits = min_max)
    pane3 <- pane3 +
      scale_fill_gradient2(target_name, high = scales::muted("red"), mid = "gray",
                           low = scales::muted("blue"), midpoint = midpoint,
                           limits = min_max)
  } else {
    pane1 <- pane1 + scale_fill_viridis_c(target_name, limits = min_max)
    pane2 <- pane2 + scale_fill_viridis_c(target_name, limits = min_max)
    pane3 <- pane3 + scale_fill_viridis_c(target_name, limits = min_max)
  }
  
  if (ncol > 1 && !is.null(min_max)) {
    legend_all <- get_legend(pane1)
    
    pane1 <- pane1 + theme(legend.position = "none")
    pane2 <- pane2 + theme(legend.position = "none")
    pane3 <- pane3 + theme(legend.position = "none")
  } else {
    legend_all <- ""
  }
  
  complete_plot <- gridExtra::arrangeGrob(
    pane1, pane2, pane3, ncol = ncol, top = title, right = legend_all
  )
  
  complete_plot
}


#' @name make_big_big_grid_plot
#' @description Make giant multifaceted heatmap grid plots for vis.R
make_big_big_grid_plot <- function(df, target_col, target_name = "RMSE",
                                   windows = c(1, 3, 7),
                                   gaps = c(0, 0, 0),
                                   divergent = FALSE, discrete = FALSE, 
                                   title = NULL, use_threshold = FALSE,
                                   minmax_symmetric = FALSE,
                                   midpoint = 0) {
  if (length(windows) != length(gaps)) stop("Length of windows and gaps must match")
  
  plots_list <- list()
  
  overall_df <- df %>% 
    filter(det_breakdown == det_breakdown_ref) %>% 
    filter(paste0(window, "_", gap) %in% paste0(windows, "_", gaps))
  
  scale_min_max <- numeric(2) 
  scale_min_max[1] <- min(overall_df[[target_col]])
  scale_min_max[2] <- max(overall_df[[target_col]])
  if (minmax_symmetric) {
    if (abs(scale_min_max[1]) > scale_min_max[2]) {
      scale_min_max[2] <- -1 * scale_min_max[1]
    } else {
      scale_min_max[1] <- -1 * scale_min_max[1]
    }
  }
  
  for (i in 1:length(windows)) {
    this_df <- overall_df %>% filter(window == windows[i], gap == gaps[i])
    plots_list[[i]] <- make_big_grid_plot(this_df, target_col, target_name, divergent,
                                          min_max = scale_min_max,
                                          use_threshold = use_threshold,
                                          midpoint = midpoint,
                                          discrete = discrete, 
                                          title = paste0(windows[i], "-", gaps[i]))
    
  }
  
  return(gridExtra::arrangeGrob(grobs = plots_list, nrow = 1, top = title))
}


#' @name make_plots_bywindow
#' @description Make the window plot for vis.R
make_plots_bywindow <- function(df, outcome_column, outcome_column_name,
                                hline = 0, summary_fn = mean, title = "") {
  
  colnames(df)[colnames(df) == outcome_column] <- "outcome"
  
  plot1 <- df %>% 
    group_by(scenario) %>%
    summarize(mean_outcome = summary_fn(outcome, na.rm = T)) %>%
    left_join(scenario_grid) %>% 
    filter(window <= 20) %>%
    filter(window + gap <= n_true_event/2) %>% 
    mutate(scenario = parse_number(scenario)) %>% 
    mutate(autocorr_str = ifelse(autocorr_str == 10, "Strong autocorr.",
                                 "No autocorr.")) %>% 
    arrange(scenario) %>% 
    ggplot(aes(window, mean_outcome, colour = as.factor(gap),
               fill = as.factor(gap), group = as.factor(gap))) +
    geom_point(show.legend = TRUE) +
    # geom_smooth() +
    facet_grid(paste0("Effort = ", n_true_event, " days")~autocorr_str) + 
    xlab("Detection window") + 
    ylab(outcome_column_name) +
    scale_x_continuous(breaks = 1:10*2) +
    theme_minimal() +
    theme(legend.position = "right", panel.grid.minor = element_blank()) +
    scale_color_viridis_d("Det. gap", end = 0.8) +
    scale_fill_viridis_d("Det. gap", end = 0.8) +
    ggtitle(title)
  
  
  df <- df %>% 
    left_join(scenario_grid) %>% 
    filter(window + gap <= n_true_event/2)
  ylim <- quantile(df$outcome, probs = c(0.005, 0.995), na.rm = T)
  # Boxplot
  plot2 <- df %>% 
    mutate(scenario = parse_number(scenario)) %>% 
    mutate(autocorr_str = ifelse(autocorr_str == 10, "Strong autocorr.",
                                 "No autocorr.")) %>% 
    arrange(scenario) %>% 
    ggplot() +
    geom_boxplot(aes(x = factor(window), y = outcome), outlier.shape = NA) +
    facet_grid(paste0("Effort = ", n_true_event, " days")~autocorr_str) + 
    xlab("Detection window") + ylab(outcome_column_name) +
    # scale_x_continuous(breaks = 1:5*2) +
    theme_minimal() +
    scale_color_viridis_d("Det prob.", end = 0.7) +
    ylim(ylim)
  
  
  if (hline != "none") {
    plot1 <- plot1 + geom_hline(yintercept = hline)
    plot2 <- plot2 + geom_hline(yintercept = hline)
  }
  
  return(list(dotplot=plot1, boxplot=plot2))
}


#' @name dethist_mtx_to_long
#' @param mtx A matrix-form detection history such as the `y` slot in an
#'            unmarkedFrameOccu object.
#' @description This function reformats a matrix-form detection history to
#'              long form.
dethist_mtx_to_long <- function(mtx, tag = NULL) {
  mtx %>% 
    as.data.frame() %>% 
    mutate(site = row_number()) %>% 
    pivot_longer(all_of(paste0("V", 1:effort_days))) %>% 
    group_by(site) %>% 
    mutate(visit = row_number()) %>% 
    ungroup() %>% 
    mutate(tag = tag)
}



#' @name get_unmarked_summary
#' @param this_fit A fitted unmarked model, as returned by occu
#' @description This function makes a data frame summarizing the results
#'    of a call to occu.
get_unmarked_summary <- function(this_fit) {
  
  cap <- capture.output(state_summary <- summary(this_fit)$state)
  cap <- capture.output(det_summary <- summary(this_fit)$det)
  
  ests <- bind_rows(
    det_summary   %>% mutate(submodel = "Detection"),
    state_summary %>% mutate(submodel = "Occupancy")
  )
  
  ests$param <- rownames(ests)
  rownames(ests) <- NULL
  ests$param[grepl("Intercept", ests$param)] <- "Intercept"
  
  return(ests)
}

