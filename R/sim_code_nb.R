# ┌───────────────────────────────────────┐
# │          Study 2 Simulation           │
# └───────────────────────────────────────┘
sim_start <- Sys.time()
# Import helper functions, constants, objects, packages, and set-up detail
source("setup_nb.R")

# ┌── Info output and logging ───────────────────────────┐
message("Total cores assigned: ", total_cores, 
        " | SimDesign core used: ", n_sim_cores, 
        " | brms cores used: ", brms_cores,
        " | brms chains: ", n_brms_chains, 
        " | threads per chain: ", n_brms_threads)
# └───────────────────────────────────────────────────────┘

# ┌── SLURM array ID logic ───────────────────────────┐
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Please provide SLURM_ARRAY_TASK_ID")
array_id <- as.integer(args[1])
# └───────────────────────────────────────────────────────┘

# ┌── Conditions ───────────────────────────────────────────┐
Design <- createDesign(n =c(10, 
                            20, 
                            50, 
                            100, 
                            250),
                       n_clusters= c(12,
                                     48,
                                     128),
                       null_prop= c(0, 0.25, 0.5, 0.75, 1),  
                       nonnull_effect = c(1, 5), 
                       sigma = 1, 
                       subset = !(null_prop == 1 & nonnull_effect %in% c(5)))
# └───────────────────────────────────────────────────────┘




# ┌── Generate ─────────────────────────────────────────────────┐
Generate <- function(condition, fixed_objects = fixed_objects) {
  
  Attach(condition) # make n, n_clusters, sigma, and sd_int directly accessible
  Attach(fixed_objects) # make intercept, MTeffect, sd_MT, cor, abd alpha directly accessible
  
  # Total sample size
  N <- n*n_clusters
  # Calculate the number of null-effect clusters
  n_null_clusters <- n_clusters * null_prop
  
  n_MT1 <- n/2
  # Generate data for multiple clusters with dynamic MTeffect
  dat <- do.call(rbind, lapply(1:n_clusters, function(cluster_id) {
    # Set MTeffect dynamically based on null_prop
    MTeffect <- ifelse(cluster_id <= n_null_clusters, 0, nonnull_effect)
    # Simulate data for the current cluster
    dat <- simulate_cluster_data(n, n_MT1, intercept, MTeffect, sigma)
    # Add cluster number
    dat$cluster <- cluster_id
    return(dat)
  }))
  
  dat
}
# └───────────────────────────────────────────────────────┘


# ┌── Analyse ─────────────────────────────────────────────────┐

Analyse <- function(condition, dat, fixed_objects=fixed_objects) {
  
  Attach(condition)
  Attach(dat)
  Attach(fixed_objects)
  
  
  # Define null and non-null cluster indices
  n_null_clusters <- n_clusters * null_prop
  if (n_null_clusters %% 1 != 0) stop("n_null_clusters must be a whole number")
  true_null <- if (n_null_clusters > 0) seq_len(n_null_clusters) else integer(0)
  true_non_null <- if (n_null_clusters < n_clusters) seq(n_null_clusters + 1, n_clusters) else integer(0)
  
  ######### TRADITIONAL APPROACHES #########
  
  traditional_res <- by(dat, dat$cluster, t_test_func)
  traditional_res <- data.frame(do.call(rbind, traditional_res))
  p.vals <- traditional_res$p.value
  
  if (length(p.vals) == 0 || all(is.na(p.vals)) || all(p.vals < 0 | p.vals > 1, na.rm = TRUE)) {
    stop("Invalid p-values.")
  }
  
  unadjusted_error <- if (length(true_null) > 0) mean(p.vals[true_null] < alpha) else NAflag
  unadjusted_family_error <- if (length(true_null) > 0) any(p.vals[true_null] < alpha) else NAflag
  unadjusted_power <- if (length(true_non_null) > 0) mean(p.vals[true_non_null] < alpha) else NAflag
  
  FWE_p.vals <- p.adjust(p.vals, method = "bonferroni")
  if (any(is.na(FWE_p.vals)) || all(FWE_p.vals < 0 | FWE_p.vals > 1, na.rm = TRUE)) {
    stop("Invalid FWE p-values.")
  }
  FWE_error <- if (length(true_null) > 0) mean(FWE_p.vals[true_null] < alpha) else NAflag
  FWE_family_error <- if (length(true_null) > 0) any(FWE_p.vals[true_null] < alpha) else NAflag
  FWE_power <- if (length(true_non_null) > 0) mean(FWE_p.vals[true_non_null] < alpha) else NAflag
  
  FDR_p.vals <- p.adjust(p.vals, method = "fdr")
  if (any(is.na(FDR_p.vals)) || all(FDR_p.vals < 0 | FDR_p.vals > 1, na.rm = TRUE)) {
    stop("Invalid FDR p-values.")
  }
  FDR_error <- if (length(true_null) > 0) mean(FDR_p.vals[true_null] < alpha) else NAflag
  FDR_family_error <- if (length(true_null) > 0) any(FDR_p.vals[true_null] < alpha) else NAflag
  FDR_power <- if (length(true_non_null) > 0) mean(FDR_p.vals[true_non_null] < alpha) else NAflag
  
  
  
  ######### FREQUENTIST MLM #########
  
  
  for (opt in c("bobyqa", "Nelder_Mead", "nloptwrap")) {
    mod <- fit_lmer_with_control(dat, opt, maxfun = 1e5)
    if (!is.null(mod) && isTRUE(mod@optinfo$conv$opt == 0)) break
  }
  if (is.null(mod) || isTRUE(mod@optinfo$conv$opt != 0)) {
    stop("All optimizers failed.")
  }
  
  
  randCIs <- process_CI_results(mod)
  if (anyNA(randCIs$UB) || anyNA(randCIs$LB) || any(is.nan(randCIs$UB)) || any(is.nan(randCIs$LB)) ||
      (mean(abs(randCIs$LB - median(randCIs$LB)) < 1e-6) >= 0.5) ||
      (mean(abs(randCIs$UB - median(randCIs$UB)) < 1e-6) >= 0.5)) {
    stop("MLM model estimation issues detected.")
  }
  
  MLM_CI_prod <- randCIs$LB * randCIs$UB
  MLM_error <- if (length(true_null) > 0) mean(MLM_CI_prod[true_null] > 0) else NAflag
  MLM_family_error <- if (length(true_null) > 0) any(MLM_CI_prod[true_null] > 0) else NAflag
  MLM_power <- if (length(true_non_null) > 0) mean(MLM_CI_prod[true_non_null] > 0) else NAflag
  
  
  
  
  
  ######### BAYESIAN MLM(s) #########
  
  # Default Prior
  start <- Sys.time()
  message("----- Starting default model fit -----")
  
  bmlm_default <- fit_brms_model(data=dat, 
                                 model_tag = prior_list$brms_default$model_tag,
                                 exe_dir = exe_dir)
  
  if (is.null(bmlm_default)) {
    message("Default model fitting failed after all attempts.")} 
  
  end <- Sys.time()
  message("----- Default model fit completed -----")
  message("Elapsed time (minutes): ", round(difftime(end, start, units = "mins"), 2))
  
  brandCIs <- process_CI_results(bmlm_default)
  
  if (any(is.na(brandCIs$UB)) || 
      any(is.na(brandCIs$LB)) || 
      any(is.nan(brandCIs$UB)) || 
      any(is.nan(brandCIs$LB)) || 
      (mean(abs(brandCIs$LB - median(brandCIs$LB)) < 1e-6) >= 0.5) ||
      (mean(abs(brandCIs$UB - median(brandCIs$UB)) < 1e-6) >= 0.5)) {
    stop("Default model estimation issues persist: NA/NaN values, identical CIs, convergence problems, or singular fit.")}
  
  bMLM_CI_prod <- brandCIs$LB * brandCIs$UB
  
  bMLM_error <- if (length(true_null) > 0) mean(bMLM_CI_prod[true_null] > 0) else NAflag
  bMLM_family_error <- if (length(true_null) > 0) any(bMLM_CI_prod[true_null] > 0) else NAflag
  bMLM_power <- if (length(true_non_null) > 0) mean(bMLM_CI_prod[true_non_null] > 0) else NAflag
  
  
  
  ############################### @@@@@ Very informative null @@@@@ ##############################
  
  # VIN Prior
  
  start <- Sys.time()
  message("----- Starting VIN model fit -----")
  
  brms_vin <- fit_brms_model(data=dat, 
                             model_tag = prior_list$very_info_nill$model_tag,
                             exe_dir   =   exe_dir)
  
  if (is.null(brms_vin)) {
    message("VIN model fitting failed after all attempts.")} 
  
  end <- Sys.time()
  message("----- VIN model fit completed -----")
  message("Elapsed time (minutes): ", round(difftime(end, start, units = "mins"), 2))
  
  brandCIs_vin <- process_CI_results(brms_vin)
  
  if (any(is.na(brandCIs_vin$UB)) || 
      any(is.na(brandCIs_vin$LB)) || 
      any(is.nan(brandCIs_vin$UB)) || 
      any(is.nan(brandCIs_vin$LB)) || 
      (mean(abs(brandCIs_vin$LB - median(brandCIs_vin$LB)) < 1e-6) >= 0.5) ||
      (mean(abs(brandCIs_vin$UB - median(brandCIs_vin$UB)) < 1e-6) >= 0.5)) {
    stop("VIN model estimation issues persist: NA/NaN values, identical CIs, convergence problems, or singular fit.")}
  
  bprod <- brandCIs_vin$LB * brandCIs_vin$UB
  
  bMLM_vin_error <- if (length(true_null) > 0) mean(bprod[true_null] > 0) else NAflag
  bMLM_vin_family_error <- if (length(true_null) > 0) any(bprod[true_null] > 0) else NAflag
  bMLM_vin_power <- if (length(true_non_null) > 0) mean(bprod[true_non_null] > 0) else NAflag
  
  
  ############################### @@@@@ Weakly informative null @@@@@ ##############################
  
  start <- Sys.time()
  message("----- Starting WIN model fit -----")
  
  brms_win <- fit_brms_model(data=dat, 
                             model_tag = prior_list$weak_info_nill$model_tag,
                             exe_dir = exe_dir)
  
  if (is.null(brms_win)) {
    message("WIN model fitting failed after all attempts.")} 
  
  end <- Sys.time()
  message("----- WIN model fit completed -----")
  message("Elapsed time (minutes): ", round(difftime(end, start, units = "mins"), 2))
  
  brandCIs_win <- process_CI_results(brms_win)
  
  
  if (any(is.na(brandCIs_win$UB)) || 
      any(is.na(brandCIs_win$LB)) || 
      any(is.nan(brandCIs_win$UB)) || 
      any(is.nan(brandCIs_win$LB)) || 
      (mean(abs(brandCIs_win$LB - median(brandCIs_win$LB)) < 1e-6) >= 0.5) ||
      (mean(abs(brandCIs_win$UB - median(brandCIs_win$UB)) < 1e-6) >= 0.5)) {
    stop("WIN model estimation issues persist: NA/NaN values, identical CIs, convergence problems, or singular fit.")}
  
  bprod_win <- brandCIs_win$LB * brandCIs_win$UB
  
  bMLM_win_error <- if (length(true_null) > 0) mean(bprod_win[true_null] > 0) else NAflag
  bMLM_win_family_error <- if (length(true_null) > 0) any(bprod_win[true_null] > 0) else NAflag
  bMLM_win_power <- if (length(true_non_null) > 0) mean(bprod_win[true_non_null] > 0) else NAflag
  
  

  
  
  
  ############################### @@@@@ HORSESHOE @@@@@ ##############################
  
  start <- Sys.time()
  message("----- Starting HS model fit -----")
  
  
  brms_hs <- fit_brms_model(data=dat, 
                            model_tag = prior_list$hs$model_tag,
                            exe_dir   = exe_dir)
  
  if (is.null(brms_hs)) {
    message("HS model fitting failed after all attempts.")
  } 
  
  end <- Sys.time()
  message("----- HS model fit completed -----")
  message("Elapsed time (minutes): ", round(difftime(end, start, units = "mins"), 2))
  
  brandCIs_hs <- process_CI_results(brms_hs)
  
  
  if (any(is.na(brandCIs_hs$UB)) || 
      any(is.na(brandCIs_hs$LB)) || 
      any(is.nan(brandCIs_hs$UB)) || 
      any(is.nan(brandCIs_hs$LB)) || 
      (mean(abs(brandCIs_hs$LB - median(brandCIs_hs$LB)) < 1e-6) >= 0.5) ||
      (mean(abs(brandCIs_hs$UB - median(brandCIs_hs$UB)) < 1e-6) >= 0.5)) {
    stop("HS model estimation issues persist: NA/NaN values, identical CIs, convergence problems, or singular fit.")}
  
  bprod_hs <- brandCIs_hs$LB * brandCIs_hs$UB
  
  bMLM_hs_error <- if (length(true_null) > 0) mean(bprod_hs[true_null] > 0) else NAflag
  bMLM_hs_family_error <- if (length(true_null) > 0) any(bprod_hs[true_null] > 0) else NAflag
  bMLM_hs_power <- if (length(true_non_null) > 0) mean(bprod_hs[true_non_null] > 0) else NAflag
  
  
  
  ############################### @@@@@ SPIKE AND SLAB @@@@@ ##############################
  
  start <- Sys.time()
  message("----- Starting SPIKE & SLAB model fit -----")
  
  brms_ss <- fit_brms_model(data=dat, 
                            model_tag = prior_list$spike_slab$model_tag,
                            exe_dir   = exe_dir)
  # Check the summary if successful
  if (is.null(brms_ss)) {
    message("SPIKE & SLAB Model fitting failed after all attempts.")} 
  
  end <- Sys.time()
  message("----- SPIKE & SLAB model fit completed -----")
  message("Elapsed time (minutes): ", round(difftime(end, start, units = "mins"), 2))
  
  
  # Recalculate random effects and confidence intervals
  brandCIs_ss <- process_CI_results(brms_ss)
  
  # Check for estimation issues (e.g., NA/NaN values, CIs w/o variability, convergence problems)
  if (any(is.na(brandCIs_ss$UB)) || 
      any(is.na(brandCIs_ss$LB)) || 
      any(is.nan(brandCIs_ss$UB)) || 
      any(is.nan(brandCIs_ss$LB)) || 
      (mean(abs(brandCIs_ss$LB - median(brandCIs_ss$LB)) < 1e-6) >= 0.5) ||
      (mean(abs(brandCIs_ss$UB - median(brandCIs_ss$UB)) < 1e-6) >= 0.5)) {
    stop("SPIKE & SLAB model estimation issues persist: NA/NaN values, identical CIs, convergence problems, or singular fit.")
  }
  
  bprod_ss <- brandCIs_ss$LB * brandCIs_ss$UB
  
  bMLM_ss_error <- if (length(true_null) > 0) mean(bprod_ss[true_null] > 0) else NAflag
  bMLM_ss_family_error <- if (length(true_null) > 0) any(bprod_ss[true_null] > 0) else NAflag
  bMLM_ss_power <- if (length(true_non_null) > 0) mean(bprod_ss[true_non_null] > 0) else NAflag
  
  
  
  ############################### @@@@@ Tau half t @@@@@ ##############################
  
  
  start <- Sys.time()
  message("----- Starting TAU half-t model fit -----")
  
  brms_tau_t <- fit_brms_model(data=dat, 
                               model_tag = prior_list$tau_sd_MT_t$model_tag,
                               exe_dir   = exe_dir)
  
  if (is.null(brms_tau_t)) {
    message("Tau half t model fitting failed after all attempts.")
  } 
  
  end <- Sys.time()
  message("----- brms_tau_t model fit completed -----")
  message("Elapsed time (minutes): ", round(difftime(end, start, units = "mins"), 2))
  
  
  brandCIs_ts <- process_CI_results(brms_tau_t)
  
  if (any(is.na(brandCIs_ts$UB)) || 
      any(is.na(brandCIs_ts$LB)) || 
      any(is.nan(brandCIs_ts$UB)) || 
      any(is.nan(brandCIs_ts$LB)) || 
      (mean(abs(brandCIs_ts$LB - median(brandCIs_ts$LB)) < 1e-6) >= 0.5) ||
      (mean(abs(brandCIs_ts$UB - median(brandCIs_ts$UB)) < 1e-6) >= 0.5)) {
    stop("brms_tau_t model estimation issues persist: NA/NaN values, identical CIs, convergence problems, or singular fit.")}
  
  bprod_ts <- brandCIs_ts$LB * brandCIs_ts$UB
  
  bMLM_tau_t_error <- if (length(true_null) > 0) mean(bprod_ts[true_null] > 0) else NAflag
  bMLM_tau_t_family_error <- if (length(true_null) > 0) any(bprod_ts[true_null] > 0) else NAflag
  bMLM_tau_t_power <- if (length(true_non_null) > 0) mean(bprod_ts[true_non_null] > 0) else NAflag
  
  
  
  ############################### @@@@@ Tau exp @@@@@ ##############################
  
  start <- Sys.time()
  message("----- Starting Tau exp model fit -----")
  
  brms_tau_exp <- fit_brms_model(data=dat, 
                                 model_tag = prior_list$tau_sd_MT_exp$model_tag,
                                 exe_dir   = exe_dir)
  
  if (is.null(brms_tau_exp)) {
    message("Tau exp model fitting failed after all attempts.")
  } 
  
  end <- Sys.time()
  message("----- brms_tau_exp model fit completed -----")
  message("Elapsed time (minutes): ", round(difftime(end, start, units = "mins"), 2))
  
  brandCIs_exp <- process_CI_results(brms_tau_exp)
  
  if (any(is.na(brandCIs_exp$UB)) || 
      any(is.na(brandCIs_exp$LB)) || 
      any(is.nan(brandCIs_exp$UB)) || 
      any(is.nan(brandCIs_exp$LB)) || 
      (mean(abs(brandCIs_exp$LB - median(brandCIs_exp$LB)) < 1e-6) >= 0.5) ||
      (mean(abs(brandCIs_exp$UB - median(brandCIs_exp$UB)) < 1e-6) >= 0.5)) {
    stop("Tau exp model estimation issues persist: NA/NaN values, identical CIs, convergence problems, or singular fit.")}
  
  bprod_exp <- brandCIs_exp$LB * brandCIs_exp$UB
  
  bMLM_tau_exp_error <- if (length(true_null) > 0) mean(bprod_exp[true_null] > 0) else NAflag
  bMLM_tau_exp_family_error <- if (length(true_null) > 0) any(bprod_exp[true_null] > 0) else NAflag
  bMLM_tau_exp_power <- if (length(true_non_null) > 0) mean(bprod_exp[true_non_null] > 0) else NAflag
  
  
  
  ############################### @@@@@ Tau log @@@@@ ##############################
  
  start <- Sys.time()
  message("----- Starting Tau log model fit -----")
  
  brms_tau_log <- fit_brms_model(data=dat, 
                                 model_tag = prior_list$tau_sd_MT_log$model_tag,
                                 exe_dir   = exe_dir)
  
  if (is.null(brms_tau_log)) {
    message("Tau exp model fitting failed after all attempts.")
  } 
  
  end <- Sys.time()
  message("----- brms_tau_log model fit completed -----")
  message("Elapsed time (minutes): ", round(difftime(end, start, units = "mins"), 2))
  
  brandCIs_log <- process_CI_results(brms_tau_log)
  
  if (any(is.na(brandCIs_log$UB)) || 
      any(is.na(brandCIs_log$LB)) || 
      any(is.nan(brandCIs_log$UB)) || 
      any(is.nan(brandCIs_log$LB)) || 
      (mean(abs(brandCIs_log$LB - median(brandCIs_log$LB)) < 1e-6) >= 0.5) ||
      (mean(abs(brandCIs_log$UB - median(brandCIs_log$UB)) < 1e-6) >= 0.5)) {
    stop("Tau log model estimation issues persist: NA/NaN values, identical CIs, convergence problems, or singular fit.")}
  
  bprod_log <- brandCIs_log$LB * brandCIs_log$UB
  
  bMLM_tau_log_error <- if (length(true_null) > 0) mean(bprod_log[true_null] > 0) else NAflag
  bMLM_tau_log_family_error <- if (length(true_null) > 0) any(bprod_log[true_null] > 0) else NAflag
  bMLM_tau_log_power <- if (length(true_non_null) > 0) mean(bprod_log[true_non_null] > 0) else NAflag
  
  
  
  
  ############################### @@@@@ Tau gamma @@@@@ ##############################
  
  start <- Sys.time()
  message("----- Starting Tau gamma model fit -----")
  
  brms_tau_gamma <- fit_brms_model(data=dat, 
                                   model_tag = prior_list$tau_sd_MT_gamma$model_tag,
                                   exe_dir   = exe_dir)
  
  if (is.null(brms_tau_gamma)) {
    message("Tau exp model fitting failed after all attempts.")
  } 
  
  end <- Sys.time()
  message("----- brms_tau_gamma model fit completed -----")
  message("Elapsed time (minutes): ", round(difftime(end, start, units = "mins"), 2))
  
  brandCIs_gamma <- process_CI_results(brms_tau_gamma)
  
  if (any(is.na(brandCIs_gamma$UB)) || 
      any(is.na(brandCIs_gamma$LB)) || 
      any(is.nan(brandCIs_gamma$UB)) || 
      any(is.nan(brandCIs_gamma$LB)) || 
      (mean(abs(brandCIs_gamma$LB - median(brandCIs_gamma$LB)) < 1e-6) >= 0.5) ||
      (mean(abs(brandCIs_gamma$UB - median(brandCIs_gamma$UB)) < 1e-6) >= 0.5)) {
    stop("Tau log model estimation issues persist: NA/NaN values, identical CIs, convergence problems, or singular fit.")}
  
  bprod_gamma <- brandCIs_gamma$LB * brandCIs_gamma$UB
  
  bMLM_tau_gamma_error <- if (length(true_null) > 0) mean(bprod_gamma[true_null] > 0) else NAflag
  bMLM_tau_gamma_family_error <- if (length(true_null) > 0) any(bprod_gamma[true_null] > 0) else NAflag
  bMLM_tau_gamma_power <- if (length(true_non_null) > 0) mean(bprod_gamma[true_non_null] > 0) else NAflag
  
  
  
  
  ret <- nc(unadjusted_error, unadjusted_family_error, unadjusted_power,
            FWE_error, FWE_family_error, FWE_power,
            FDR_error, FDR_family_error, FDR_power,
            
            MLM_error, MLM_family_error, MLM_power,
            
            bMLM_error, bMLM_family_error, bMLM_power, #1
            
            bMLM_vin_error, bMLM_vin_family_error, bMLM_vin_power, #2
            bMLM_win_error, bMLM_win_family_error, bMLM_win_power, #3
         
            
            bMLM_hs_error, bMLM_hs_family_error, bMLM_hs_power,    #4
            
            bMLM_ss_error, bMLM_ss_family_error, bMLM_ss_power,    #5
            
            bMLM_tau_t_error, bMLM_tau_t_family_error, bMLM_tau_t_power,    #6
            bMLM_tau_exp_error, bMLM_tau_exp_family_error, bMLM_tau_exp_power, #7
            bMLM_tau_log_error, bMLM_tau_log_family_error, bMLM_tau_log_power, #8
            bMLM_tau_gamma_error, bMLM_tau_gamma_family_error, bMLM_tau_gamma_power #19
  )
  
  message("Memory used (MB): ", round(sum(gc()[,2]), 2))
  return(ret)
  
}

# └──────────────────────────────────────────────────────────────┘

# ┌── Summarise ─────────────────────────────────────────────────┐
Summarise <- function(condition, results, fixed_objects = fixed_objects) {
  
  ret <- c(
    unadjusted_error        = mean(results$unadjusted_error),
    unadjusted_family_error = mean(results$unadjusted_family_error),
    unadjusted_power        = mean(results$unadjusted_power),
    
    FWE_error               = mean(results$FWE_error),
    FWE_family_error        = mean(results$FWE_family_error),
    FWE_power               = mean(results$FWE_power),
    
    FDR_error               = mean(results$FDR_error),
    FDR_family_error        = mean(results$FDR_family_error),
    FDR_power               = mean(results$FDR_power),
    
    MLM_error               = mean(results$MLM_error),
    MLM_family_error        = mean(results$MLM_family_error),
    MLM_power               = mean(results$MLM_power),
    # 
    bMLM_error              = mean(results$bMLM_error),
    bMLM_family_error       = mean(results$bMLM_family_error),
    bMLM_power              = mean(results$bMLM_power),
    
    bMLM_vin_error          = mean(results$bMLM_vin_error),
    bMLM_vin_family_error   = mean(results$bMLM_vin_family_error),
    bMLM_vin_power          = mean(results$bMLM_vin_power),
    
    bMLM_win_error          = mean(results$bMLM_win_error),
    bMLM_win_family_error   = mean(results$bMLM_win_family_error),
    bMLM_win_power          = mean(results$bMLM_win_power),
  
    
    bMLM_hs_error          = mean(results$bMLM_hs_error),
    bMLM_hs_family_error   = mean(results$bMLM_hs_family_error),
    bMLM_hs_power          = mean(results$bMLM_hs_power),
    
    
    bMLM_ss_error          = mean(results$bMLM_ss_error),
    bMLM_ss_family_error   = mean(results$bMLM_ss_family_error),
    bMLM_ss_power          = mean(results$bMLM_ss_power),
    
    
    bMLM_tau_t_error          = mean(results$bMLM_tau_t_error),
    bMLM_tau_t_family_error   = mean(results$bMLM_tau_t_family_error),
    bMLM_tau_t_power          = mean(results$bMLM_tau_t_power),
    
    bMLM_tau_exp_error          = mean(results$bMLM_tau_exp_error),
    bMLM_tau_exp_family_error   = mean(results$bMLM_tau_exp_family_error),
    bMLM_tau_exp_power          = mean(results$bMLM_tau_exp_power),
    
    bMLM_tau_log_error        = mean(results$bMLM_tau_log_error), 
    bMLM_tau_log_family_error = mean(results$bMLM_tau_log_family_error),
    bMLM_tau_log_power        = mean(results$bMLM_tau_log_power),
    
    bMLM_tau_gamma_error        = mean(results$bMLM_tau_gamma_error), 
    bMLM_tau_gamma_family_error = mean(results$bMLM_tau_gamma_family_error), 
    bMLM_tau_gamma_power        = mean(results$bMLM_tau_gamma_power) 
    
  )
 
  ret
}
# └──────────────────────────────────────────────────────────────┘






# ┌── Simulation global setup ───────────────────────────────────┐
#reps_per_cond <- 3L   # test small sim/Total reps per condition
#n_jobs <- ceiling(nrow(Design) / fo$chunk_size) # Number of SLURM array jobs, 
iseed <- 285544376  # Initial seed (use genSeeds() once)
# └──────────────────────────────────────────────────────────────┘


# ┌── Simulation local setup ───────────────────────────┐
#replications <- rep(reps_per_cond, n_jobs)  
arrayID <- SimDesign::getArrayID(type = 'slurm')
#ix <- array2row(arrayID, k = fo$chunk_size, Design)
# └─────────────────────────────────────────────────────┘


Design_Exp <- expandDesign(Design, repeat_conditions = 20L)

# ┌── Run simulation! ───────────────────────────┐
res <- runArraySimulation(
  design = Design_Exp,
  fixed_objects = fixed_objects,
  replications = 50L,
  generate = Generate,
  analyse = Analyse,
  summarise = Summarise,
  iseed = iseed,
  arrayID = arrayID,
  
  max_errors = 100L,
  
  parallel = TRUE, 
 
  #################### array2row = array2row(arrayID, k = fo$chunk_size, Design_slowExp),
  
  
  ncores = n_sim_cores,
  control = list(include_replication_index = TRUE),
  store_results = FALSE, # Mmmm... May want to change to FALSE to improve RAM
  filename = "Good_Design_Exp20",
  dirname = if(exists("results_dir")) {here::here(results_dir)} else {"/project/def-cribbie/udialter/sim_results/expanded/"})
# └───────────────────────────────────────────────┘

# --- Logging ---
#cat("Array ID =", arrayID, "\nConditions (incl.):", min(ix), "to", max(ix) )
message("Memory used (GB): ", round(sum(gc()[,2]) / 1024, 2))
sim_end <- Sys.time()
message("TOTAL JOB TIME (minutes): ", round(difftime(sim_end, sim_start, units = "mins"), 2))

