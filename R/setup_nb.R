# stuffineed3.R
# ┌── All constants & objects ───────────────────────────────────────────────────────────────────────────────┐
# Core layout 
total_cores    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 192))

  # PARALLEL MODE: spread across available cores
  brms_cores     <- 4L   # keep brms single-process
  n_brms_chains  <- 4L
  n_brms_threads <- 1L   # threads per chain
  n_sim_cores    <- total_cores %/% (brms_cores * n_brms_threads)
  stopifnot(total_cores %% (brms_cores * n_brms_threads) == 0)
  stopifnot(n_sim_cores >= 1)


# Threads / env
Sys.setenv(STAN_NUM_THREADS = n_brms_threads,
           OPENBLAS_NUM_THREADS = 1L,
           BLAS_NUM_THREADS     = 1L,
           MKL_NUM_THREADS = 1L)

message(sprintf("Layout: sims=%d | chains=%d (seq) | threads=%d ⇒ CPUs=%d",
                n_sim_cores, n_brms_chains, n_brms_threads,
                n_sim_cores * brms_cores * n_brms_threads))


idle <- total_cores - n_sim_cores * brms_cores * n_brms_threads
if (idle > 0) message(sprintf("%d core(s) will idle (integer division).", idle))



# persistent location
cache_dir <- "/home/udialter/brms_cache"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
# Scratch (per-node) tmpdir for general R tempfiles
tmp_path <- Sys.getenv("TMPDIR", unset = tempdir())
Sys.setenv(TMPDIR = tmp_path, TMP = tmp_path, TEMP = tmp_path)

exe_dir <- "/home/udialter/brms_exec"
dir.create(exe_dir, showWarnings = FALSE, recursive = TRUE)

options(cmdstanr_write_stan_file_dir = cache_dir)
options(cmdstanr_write_dir = cache_dir)          # both for safety
Sys.setenv(CMDSTANR_WRITE_STAN_FILE_DIR = cache_dir)
Sys.setenv(CMDSTANR_WRITE_DIR = cache_dir)

options(brms.cache = cache_dir)
Sys.setenv(BRMS_CACHE = cache_dir)
options(mc.cores = n_sim_cores)


results_dir <- "/scratch/udialter/sim_results/"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# Expose for your wrappers
assign("n_brms_chains",  n_brms_chains,  envir = .GlobalEnv)
assign("n_brms_threads", n_brms_threads, envir = .GlobalEnv)

# └───────────────────────────────────────────────────────────────────────────────────────────────────┘


# ┌── Load libraries & CmdStan path ───────────────────────────┐
library(Matrix)
library(here)
library(tidyverse)
library(lme4)
library(assertthat)
library(snow)
library(SimDesign)
#rebuild_cmdstan()
library(cmdstanr)
library(brms)
library(mixedup)


# Use prebuilt CmdStan (from sbatch), with a safe fallback
cs <- Sys.getenv("CMDSTAN")
if (!nzchar(cs)) cs <- "/home/udialter/cmdstan-2.36.0"
cmdstanr::set_cmdstan_path(cs)

options(brms.backend = "cmdstanr")
options(brms.backend.verbose = TRUE)

# (tmp_path + cmdstanr_write_dir already set above)
message("CmdStan path: ", cmdstanr::cmdstan_path())
message("Write dir:    ", getOption("cmdstanr_write_dir"))

if (!dir.exists(cmdstanr::cmdstan_path())) {
  stop("CmdStan not found at: ", cmdstanr::cmdstan_path(),
       "\nSet CMDSTAN in sbatch or adjust fallback path.")
}

# └────────────────────────────────────────────────────────────┘


# ┌── Define Local Functions ───────────────────────────┐

# array2row <- function(arrayID, k=1L) {
#   start <- 1L + (arrayID - 1L)*k
#   end   <- min(start + k - 1L, nrow(Design))
#   start:end
# }

# General function for simulating data with two groups and regression error variance
simulate_cluster_data <- function(n, n_MT1, intercept, MTeffect, sigma) {
  # Input validation
  if (!is.numeric(n) | !is.numeric(n_MT1) | !is.numeric(intercept) | !is.numeric(MTeffect) | !is.numeric(sigma))stop("All inputs must be numeric.")
  if (n_MT1 >= n) stop("n_MT1 must be less than n.")
  #if (n < 2 * min_group_size) stop("n must be at least twice the minimum group size.")
  #if (n_MT1 < min_group_size | (n - n_MT1) < min_group_size) stop("Each group must have at least min_group_size participants.")
  # Define group sizes
  group_sizes <- c(n - n_MT1, n_MT1)
  # Create a vector of group membership (MT)
  MT <- rep(0:1, times = group_sizes)
  # Calculate group effect
  MTeffect <- MTeffect * MT
  # Generate outcome variable with random errors
  ment_health <- intercept + MTeffect + rnorm(n, mean = 0, sd = sigma)
  # Combine into a data frame
  data <- data.frame(MT = MT, ment_health = ment_health)
  return(data)
}


# Defining local function for collecting t-test results quickly
t_test_func <- function(subset_data) {
  result <- t.test(ment_health ~ MT, data = subset_data)
  c(p.value = result$p.value,
    mean_diff = result$estimate[2] - result$estimate[1],
    up_ci = -result$conf.int[1],
    low_ci = -result$conf.int[2])
}


# Function to fit the model and handle errors/warnings
fit_lmer_with_control <- function(data, optimizer, maxfun = 1e5) {
  tryCatch({
    mod <- lme4::lmer(ment_health ~ MT + (1 + MT | cluster), data = data,
                      control = lmerControl(optimizer = optimizer, 
                                            optCtrl = list(maxfun = maxfun)))
    # Check for singular fit and log a message if detected
    if (isSingular(mod, tol = 1e-4)) {
      message(paste("Singular fit detected with optimizer:", optimizer))
    }
    return(mod)
  }, warning = function(w) {
    message(paste("Warning with optimizer:", optimizer, "-", conditionMessage(w)))
    return(NULL)  # Return NULL if a warning occurs
  }, error = function(e) {
    message(paste("Error with optimizer:", optimizer, "-", conditionMessage(e)))
    return(NULL)  # Return NULL if an error occurs
  })
}



process_CI_results <- function(mod, group = "cluster", effect = "MT") {
 
  if (is.null(mod) || !inherits(mod, c("brmsfit", "lmerMod"))) {
    stop("process_CI_results(): object is not a valid brmsfit or lmerMod")
  }
  
  
  
   if (inherits(mod, "brmsfit")) {
    ## Posterior quantile intervals
    co <- coef(mod)[[group]] # 3D array: [groups, effects, stats]
    coMT <- co[,,effect] # MT only, no intercept
    out <- data.frame(cluster = 1:nrow(coMT),
                      value   = coMT[,"Estimate"], 
                      se      = coMT[,"Est.Error"], 
                      LB      = coMT[,"Q2.5"], 
                      UB      =  coMT[,"Q97.5"])
  } else {
    ## Wald intervals from lmer
    r <- mixedup::extract_random_coefs(mod)
    r <- r[r$effect == effect, ]
    out <- data.frame(
      cluster = r$group,
      value   = r$value,
      se      = r$se,
      LB      = r$lower_2.5,
      UB      = r$upper_97.5
    )
  }
  out
}


fit_brms_with_control <- function(data,
                                  model_tag,
                                  formula = ment_health ~ MT + (1 + MT | cluster),
                                  use_template = TRUE,
                                  exe_dir,
                                  cores = fixed_objects$brms_cores,
                                  chains = fixed_objects$n_brms_chains,
                                  iter = fixed_objects$iter,
                                  warmup = fixed_objects$warmup,
                                  ...) {
  stopifnot(is.character(model_tag), nzchar(model_tag))
  
  fit <- NULL
  
  if (isTRUE(use_template)) {
    # reload the persistent template
    
    tpl_base <- file.path(exe_dir, paste0("tpl_", model_tag))
    tpl_rds  <- paste0(tpl_base, ".rds")
    if (!file.exists(tpl_rds)) stop("Template not found: ", tpl_rds)
    
    tpl <- brm(
      formula = formula,
      backend = "cmdstanr",
      file    = tpl_base,
      file_refit = "never"
    )
    
    if (!inherits(tpl, "brmsfit")) {
      message("template is not brmsfit: ", model_tag)
    } else {
      message("Using template for ", model_tag)
      fit <- tryCatch(
        update(
          tpl,
          newdata    = data,
          recompile  = FALSE,
          file_refit = "never",
          cores      = cores,
          chains     = chains,
          iter       = iter,
          warmup     = warmup,
          refresh  = 50, silent = 2,
          ...
        ),
        error = function(e) {
          message("update() failed: ", conditionMessage(e))
          NULL
        }
      )
    }
  }
  
  fit
}



has_bad_rhat <- function(model, monitor_params = c("b_MT", "sd_cluster__MT")) {
  if (is.null(model)) return(TRUE)
  rhat_vals <- tryCatch(brms::rhat(model), error = function(e) NA_real_)
  monitored <- rhat_vals[names(rhat_vals) %in% monitor_params]
  any(is.na(monitored) | monitored > 1.1)
}



fit_brms_model <- function(data,
                           model_tag, 
                           exe_dir,
                           ...) {
  
  
  bmod <- fit_brms_with_control(data=data, 
                                model_tag=model_tag, 
                                exe_dir = exe_dir, 
                                ...)
  
  if (is.null(bmod) || has_bad_rhat(bmod)) {
    
    bmod <- NULL
    stop("All model fitting attempts failed or had poor convergence.")
    
  }
  return(bmod)
}
# └───────────────────────────────────────────────────────┘

# ┌── Define Priors & Stan code ───────────────────────────┐

# Spike and slab (equal weights 0.5 / 0.5)
stanvars_ss <- stanvar(scode = "
  real spike_slab_prior(real x) {
    real lp_spike = normal_lpdf(x | 0, 0.1);
    real lp_slab  = normal_lpdf(x | 0, 5);
    return log_mix(0.5, lp_spike, lp_slab);
  }
", block = "functions") +
  stanvar(scode = "target += spike_slab_prior(b[1]);", block = "model")




prior_list <- list(
  
  brms_default     = list(prior = NULL, stanvars = NULL, model_tag="default"),
  very_info_nill   = list(prior = prior(normal(0, 0.5), class = "b", coef = "MT"), stanvars = NULL, model_tag= "VIN"),
  weak_info_nill   = list(prior = prior(normal(0, 3), class = "b", coef = "MT"), stanvars = NULL , model_tag= "WIN"),

  
  hs               = list(prior = prior(horseshoe(df = 1, par_ratio = 0.5), class = "b"), stanvars = NULL, model_tag= "HS"),     # Note no coef = "MT"
  spike_slab       = list(prior = prior("", class = "b", coef = "MT"), stanvars = stanvars_ss, model_tag= "SS"),
  
  
  tau_sd_MT_t      = list(prior = prior(student_t(3, 0, 1), class = "sd", group = "cluster", coef = "MT"), stanvars = NULL, model_tag= "tau_t"),
  tau_sd_MT_exp    = list(prior = prior(exponential(2),     class = "sd", group = "cluster", coef = "MT"), stanvars = NULL, model_tag= "tau_exp"),
  tau_sd_MT_log    = list(prior = prior(lognormal(log(2.5), 0.20), class = "sd", group = "cluster", coef = "MT"), stanvars = NULL, model_tag= "tau_log"),
  tau_sd_MT_gamma  = list(prior = prior(gamma(2.5, 1),                 class = "sd", group = "cluster", coef = "MT"), stanvars = NULL, model_tag= "tau_gam"))
# └───────────────────────────────────────────────────────┘



# ┌── Bundle everything into fixed_objects ───────────────────────────┐
fo <- fixed_objects <- list(
  
  total_cores              = total_cores,
  n_sim_cores              = n_sim_cores,
  brms_cores               = brms_cores,
  n_brms_chains            = n_brms_chains,
  n_brms_threads           = n_brms_threads,
  
  exe_dir                  = exe_dir,
  results_dir              = results_dir,
  
 # array2row               = array2row,
  simulate_cluster_data   = simulate_cluster_data,
  t_test_func             = t_test_func,
  fit_lmer_with_control   = fit_lmer_with_control,
  process_CI_results      = process_CI_results, 
  fit_brms_with_control   = fit_brms_with_control,
  has_bad_rhat            = has_bad_rhat, 
  fit_brms_model          = fit_brms_model,
  stanvars_ss             = stanvars_ss,
  
  prior_list              = prior_list,
  
  adapt_delta             = 0.8, 
  max_treedepth           = 10L,
  warmup                  = 500L,
  iter                    = 1500L,
  
  chunk_size = 1L,
  intercept = 0L,
  alpha = 0.05,
  NAflag = -999311999)
# └───────────────────────────────────────────────────────┘

