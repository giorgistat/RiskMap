<<<<<<< HEAD

=======
>>>>>>> 2ec0b26 (Completion of coding for the "dsgm" function)
##' @title Convert penalty list to TMB format
##' @keywords internal
convert_penalty_to_tmb <- function(penalty) {

  # Default values (no penalty)
  tmb_penalty <- list(
    use_alpha_penalty = 0,
    alpha_penalty_type = 1,
    alpha_param1 = 2,
    alpha_param2 = 2,
    use_gamma_penalty = 0,
    gamma_penalty_type = 1,
    gamma_param1 = 2,
    gamma_param2 = 1
  )

  if(!is.null(penalty)) {

    # =========================================================================
    # ALPHA PENALTY DETECTION
    # =========================================================================

    # Check for parameter-based penalty (NEW FORMAT)
    if(!is.null(penalty$alpha_param1) && !is.null(penalty$alpha_param2)) {
      tmb_penalty$use_alpha_penalty <- 1
      tmb_penalty$alpha_penalty_type <- 1  # Beta distribution
      tmb_penalty$alpha_param1 <- penalty$alpha_param1
      tmb_penalty$alpha_param2 <- penalty$alpha_param2
    }
    # Check for function-based penalty (OLD FORMAT)
    else if(!is.null(penalty$alpha) && is.function(penalty$alpha)) {
      tmb_penalty$use_alpha_penalty <- 1
      tmb_penalty$alpha_penalty_type <- 1

      # Try to infer Beta parameters from function
      # For Beta(a, b): -(a-1)*log(x) - (b-1)*log(1-x)
      test_val <- 0.5
      pen_val <- penalty$alpha(test_val)

      # If we have gradient function, use it to infer parameters
      if(!is.null(penalty$alpha_grad)) {
        grad_val <- penalty$alpha_grad(test_val)
        # At x=0.5: grad = -(a-1)/0.5 + (b-1)/0.5
        # If symmetric: a = b
        a_minus_1 <- -grad_val / 2
        tmb_penalty$alpha_param1 <- a_minus_1 + 1
        tmb_penalty$alpha_param2 <- a_minus_1 + 1
      } else {
        # Default to Beta(2,2)
        warning("Could not infer Beta parameters from alpha function. Using Beta(2,2).")
        tmb_penalty$alpha_param1 <- 2
        tmb_penalty$alpha_param2 <- 2
      }
    }

    # =========================================================================
    # GAMMA_W PENALTY DETECTION
    # =========================================================================

    # Check for Gamma distribution (NEW FORMAT)
    if(!is.null(penalty$gamma_type) && penalty$gamma_type == "gamma") {
      tmb_penalty$use_gamma_penalty <- 1
      tmb_penalty$gamma_penalty_type <- 1  # Gamma distribution
      tmb_penalty$gamma_param1 <- penalty$gamma_shape
      tmb_penalty$gamma_param2 <- penalty$gamma_rate
    }
    # Check for Normal distribution (NEW FORMAT)
    else if(!is.null(penalty$gamma_type) && penalty$gamma_type == "normal") {
      tmb_penalty$use_gamma_penalty <- 1
      tmb_penalty$gamma_penalty_type <- 2  # Normal distribution
      tmb_penalty$gamma_param1 <- penalty$gamma_mean
      tmb_penalty$gamma_param2 <- penalty$gamma_sd
    }
    # Check for direct Gamma parameters (ALTERNATIVE NEW FORMAT)
    else if(!is.null(penalty$gamma_shape) && !is.null(penalty$gamma_rate)) {
      tmb_penalty$use_gamma_penalty <- 1
      tmb_penalty$gamma_penalty_type <- 1  # Gamma distribution
      tmb_penalty$gamma_param1 <- penalty$gamma_shape
      tmb_penalty$gamma_param2 <- penalty$gamma_rate
    }
    # Check for direct Normal parameters (ALTERNATIVE NEW FORMAT)
    else if(!is.null(penalty$gamma_mean) && !is.null(penalty$gamma_sd)) {
      tmb_penalty$use_gamma_penalty <- 1
      tmb_penalty$gamma_penalty_type <- 2  # Normal distribution
      tmb_penalty$gamma_param1 <- penalty$gamma_mean
      tmb_penalty$gamma_param2 <- penalty$gamma_sd
    }
    # Check for function-based penalty (OLD FORMAT)
    else if(!is.null(penalty$gamma) && is.function(penalty$gamma)) {
      tmb_penalty$use_gamma_penalty <- 1

      # Default to Gamma(2,1) if can't infer
      warning("Using function-based gamma penalty. Defaulting to Gamma(2,1).")
      tmb_penalty$gamma_penalty_type <- 1
      tmb_penalty$gamma_param1 <- 2
      tmb_penalty$gamma_param2 <- 1
    }
  }

  return(tmb_penalty)
}

##' @title Fit DSGM using TMB
##' @description MCML estimation using TMB for automatic differentiation
##' @keywords internal
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats nlminb
dsgm_fit_tmb <- function(y_prev, intensity_data, D, coords, ID_coords,
                         int_mat, survey_times_data, mda_times,
                         par0, cov_offset, fix_alpha_W, fix_gamma_W,
                         penalty, S_samples_obj,
                         use_hessian_refinement = TRUE,
                         messages = TRUE) {
  # Extract data
  n <- length(y_prev)
  n_loc <- nrow(coords)
  S_samples <- S_samples_obj$S_samples
<<<<<<< HEAD
  dist_mat <- as.matrix(dist(coords))
  pos_idx <- which(y_prev == 1) - 1  # 0-indexed

=======
  pos_idx <- which(y_prev == 1) - 1  # 0-indexed

  # =============================================================================
  # OPTIMIZATION 1: SPARSE MDA REPRESENTATION
  # =============================================================================
  if(messages) message("Preprocessing sparse MDA matrix...")

  mda_sparse_idx <- which(int_mat > 0, arr.ind = TRUE)

  if(nrow(mda_sparse_idx) > 0) {
    mda_i <- as.integer(mda_sparse_idx[, 1] - 1)  # 0-indexed
    mda_j <- as.integer(mda_sparse_idx[, 2] - 1)  # 0-indexed
    mda_coverage <- as.numeric(int_mat[mda_sparse_idx])
    n_mda_pairs <- length(mda_i)
  } else {
    # No MDA exposure at all
    mda_i <- integer(0)
    mda_j <- integer(0)
    mda_coverage <- numeric(0)
    n_mda_pairs <- 0L
  }

  if(messages) {
    sparsity <- 100 * (1 - n_mda_pairs / (n * length(mda_times)))
    message(sprintf("  MDA matrix sparsity: %.1f%% (reduced from %d to %d entries)",
                    sparsity, n * length(mda_times), n_mda_pairs))
  }

  # =============================================================================
  # OPTIMIZATION 2: COMPRESS DISTANCE MATRIX
  # =============================================================================
  if(messages) message("Compressing distance matrix...")

  # Store as vector (internally dist() already does this)
  dist_vec <- as.numeric(dist(coords))

  if(messages) {
    full_size <- n_loc * n_loc * 8 / (1024^2)  # MB
    compressed_size <- length(dist_vec) * 8 / (1024^2)  # MB
    message(sprintf("  Distance matrix: %.2f MB -> %.2f MB (%.1f%% reduction)",
                    full_size, compressed_size, 100 * (1 - compressed_size/full_size)))
  }

>>>>>>> 2ec0b26 (Completion of coding for the "dsgm" function)
  tmb_penalty <- convert_penalty_to_tmb(penalty)

  # =============================================================================
  # STEP 1: COMPUTE DENOMINATOR USING TMB AT θ₀
  # =============================================================================

  # Data for denominator computation
  data_denom <- list(
    y_prev = y_prev,
    intensity_data = intensity_data,
    pos_idx = pos_idx,
    D = D,
    cov_offset = cov_offset,
    S_samples = S_samples,
    ID_coords = ID_coords - 1,
<<<<<<< HEAD
    dist_mat = dist_mat,
    survey_times = survey_times_data,
    mda_times = mda_times,
    int_mat = int_mat,
    use_alpha_penalty = 0,  # No penalty for denominator
=======
    dist_vec = dist_vec,
    n_loc = as.integer(n_loc),
    survey_times = survey_times_data,
    mda_times = mda_times,
    mda_i = mda_i,
    mda_j = mda_j,
    mda_coverage = mda_coverage,
    n_mda_pairs = as.integer(n_mda_pairs),
    use_alpha_penalty = 0,
>>>>>>> 2ec0b26 (Completion of coding for the "dsgm" function)
    alpha_penalty_type = 1,
    alpha_param1 = 2,
    alpha_param2 = 2,
    use_gamma_penalty = 0,
    gamma_penalty_type = 1,
    gamma_param1 = 2,
    gamma_param2 = 1,
<<<<<<< HEAD
    compute_denominator_only = 1,  # SPECIAL FLAG
    log_denominator_vals = numeric(nrow(S_samples))  # Dummy, not used
=======
    compute_denominator_only = 1,
    log_denominator_vals = numeric(nrow(S_samples))
>>>>>>> 2ec0b26 (Completion of coding for the "dsgm" function)
  )

  # Parameters at θ₀
  params_denom <- list(
    beta = par0$beta,
    log_k = log(par0$k),
    log_rho = log(par0$rho),
    logit_alpha = qlogis(par0$alpha_W),
    log_gamma = log(par0$gamma_W),
    log_sigma2 = log(par0$sigma2),
    log_phi = log(par0$phi)
  )

  # Build TMB object for denominator
  obj_denom <- TMB::MakeADFun(
    data = data_denom,
    parameters = params_denom,
    DLL = "RiskMap",
    silent = TRUE
  )

  # Extract denominator values
<<<<<<< HEAD
  obj_denom$fn()  # Trigger computation
=======
  obj_denom$fn()
>>>>>>> 2ec0b26 (Completion of coding for the "dsgm" function)
  log_denominator_vals <- obj_denom$report()$log_f_vals


  # =============================================================================
  # STEP 2: BUILD MAIN TMB OBJECT FOR OPTIMIZATION
  # =============================================================================

  if(messages) {
    message("Building TMB objective for optimization...")
  }

  data_list <- list(
    y_prev = y_prev,
    intensity_data = intensity_data,
    pos_idx = pos_idx,
    D = D,
    cov_offset = cov_offset,
    S_samples = S_samples,
    ID_coords = ID_coords - 1,
<<<<<<< HEAD
    dist_mat = dist_mat,
    survey_times = survey_times_data,
    mda_times = mda_times,
    int_mat = int_mat,
=======
    dist_vec = dist_vec,
    n_loc = as.integer(n_loc),
    survey_times = survey_times_data,
    mda_times = mda_times,
    mda_i = mda_i,
    mda_j = mda_j,
    mda_coverage = mda_coverage,
    n_mda_pairs = as.integer(n_mda_pairs),
>>>>>>> 2ec0b26 (Completion of coding for the "dsgm" function)
    use_alpha_penalty = tmb_penalty$use_alpha_penalty,
    alpha_penalty_type = tmb_penalty$alpha_penalty_type,
    alpha_param1 = tmb_penalty$alpha_param1,
    alpha_param2 = tmb_penalty$alpha_param2,
    use_gamma_penalty = tmb_penalty$use_gamma_penalty,
    gamma_penalty_type = tmb_penalty$gamma_penalty_type,
    gamma_param1 = tmb_penalty$gamma_param1,
    gamma_param2 = tmb_penalty$gamma_param2,
<<<<<<< HEAD
    compute_denominator_only = 0,  # Normal mode
    log_denominator_vals = log_denominator_vals  # Use computed values
=======
    compute_denominator_only = 0,
    log_denominator_vals = log_denominator_vals
>>>>>>> 2ec0b26 (Completion of coding for the "dsgm" function)
  )

  parameters <- list(
    beta = par0$beta,
    log_k = log(par0$k),
    log_rho = log(par0$rho),
    logit_alpha = qlogis(par0$alpha_W),
    log_gamma = log(par0$gamma_W),
    log_sigma2 = log(par0$sigma2),
    log_phi = log(par0$phi)
  )

  # Map for fixed parameters
  map_list <- list()
  if(!is.null(fix_alpha_W)) {
    parameters$logit_alpha <- qlogis(fix_alpha_W)
    map_list$logit_alpha <- factor(NA)
  }
  if(!is.null(fix_gamma_W)) {
    parameters$log_gamma <- log(fix_gamma_W)
    map_list$log_gamma <- factor(NA)
  }

  # Build main object
  obj <- TMB::MakeADFun(
    data = data_list,
    parameters = parameters,
    map = if(length(map_list) > 0) map_list else NULL,
    DLL = "RiskMap",
    silent = !messages
  )

  # =============================================================================
  # SANITY CHECK
  # =============================================================================

  obj_at_par0 <- obj$fn(obj$par)

  # Expected penalty
  expected_penalty <- 0
  if(tmb_penalty$use_alpha_penalty && tmb_penalty$alpha_penalty_type == 1) {
    expected_penalty <- expected_penalty -
      (tmb_penalty$alpha_param1 - 1) * log(par0$alpha_W) -
      (tmb_penalty$alpha_param2 - 1) * log(1 - par0$alpha_W)
  }
  if(tmb_penalty$use_gamma_penalty) {
    if(tmb_penalty$gamma_penalty_type == 1) {
      expected_penalty <- expected_penalty -
        (tmb_penalty$gamma_param1 - 1) * log(par0$gamma_W) +
        tmb_penalty$gamma_param2 * par0$gamma_W
    } else if(tmb_penalty$gamma_penalty_type == 2) {
      diff <- par0$gamma_W - tmb_penalty$gamma_param1
      expected_penalty <- expected_penalty +
        0.5 * diff^2 / (tmb_penalty$gamma_param2^2)
    }
  }

  if(abs(obj_at_par0 - expected_penalty) > 0.1) {
    stop("SANITY CHECK FAILED!")
  }

  # Stage 1: Gradient only
  if(messages) message("Stage 1: Optimization with gradient...")

<<<<<<< HEAD
  opt1 <- nlminb(obj$par, obj$fn, obj$gr,
                 control = list(eval.max = 500, iter.max = 250,
                                trace = ifelse(messages, 1, 0)))
=======
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(eval.max = 500, iter.max = 250,
                               trace = ifelse(messages, 1, 0)))
>>>>>>> 2ec0b26 (Completion of coding for the "dsgm" function)

  # Standard errors
  if(messages) message("Computing standard errors...")
  sdr <- TMB::sdreport(obj)

  par_est <- summary(sdr, "report")

  result <- list(
    params = list(
      beta = obj$env$parList()$beta,
      k = par_est["k", "Estimate"],
      rho = par_est["rho", "Estimate"],
      alpha_W = par_est["alpha_W", "Estimate"],
      gamma_W = par_est["gamma_W", "Estimate"],
      sigma2 = par_est["sigma2", "Estimate"],
      phi = par_est["phi", "Estimate"]
    ),
    params_se = list(
      beta = summary(sdr, "fixed")[grepl("beta", rownames(summary(sdr, "fixed"))), "Std. Error"],
      k = par_est["k", "Std. Error"],
      rho = par_est["rho", "Std. Error"],
      alpha_W = par_est["alpha_W", "Std. Error"],
      gamma_W = par_est["gamma_W", "Std. Error"],
      sigma2 = par_est["sigma2", "Std. Error"],
      phi = par_est["phi", "Std. Error"]
    ),
    convergence = opt$convergence,
    log_likelihood = -opt$objective,
    message = opt$message,
    iterations = opt$iterations,
    evaluations = opt$evaluations,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr,
    posterior_samples = S_samples_obj
  )

  return(result)
}
