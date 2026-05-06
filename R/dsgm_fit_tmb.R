##' @title Compile and load the dsgm_mdiag TMB template on first use
##' @keywords internal
.load_dsgm_mdiag <- function(messages = TRUE) {
  if (isTRUE(.dsgm_cache$mdiag_loaded)) return(invisible(NULL))

  cpp_file <- system.file("tmb/dsgm_mdiag.cpp", package = "RiskMap")
  if (!file.exists(cpp_file))
    stop("TMB template not found: inst/tmb/dsgm_mdiag.cpp")

  if (messages)
    message("Compiling dsgm_mdiag TMB template (first use this session)...")

  TMB::compile(cpp_file, silent = !messages)
  dyn.load(TMB::dynlib(tools::file_path_sans_ext(cpp_file)))

  .dsgm_cache$mdiag_loaded <- TRUE
  if (messages) message("dsgm_mdiag compiled and loaded.")
  invisible(NULL)
}


##' @title Convert penalty list to TMB format
##' @keywords internal
convert_penalty_to_tmb <- function(penalty) {

  tmb_penalty <- list(
    use_alpha_penalty  = 0,
    alpha_penalty_type = 1,
    alpha_param1       = 2,
    alpha_param2       = 2,
    use_gamma_penalty  = 0,
    gamma_penalty_type = 1,
    gamma_param1       = 2,
    gamma_param2       = 1,
    use_rho_penalty    = 0,   # <- missing
    rho_penalty_type   = 3,
    rho_param1         = 0,
    rho_param2         = 1
  )

  if (is.null(penalty)) return(tmb_penalty)

  # ===========================================================================
  # ALPHA PENALTY
  # ===========================================================================
  if (!is.null(penalty$alpha_param1) && !is.null(penalty$alpha_param2)) {
    tmb_penalty$use_alpha_penalty  <- 1
    tmb_penalty$alpha_penalty_type <- 1
    tmb_penalty$alpha_param1       <- penalty$alpha_param1
    tmb_penalty$alpha_param2       <- penalty$alpha_param2

  } else if (!is.null(penalty$alpha) && is.function(penalty$alpha)) {
    tmb_penalty$use_alpha_penalty  <- 1
    tmb_penalty$alpha_penalty_type <- 1
    if (!is.null(penalty$alpha_grad)) {
      grad_val  <- penalty$alpha_grad(0.5)
      a_minus_1 <- -grad_val / 2
      tmb_penalty$alpha_param1 <- a_minus_1 + 1
      tmb_penalty$alpha_param2 <- a_minus_1 + 1
    } else {
      warning("Could not infer Beta parameters from alpha function. Using Beta(2,2).")
      tmb_penalty$alpha_param1 <- 2
      tmb_penalty$alpha_param2 <- 2
    }
  }

  # ===========================================================================
  # GAMMA_W PENALTY
  # ===========================================================================
  if (!is.null(penalty$gamma_type)) {

    tmb_penalty$use_gamma_penalty <- 1

    if (penalty$gamma_type == "gamma") {
      tmb_penalty$gamma_penalty_type <- 1
      tmb_penalty$gamma_param1       <- penalty$gamma_shape
      tmb_penalty$gamma_param2       <- penalty$gamma_rate

    } else if (penalty$gamma_type == "normal") {
      tmb_penalty$gamma_penalty_type <- 2
      tmb_penalty$gamma_param1       <- penalty$gamma_mean
      tmb_penalty$gamma_param2       <- penalty$gamma_sd

    } else if (penalty$gamma_type == "lognormal") {
      # Log-Normal: Normal(mean, sd) on log(gamma_W)
      # gradient at theta_0 is O(1) vs O(gamma_W) for Gamma prior
      tmb_penalty$gamma_penalty_type <- 3
      tmb_penalty$gamma_param1       <- penalty$gamma_mean  # mu on log scale
      tmb_penalty$gamma_param2       <- penalty$gamma_sd    # sd on log scale

    } else {
      warning(sprintf("Unknown gamma_type '%s'. No gamma penalty applied.", penalty$gamma_type))
      tmb_penalty$use_gamma_penalty <- 0
    }

  } else if (!is.null(penalty$gamma_shape) && !is.null(penalty$gamma_rate)) {
    tmb_penalty$use_gamma_penalty  <- 1
    tmb_penalty$gamma_penalty_type <- 1
    tmb_penalty$gamma_param1       <- penalty$gamma_shape
    tmb_penalty$gamma_param2       <- penalty$gamma_rate

  } else if (!is.null(penalty$gamma_mean) && !is.null(penalty$gamma_sd)) {
    warning("gamma_type not specified; defaulting to 'lognormal'. Set gamma_type explicitly.")
    tmb_penalty$use_gamma_penalty  <- 1
    tmb_penalty$gamma_penalty_type <- 3        # lognormal — consistent with convention
    tmb_penalty$gamma_param1       <- penalty$gamma_mean
    tmb_penalty$gamma_param2       <- penalty$gamma_sd
  } else if (!is.null(penalty$gamma) && is.function(penalty$gamma)) {
    warning("Using function-based gamma penalty. Defaulting to Gamma(2,1).")
    tmb_penalty$use_gamma_penalty  <- 1
    tmb_penalty$gamma_penalty_type <- 1
    tmb_penalty$gamma_param1       <- 2
    tmb_penalty$gamma_param2       <- 1
  }

  # ===========================================================================
  # RHO PENALTY
  # ===========================================================================

  if (!is.null(penalty$rho_type)) {
    tmb_penalty$use_rho_penalty <- 1
    if (penalty$rho_type == "gamma") {
      tmb_penalty$rho_penalty_type <- 1
      tmb_penalty$rho_param1       <- penalty$rho_shape
      tmb_penalty$rho_param2       <- penalty$rho_rate
    } else if (penalty$rho_type == "normal") {
      tmb_penalty$rho_penalty_type <- 2
      tmb_penalty$rho_param1       <- penalty$rho_mean
      tmb_penalty$rho_param2       <- penalty$rho_sd
    } else if (penalty$rho_type == "lognormal") {
      tmb_penalty$rho_penalty_type <- 3
      tmb_penalty$rho_param1       <- penalty$rho_mean  # mean of log(rho)
      tmb_penalty$rho_param2       <- penalty$rho_sd    # sd of log(rho)
    }
  } else if (!is.null(penalty$rho_mean) && !is.null(penalty$rho_sd)) {
    # Convenience: infer log-normal if mean/sd given without type
    tmb_penalty$use_rho_penalty  <- 1
    tmb_penalty$rho_penalty_type <- 3
    tmb_penalty$rho_param1       <- penalty$rho_mean
    tmb_penalty$rho_param2       <- penalty$rho_sd
  }

  return(tmb_penalty)
}


##' @title Fit DSGM using TMB
##' @description MCML estimation using TMB for automatic differentiation.
##' @param intensity_family Integer; 0 = shifted Gamma (default), 1 = zero-truncated NegBin.
##' @keywords internal
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats nlminb
dsgm_fit_tmb <- function(y_prev            = NULL,
                         intensity_data    = NULL,
                         D,
                         coords,
                         ID_coords,
                         int_mat           = NULL,
                         survey_times_data = NULL,
                         mda_times         = NULL,
                         par0,
                         cov_offset,
                         fix_alpha_W       = NULL,
                         fix_gamma_W       = NULL,
                         penalty           = NULL,
                         S_samples_obj,
                         # lf_mdiag-specific
                         model             = "sth",
                         y_counts          = NULL,
                         units_m           = NULL,
                         which_diag        = NULL,
                         gamma_sens        = 0.97,
                         fix_k             = NULL,
                         fix_tau2          = NULL,
                         use_mda           = TRUE,
                         intensity_family  = 0L,    # 0 = shifted Gamma, 1 = zero-trunc NegBin
                         use_hessian_refinement = TRUE,
                         messages          = TRUE) {

  model <- match.arg(model, c("sth", "lf_mdiag"))

  if (model == "lf_mdiag") {
    return(.dsgm_fit_tmb_lf_mdiag(
      y_counts          = y_counts,
      units_m           = units_m,
      which_diag        = which_diag,
      D                 = D,
      coords            = coords,
      ID_coords         = ID_coords,
      cov_offset        = cov_offset,
      gamma_sens        = gamma_sens,
      fix_k             = fix_k,
      fix_tau2          = fix_tau2,
      use_mda           = use_mda,
      int_mat           = int_mat,
      survey_times_data = survey_times_data,
      mda_times         = mda_times,
      fix_alpha_W       = fix_alpha_W,
      fix_gamma_W       = fix_gamma_W,
      penalty           = penalty,
      par0              = par0,
      S_samples_obj     = S_samples_obj,
      messages          = messages
    ))
  }

  # ---------------------------------------------------------------------------
  # STH branch
  # ---------------------------------------------------------------------------
  n         <- length(y_prev)
  n_loc     <- nrow(coords)
  S_samples <- S_samples_obj$S_samples
  pos_idx   <- which(y_prev == 1) - 1L   # 0-indexed

  if (messages) message("Preprocessing sparse MDA matrix...")
  mda_sparse_idx <- which(int_mat > 0, arr.ind = TRUE)
  if (nrow(mda_sparse_idx) > 0) {
    mda_i        <- as.integer(mda_sparse_idx[, 1] - 1)
    mda_j        <- as.integer(mda_sparse_idx[, 2] - 1)
    mda_coverage <- as.numeric(int_mat[mda_sparse_idx])
    n_mda_pairs  <- length(mda_i)
  } else {
    mda_i <- integer(0); mda_j <- integer(0)
    mda_coverage <- numeric(0); n_mda_pairs <- 0L
  }
  if (messages) {
    sparsity <- 100 * (1 - n_mda_pairs / (n * length(mda_times)))
    message(sprintf("  MDA matrix sparsity: %.1f%% (reduced from %d to %d entries)",
                    sparsity, n * length(mda_times), n_mda_pairs))
  }

  if (messages) message("Compressing distance matrix...")
  dist_vec <- as.numeric(dist(coords))
  if (messages) {
    full_size       <- n_loc * n_loc * 8 / (1024^2)
    compressed_size <- length(dist_vec) * 8 / (1024^2)
    message(sprintf("  Distance matrix: %.2f MB -> %.2f MB (%.1f%% reduction)",
                    full_size, compressed_size,
                    100 * (1 - compressed_size / full_size)))
  }

  tmb_penalty <- convert_penalty_to_tmb(penalty)

  # Helper: build TMB data list, varying only the IS fields
  make_data <- function(compute_denom, log_denom_vals) {
    list(
      y_prev                   = y_prev,
      intensity_data           = intensity_data,
      pos_idx                  = pos_idx,
      D                        = D,
      cov_offset               = cov_offset,
      S_samples                = S_samples,
      ID_coords                = ID_coords - 1L,
      dist_vec                 = dist_vec,
      n_loc                    = as.integer(n_loc),
      survey_times             = survey_times_data,
      mda_times                = mda_times,
      mda_i                    = mda_i,
      mda_j                    = mda_j,
      mda_coverage             = mda_coverage,
      n_mda_pairs              = as.integer(n_mda_pairs),
      use_alpha_penalty        = tmb_penalty$use_alpha_penalty,
      alpha_penalty_type       = tmb_penalty$alpha_penalty_type,
      alpha_param1             = tmb_penalty$alpha_param1,
      alpha_param2             = tmb_penalty$alpha_param2,
      use_gamma_penalty        = tmb_penalty$use_gamma_penalty,
      gamma_penalty_type       = tmb_penalty$gamma_penalty_type,
      gamma_param1             = tmb_penalty$gamma_param1,
      gamma_param2             = tmb_penalty$gamma_param2,
      use_rho_penalty          = tmb_penalty$use_rho_penalty,
      rho_penalty_type         = tmb_penalty$rho_penalty_type,
      rho_param1               = tmb_penalty$rho_param1,
      rho_param2               = tmb_penalty$rho_param2,
      compute_denominator_only = as.integer(compute_denom),
      log_denominator_vals     = log_denom_vals,
      intensity_family         = as.integer(intensity_family)
    )
  }

  parameters <- list(
    beta        = par0$beta,
    log_k       = log(par0$k),
    log_rho     = log(par0$rho),
    logit_alpha = qlogis(par0$alpha_W),
    log_gamma   = log(par0$gamma_W),
    log_sigma2  = log(par0$sigma2),
    log_phi     = log(par0$phi)
  )

  map_list <- list()
  if (!use_mda) {
    map_list$logit_alpha <- factor(NA)
    map_list$log_gamma   <- factor(NA)
  } else {
    if (!is.null(fix_alpha_W)) {
      parameters$logit_alpha <- qlogis(fix_alpha_W)
      map_list$logit_alpha   <- factor(NA)
    }
    if (!is.null(fix_gamma_W)) {
      parameters$log_gamma <- log(fix_gamma_W)
      map_list$log_gamma   <- factor(NA)
    }
  }
  tmb_map <- if (length(map_list) > 0) map_list else NULL

  # --- Pass 1: denominator at theta_0 (no penalties) ---
  obj_d <- TMB::MakeADFun(
    data       = make_data(1L, numeric(nrow(S_samples))),
    parameters = parameters,
    map        = tmb_map,
    DLL        = "RiskMap",
    silent     = TRUE
  )
  obj_d$fn()
  log_denominator_vals <- obj_d$report()$log_f_vals

  # --- Pass 2: MCML objective ---
  if (messages) message("Building TMB objective for optimization...")
  obj <- TMB::MakeADFun(
    data       = make_data(0L, log_denominator_vals),
    parameters = parameters,
    map        = tmb_map,
    DLL        = "RiskMap",
    silent     = !messages
  )

  # --- Sanity check: NLL at theta_0 should equal penalty only ---
  obj_at_par0      <- obj$fn(obj$par)
  expected_penalty <- 0
  if (tmb_penalty$use_alpha_penalty == 1 && tmb_penalty$alpha_penalty_type == 1) {
    expected_penalty <- expected_penalty -
      (tmb_penalty$alpha_param1 - 1) * log(par0$alpha_W) -
      (tmb_penalty$alpha_param2 - 1) * log(1 - par0$alpha_W)
  }
  if (tmb_penalty$use_gamma_penalty == 1) {
    if (tmb_penalty$gamma_penalty_type == 1) {
      expected_penalty <- expected_penalty -
        (tmb_penalty$gamma_param1 - 1) * log(par0$gamma_W) +
        tmb_penalty$gamma_param2 * par0$gamma_W
    } else if (tmb_penalty$gamma_penalty_type == 2) {
      d <- par0$gamma_W - tmb_penalty$gamma_param1
      expected_penalty <- expected_penalty + 0.5 * d^2 / tmb_penalty$gamma_param2^2
    } else if (tmb_penalty$gamma_penalty_type == 3) {
      # Log-Normal: penalty lives on log(gamma_W)
      d <- log(par0$gamma_W) - tmb_penalty$gamma_param1
      expected_penalty <- expected_penalty + 0.5 * d^2 / tmb_penalty$gamma_param2^2
    }
  }
  if (tmb_penalty$use_rho_penalty == 1) {
    if (tmb_penalty$rho_penalty_type == 1) {
      expected_penalty <- expected_penalty -
        (tmb_penalty$rho_param1 - 1) * log(par0$rho) +
        tmb_penalty$rho_param2 * par0$rho
    } else if (tmb_penalty$rho_penalty_type == 2) {
      d <- par0$rho - tmb_penalty$rho_param1
      expected_penalty <- expected_penalty + 0.5 * d^2 / tmb_penalty$rho_param2^2
    } else if (tmb_penalty$rho_penalty_type == 3) {
      d <- log(par0$rho) - tmb_penalty$rho_param1
      expected_penalty <- expected_penalty + 0.5 * d^2 / tmb_penalty$rho_param2^2
    }
  }
  if (abs(obj_at_par0 - expected_penalty) > 0.1)
    stop("SANITY CHECK FAILED: NLL at theta_0 != penalty. Check denominator computation.")

  if (messages) message("Optimising...")
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(eval.max = 500,
                               iter.max = 250,
                               step.min = 1e-4,
                               step.max = 0.5,
                               trace    = ifelse(messages, 1, 0)))

  if (messages) message("Computing standard errors...")
  sdr     <- TMB::sdreport(obj)
  par_est <- summary(sdr, "report")

  result <- list(
    params = c(
      list(
        beta   = obj$env$parList()$beta,
        k      = par_est["k",      "Estimate"],
        rho    = par_est["rho",    "Estimate"],
        sigma2 = par_est["sigma2", "Estimate"],
        phi    = par_est["phi",    "Estimate"]
      ),
      if (use_mda) list(
        alpha_W = par_est["alpha_W", "Estimate"],
        gamma_W = par_est["gamma_W", "Estimate"]
      )
    ),
    params_se = c(
      list(
        beta   = summary(sdr, "fixed")[grepl("beta", rownames(summary(sdr, "fixed"))), "Std. Error"],
        k      = par_est["k",      "Std. Error"],
        rho    = par_est["rho",    "Std. Error"],
        sigma2 = par_est["sigma2", "Std. Error"],
        phi    = par_est["phi",    "Std. Error"]
      ),
      if (use_mda) list(
        alpha_W = par_est["alpha_W", "Std. Error"],
        gamma_W = par_est["gamma_W", "Std. Error"]
      )
    ),
    convergence      = opt$convergence,
    log_likelihood   = -opt$objective,
    message          = opt$message,
    iterations       = opt$iterations,
    evaluations      = opt$evaluations,
    tmb_obj          = obj,
    tmb_opt          = opt,
    tmb_sdr          = sdr,
    posterior_samples = S_samples_obj
  )
  return(result)
}


# =============================================================================
# Internal: TMB-MCML engine for the LF multi-diagnostic model
# =============================================================================

##' @keywords internal
.dsgm_fit_tmb_lf_mdiag <- function(y_counts, units_m, which_diag, D, coords,
                                   ID_coords, cov_offset, gamma_sens,
                                   fix_k, fix_tau2 = NULL, use_mda,
                                   int_mat, survey_times_data, mda_times,
                                   fix_alpha_W, fix_gamma_W, penalty,
                                   par0, S_samples_obj, messages) {

  .load_dsgm_mdiag(messages = messages)

  n         <- length(y_counts)
  n_loc     <- nrow(coords)
  S_samples <- S_samples_obj$S_samples
  n_samples <- nrow(S_samples)
  dist_vec  <- as.numeric(dist(coords))

  tmb_penalty <- convert_penalty_to_tmb(penalty)

  if (use_mda) {
    if (messages) message("  Preprocessing MDA matrix (lf_mdiag)...")
    mda_sparse   <- which(int_mat > 0, arr.ind = TRUE)
    mda_i        <- as.integer(mda_sparse[, 1] - 1L)
    mda_j        <- as.integer(mda_sparse[, 2] - 1L)
    mda_coverage <- as.numeric(int_mat[mda_sparse])
    n_mda_pairs  <- length(mda_i)
  } else {
    mda_i <- integer(0); mda_j <- integer(0)
    mda_coverage <- numeric(0); n_mda_pairs <- 0L
    survey_times_data <- rep(0, n); mda_times <- 0
  }

  tmb_data <- list(
    y                    = as.integer(y_counts),
    units_m              = as.integer(units_m),
    is_mf                = as.integer(which_diag),
    D                    = D,
    cov_offset           = cov_offset,
    S_samples            = S_samples,
    ID_coords            = as.integer(ID_coords - 1L),
    dist_vec             = dist_vec,
    n_loc                = as.integer(n_loc),
    gamma_sens           = gamma_sens,
    fix_omega            = as.integer(!is.null(fix_k)),
    fixed_omega_val      = if (!is.null(fix_k)) as.numeric(fix_k) else 0.0,
    use_mda              = as.integer(use_mda),
    survey_times         = survey_times_data,
    mda_times            = mda_times,
    mda_i                = mda_i,
    mda_j                = mda_j,
    mda_coverage         = mda_coverage,
    n_mda_pairs          = as.integer(n_mda_pairs),
    use_alpha_W_penalty  = as.integer(tmb_penalty$use_alpha_penalty),
    alpha_W_penalty_type = as.integer(tmb_penalty$alpha_penalty_type),
    alpha_W_param1       = tmb_penalty$alpha_param1,
    alpha_W_param2       = tmb_penalty$alpha_param2,
    use_gamma_W_penalty  = as.integer(tmb_penalty$use_gamma_penalty),
    gamma_W_penalty_type = as.integer(tmb_penalty$gamma_penalty_type),
    gamma_W_param1       = tmb_penalty$gamma_param1,
    gamma_W_param2       = tmb_penalty$gamma_param2
  )

  tmb_params <- list(
    beta          = par0$beta,
    log_sigma2    = log(par0$sigma2),
    log_phi       = log(par0$phi),
    log_nu2       = log(max(if (!is.null(fix_tau2)) fix_tau2 else par0$tau2,
                            1e-10) / par0$sigma2),
    log_omega     = log(if (!is.null(fix_k)) fix_k else par0$k),
    log_alpha     = log(par0$rho),
    logit_alpha_W = if (use_mda) qlogis(par0$alpha_W) else 0.0,
    log_gamma_W   = if (use_mda) log(par0$gamma_W)    else 0.0
  )

  map_list <- list()
  if (!is.null(fix_tau2)) map_list$log_nu2   <- factor(NA)
  if (!is.null(fix_k))    map_list$log_omega <- factor(NA)
  if (!use_mda) {
    map_list$logit_alpha_W <- factor(NA)
    map_list$log_gamma_W   <- factor(NA)
  } else {
    if (!is.null(fix_alpha_W)) {
      tmb_params$logit_alpha_W <- qlogis(fix_alpha_W)
      map_list$logit_alpha_W   <- factor(NA)
    }
    if (!is.null(fix_gamma_W)) {
      tmb_params$log_gamma_W <- log(fix_gamma_W)
      map_list$log_gamma_W   <- factor(NA)
    }
  }
  make_map <- function() if (length(map_list) > 0) map_list else NULL

  if (messages) message("  Computing IS denominator at theta_0 (lf_mdiag)...")
  obj_d <- TMB::MakeADFun(
    data       = c(tmb_data, list(compute_denominator_only = 1L,
                                  log_denominator_vals = numeric(n_samples))),
    parameters = tmb_params,
    map        = make_map(),
    DLL        = "dsgm_mdiag",
    silent     = TRUE
  )
  obj_d$fn()
  log_denom <- obj_d$report()$log_f_vals

  if (messages) message("  Building MCML objective (lf_mdiag)...")
  obj <- TMB::MakeADFun(
    data       = c(tmb_data, list(compute_denominator_only = 0L,
                                  log_denominator_vals = log_denom)),
    parameters = tmb_params,
    map        = make_map(),
    DLL        = "dsgm_mdiag",
    silent     = !messages
  )

  if (messages) message("  Optimising...")
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(eval.max = 1000, iter.max = 500,
                               trace = ifelse(messages, 1, 0)))

  if (messages) message("  Computing standard errors...")
  sdr     <- TMB::sdreport(obj)
  par_rep <- summary(sdr, "report")
  fix_par <- summary(sdr, "fixed")
  b_rows  <- grepl("^beta", rownames(fix_par))

  ge  <- function(nm) par_rep[nm, "Estimate"]
  gse <- function(nm) par_rep[nm, "Std. Error"]

  params <- list(
    beta   = fix_par[b_rows, "Estimate"],
    sigma2 = ge("sigma2"),
    phi    = ge("phi"),
    tau2   = if (!is.null(fix_tau2)) fix_tau2 else ge("tau2"),
    k      = if (!is.null(fix_k)) fix_k else ge("omega"),
    rho    = ge("alpha")
  )
  params_se <- list(
    beta   = fix_par[b_rows, "Std. Error"],
    sigma2 = gse("sigma2"),
    phi    = gse("phi"),
    tau2   = if (!is.null(fix_tau2)) NA else gse("tau2"),
    k      = if (!is.null(fix_k)) NA else gse("omega"),
    rho    = gse("alpha")
  )
  if (use_mda) {
    params$alpha_W    <- ge("alpha_W");  params_se$alpha_W <- gse("alpha_W")
    params$gamma_W    <- ge("gamma_W");  params_se$gamma_W <- gse("gamma_W")
  }

  list(
    params            = params,
    params_se         = params_se,
    convergence       = opt$convergence,
    log_likelihood    = -opt$objective,
    message           = opt$message,
    iterations        = opt$iterations,
    evaluations       = opt$evaluations,
    gamma_sens        = gamma_sens,
    fix_k             = fix_k,
    use_mda           = use_mda,
    tmb_obj           = obj,
    tmb_opt           = opt,
    tmb_sdr           = sdr,
    posterior_samples = S_samples_obj
  )
}
