##' @title Detect available Stan backend
##' @keywords internal
detect_stan_backend <- function(messages = TRUE) {
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    ok <- tryCatch({ cmdstanr::cmdstan_version(); TRUE }, error = function(e) FALSE)
    if (ok) {
      if (messages) message("Using cmdstanr backend")
      return("cmdstanr")
    }
  }
  if (requireNamespace("rstan", quietly = TRUE)) {
    if (messages) message("Using rstan backend (cmdstanr unavailable)")
    return("rstan")
  }
  stop("Neither rstan nor cmdstanr is available. Install cmdstanr for persistent caching across sessions.")
}


##' @title Compile and cache a Stan model
##' @keywords internal
get_stan_model <- function(model = c("sth", "lf"), backend = NULL,
                           messages = TRUE) {
  model <- match.arg(model)

  cache_model   <- paste0("stan_model_",   model)
  cache_backend <- paste0("stan_backend_", model)

  if (is.null(backend)) {
    if (!is.null(.dsgm_cache[[cache_backend]])) {
      backend <- .dsgm_cache[[cache_backend]]
    } else {
      backend <- detect_stan_backend(messages = messages)
    }
  }

  if (!is.null(.dsgm_cache[[cache_model]]) &&
      identical(.dsgm_cache[[cache_backend]], backend)) {
    if (messages) message("Using cached Stan model (", model, ")")
    return(list(model = .dsgm_cache[[cache_model]], backend = backend))
  }

  stan_file <- switch(model,
                      sth = system.file("stan/dsgm_spatial.stan", package = "RiskMap"),
                      lf  = system.file("stan/dsgm_mdiag.stan",   package = "RiskMap")
  )
  if (!file.exists(stan_file))
    stop("Stan model file not found: ", basename(stan_file))

  cache_dir <- tools::R_user_dir("RiskMap", "cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  if (backend == "rstan") {
    cached_stan <- file.path(cache_dir, basename(stan_file))
    if (!file.exists(cached_stan) ||
        file.mtime(stan_file) > file.mtime(cached_stan)) {
      file.copy(stan_file, cached_stan, overwrite = TRUE)
    }
    rds_path <- sub("\\.stan$", ".rds", cached_stan)
    if (file.exists(rds_path)) {
      if (messages) message("Loading pre-compiled rstan model (", model, ")")
      compiled <- rstan::stan_model(file = cached_stan,
                                    model_name = paste0("dsgm_", model),
                                    verbose = FALSE, auto_write = TRUE)
    } else {
      if (messages) message("Compiling Stan model '", basename(stan_file), "'...")
      compiled <- rstan::stan_model(file = cached_stan,
                                    model_name = paste0("dsgm_", model),
                                    verbose = FALSE, auto_write = TRUE)
      if (messages) message("  Saved to: ", rds_path)
    }
  } else {
    exe_path <- file.path(cache_dir,
                          paste0("dsgm_", model,
                                 if (.Platform$OS.type == "windows") ".exe" else ""))
    exe_is_fresh <- file.exists(exe_path) &&
      file.mtime(exe_path) >= file.mtime(stan_file)
    if (exe_is_fresh) {
      if (messages) message("Loading pre-compiled cmdstanr model (", model, ")")
      compiled <- cmdstanr::cmdstan_model(stan_file = stan_file,
                                          exe_file  = exe_path,
                                          compile   = FALSE)
    } else {
      if (messages) message("Compiling Stan model '", basename(stan_file), "'...")
      compiled <- cmdstanr::cmdstan_model(stan_file = stan_file,
                                          exe_file  = exe_path,
                                          compile   = TRUE)
      if (messages) message("  Saved to: ", exe_path)
    }
  }

  .dsgm_cache[[cache_model]]   <- compiled
  .dsgm_cache[[cache_backend]] <- backend
  if (messages) message("Stan model ready (", model, ")")
  return(list(model = compiled, backend = backend))
}


# =============================================================================
# Internal Stan dispatch helpers
# =============================================================================

.run_stan <- function(stan_model, backend, stan_data,
                      n_samples, n_warmup, n_chains, n_cores,
                      adapt_delta, max_treedepth, messages) {
  if (backend == "rstan") {
    rstan::sampling(
      stan_model,
      data    = stan_data,
      iter    = n_samples + n_warmup,
      warmup  = n_warmup,
      chains  = n_chains,
      cores   = n_cores,
      control = list(adapt_delta = adapt_delta,
                     max_treedepth = max_treedepth),
      refresh = ifelse(messages, max(1, (n_samples + n_warmup) %/% 10), 0),
      show_messages = messages,
      verbose       = messages
    )
  } else {
    stan_model$sample(
      data            = stan_data,
      iter_warmup     = n_warmup,
      iter_sampling   = n_samples,
      chains          = n_chains,
      parallel_chains = n_cores,
      adapt_delta     = adapt_delta,
      max_treedepth   = max_treedepth,
      refresh = ifelse(messages, max(1, (n_samples + n_warmup) %/% 10), 0),
      show_messages   = messages,
      show_exceptions = messages
    )
  }
}

.extract_S <- function(fit, backend, messages) {
  if (backend == "rstan") {
    S <- rstan::extract(fit, pars = "S")$S
    if (messages) {
      n_div <- sum(rstan::get_num_divergent(fit))
      n_mtr <- sum(rstan::get_num_max_treedepth(fit))
      message(sprintf("Extracted %d samples (%d locations) | divergent: %d | max-treedepth: %d",
                      nrow(S), ncol(S), n_div, n_mtr))
      if (n_div > 0) warning("Divergent transitions detected. Consider increasing adapt_delta.")
    }
  } else {
    S <- fit$draws("S", format = "matrix")
    if (messages) {
      diag  <- fit$diagnostic_summary()
      n_div <- sum(diag$num_divergent)
      n_mtr <- sum(diag$num_max_treedepth)
      message(sprintf("Extracted %d samples (%d locations) | divergent: %d | max-treedepth: %d",
                      nrow(S), ncol(S), n_div, n_mtr))
      if (n_div > 0) warning("Divergent transitions detected. Consider increasing adapt_delta.")
    }
  }
  S
}


# =============================================================================
# Stan samplers
# =============================================================================

##' @title Sample spatial process for the STH model
##' @param intensity_family Integer; 0 = shifted Gamma (default), 1 = zero-truncated NegBin.
##' @keywords internal
sample_spatial_process_stan <- function(y_prev,
                                        intensity_data,
                                        D,
                                        coords,
                                        ID_coords,
                                        int_mat,
                                        survey_times_data,
                                        mda_times,
                                        par,
                                        n_samples        = 1000,
                                        n_warmup         = 1000,
                                        n_chains         = 4,
                                        n_cores          = 4,
                                        adapt_delta      = 0.8,
                                        max_treedepth    = 10,
                                        intensity_family = 0L,
                                        backend          = NULL,
                                        messages         = TRUE) {

  n       <- length(y_prev)
  n_loc   <- nrow(coords)
  p       <- ncol(D)
  pos_idx <- which(y_prev == 1)
  n_pos   <- length(pos_idx)

  mda_impact <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                                   par$alpha_W, par$gamma_W, kappa = 1)

  stan_data <- list(
    n                = n,
    n_loc            = n_loc,
    n_pos            = n_pos,
    p                = p,
    y                = y_prev,
    C_pos            = intensity_data,
    C_pos_int        = as.integer(intensity_data),
    pos_idx          = pos_idx,
    ID_coords        = ID_coords,
    D_mat            = as.matrix(dist(coords)),
    eta_fixed        = as.numeric(D %*% par$beta),
    mda_impact       = mda_impact,
    k                = par$k,
    rho              = par$rho,
    sigma2           = par$sigma2,
    phi              = par$phi,
    intensity_family = as.integer(intensity_family)
  )

  mod     <- get_stan_model(model = "sth", backend = backend, messages = messages)
  backend <- mod$backend

  if (messages)
    message(sprintf("Sampling %d iter (%d warmup), %d chain(s) [sth, family=%s]...",
                    n_samples + n_warmup, n_warmup, n_chains,
                    ifelse(intensity_family == 0L, "shifted Gamma", "zero-trunc NegBin")))

  fit       <- .run_stan(mod$model, backend, stan_data,
                         n_samples, n_warmup, n_chains, n_cores,
                         adapt_delta, max_treedepth, messages)
  S_samples <- .extract_S(fit, backend, messages)

  result <- list(S_samples = S_samples, stan_fit = fit,
                 n_samples = nrow(S_samples), n_loc = n_loc,
                 coords = coords, par = par, backend = backend)
  class(result) <- "dsgm_spatial_samples"
  result
}


##' @title Sample spatial process for the LF multi-diagnostic model
##' @keywords internal
sample_spatial_process_stan_lf <- function(y_counts,
                                           units_m,
                                           which_diag,
                                           D,
                                           coords,
                                           ID_coords,
                                           par,
                                           mda_impact    = NULL,
                                           n_samples     = 1000,
                                           n_warmup      = 1000,
                                           n_chains      = 4,
                                           n_cores       = 4,
                                           adapt_delta   = 0.8,
                                           max_treedepth = 10,
                                           backend       = NULL,
                                           messages      = TRUE) {

  n     <- length(y_counts)
  n_loc <- nrow(coords)
  p     <- ncol(D)

  if (is.null(mda_impact)) mda_impact <- rep(1.0, n)
  use_mda_flag <- as.integer(!all(mda_impact == 1.0))

  stan_data <- list(
    n          = n,
    n_loc      = n_loc,
    p          = p,
    y          = as.integer(y_counts),
    units_m    = as.integer(units_m),
    is_mf      = as.integer(which_diag),
    ID_coords  = as.integer(ID_coords),
    D_mat      = as.matrix(dist(coords)),
    eta_fixed  = as.numeric(D %*% par$beta),
    mda_impact = as.numeric(mda_impact),
    use_mda    = use_mda_flag,
    omega      = par$k,
    alpha      = par$rho,
    gamma_sens = par$gamma_sens,
    sigma2     = par$sigma2,
    phi        = par$phi
  )

  mod     <- get_stan_model(model = "lf", backend = backend, messages = messages)
  backend <- mod$backend

  if (messages)
    message(sprintf("Sampling %d iter (%d warmup), %d chain(s) [lf_mdiag]...",
                    n_samples + n_warmup, n_warmup, n_chains))

  fit       <- .run_stan(mod$model, backend, stan_data,
                         n_samples, n_warmup, n_chains, n_cores,
                         adapt_delta, max_treedepth, messages)
  S_samples <- .extract_S(fit, backend, messages)

  result <- list(S_samples = S_samples, stan_fit = fit,
                 n_samples = nrow(S_samples), n_loc = n_loc,
                 coords = coords, par = par, backend = backend)
  class(result) <- "dsgm_spatial_samples"
  result
}


# =============================================================================
# Utility functions for spatial samples
# =============================================================================

##' @title Thin spatial process samples
##' @keywords internal
thin_spatial_samples <- function(spatial_samples, thin = 10) {
  keep <- seq(1, nrow(spatial_samples$S_samples), by = thin)
  spatial_samples$S_samples <- spatial_samples$S_samples[keep, , drop = FALSE]
  spatial_samples$n_samples  <- nrow(spatial_samples$S_samples)
  spatial_samples
}

##' @title Effective sample size for spatial process
##' @keywords internal
compute_spatial_ess <- function(spatial_samples) {
  if (!requireNamespace("coda", quietly = TRUE))
    stop("Package 'coda' is required for ESS calculation")
  sapply(seq_len(spatial_samples$n_loc), function(i)
    coda::effectiveSize(coda::as.mcmc(spatial_samples$S_samples[, i])))
}

##' @title Print method for dsgm_spatial_samples
##' @keywords internal
print.dsgm_spatial_samples <- function(x, ...) {
  cat("DSGM Spatial Process Samples\n")
  cat(sprintf("  Samples   : %d\n", x$n_samples))
  cat(sprintf("  Locations : %d\n", x$n_loc))
  cat(sprintf("  Matrix    : %d x %d\n", nrow(x$S_samples), ncol(x$S_samples)))
  if (requireNamespace("coda", quietly = TRUE)) {
    ess <- compute_spatial_ess(x)
    cat(sprintf("  ESS       : %.0f - %.0f  (median %.0f)\n",
                min(ess), max(ess), median(ess)))
  }
  invisible(x)
}


# =============================================================================
# Initial value functions
# =============================================================================

##' @title Compute initial parameter values for the STH model
##' @keywords internal
dsgm_initial_value <- function(y_prev, intensity_data, D, coords, ID_coords,
                               int_mat, survey_times_data, mda_times,
                               penalty, fix_alpha_W, fix_gamma_W,
                               start_pars, messages) {

  n <- length(y_prev)
  p <- ncol(D)

  .pen <- function(alpha_W, gamma_W, rho) {
    pen <- 0
    if (!is.null(penalty)) {
      if (is.null(fix_alpha_W)) {
        if (is.function(penalty$alpha))
          pen <- pen + penalty$alpha(alpha_W)
        else if (!is.null(penalty$alpha_param1))
          pen <- pen - (penalty$alpha_param1 - 1) * log(alpha_W) -
            (penalty$alpha_param2 - 1) * log(1 - alpha_W)
      }
      if (is.null(fix_gamma_W)) {
        if (is.function(penalty$gamma)) {
          pen <- pen + penalty$gamma(gamma_W)
        } else if (!is.null(penalty$gamma_type)) {
          if (penalty$gamma_type == "gamma")
            pen <- pen - (penalty$gamma_shape - 1) * log(gamma_W) + penalty$gamma_rate * gamma_W
          else if (penalty$gamma_type == "lognormal") {
            d <- log(gamma_W) - penalty$gamma_mean
            pen <- pen + 0.5 * d^2 / penalty$gamma_sd^2
          } else {
            d <- gamma_W - penalty$gamma_mean
            pen <- pen + 0.5 * d^2 / penalty$gamma_sd^2
          }
        } else if (!is.null(penalty$gamma_shape)) {
          pen <- pen - (penalty$gamma_shape - 1) * log(gamma_W) + penalty$gamma_rate * gamma_W
        } else if (!is.null(penalty$gamma_mean)) {
          d <- gamma_W - penalty$gamma_mean
          pen <- pen + 0.5 * d^2 / penalty$gamma_sd^2
        }
      }
    }
    if (!is.null(penalty$rho_type)) {
      if (penalty$rho_type == "gamma") {
        pen <- pen - (penalty$rho_mean - 1) * log(rho) + penalty$rho_sd * rho
      } else if (penalty$rho_type == "normal") {
        d <- rho - penalty$rho_mean
        pen <- pen + 0.5 * d^2 / penalty$rho_sd^2
      } else if (penalty$rho_type == "lognormal") {
        d <- log(rho) - penalty$rho_mean
        pen <- pen + 0.5 * d^2 / penalty$rho_sd^2
      }
    }
    pen
  }

  nll <- function(par) {
    beta    <- par[1:p]
    k       <- exp(par[p + 1])
    rho     <- exp(par[p + 2])
    idx     <- p + 3
    alpha_W <- if (is.null(fix_alpha_W)) { v <- plogis(par[idx]); idx <<- idx + 1; v } else fix_alpha_W
    gamma_W <- if (is.null(fix_gamma_W))   exp(par[idx])                               else fix_gamma_W

    mu_W <- exp(D %*% beta) *
      compute_mda_effect(survey_times_data, mda_times, int_mat,
                         alpha_W, gamma_W, kappa = 1)

    pr <- pmax(pmin(1 - (k / (k + mu_W * (1 - exp(-rho))))^k, 1 - 1e-10), 1e-10)
    ll <- sum(log(1 - pr[y_prev == 0])) + sum(log(pr[y_prev == 1]))

    pos   <- which(y_prev == 1)
    mupos <- mu_W[pos]; prpos <- pr[pos]
    mu_C  <- rho * mupos / prpos
    s2_C  <- pmax((rho * mupos * (1 + rho)) / prpos +
                    (rho^2 * mupos^2 / prpos) * (1 / k + 1 - 1 / prpos), 1e-10)
    kC    <- pmax((mu_C - 1)^2 / s2_C, 1e-10)
    tC    <- pmax(s2_C / (mu_C - 1), 1e-10)
    li    <- dgamma(intensity_data - 1, shape = kC, scale = tC, log = TRUE)
    li[!is.finite(li)] <- -1e10
    ll <- ll + sum(li)

    nll <- -(ll - .pen(alpha_W, gamma_W, rho))
    if (!is.finite(nll) || nll > 1e10) 1e10 else nll
  }

  mc1 <- mean(intensity_data - 1, na.rm = TRUE)
  vc1 <- var(intensity_data - 1,  na.rm = TRUE)
  k0   <- if (!is.null(start_pars$k)) start_pars$k else
    if (!is.na(vc1) && vc1 > 0 && mc1 > 0) max(mc1^2 / vc1, 0.1) else 0.5
  rho0 <- if (!is.null(start_pars$rho))     start_pars$rho     else 1
  aW0  <- if (!is.null(start_pars$alpha_W)) start_pars$alpha_W else 0.5
  gW0  <- if (!is.null(start_pars$gamma_W)) start_pars$gamma_W else 2.0
  b0   <- if (!is.null(start_pars$beta)) start_pars$beta else {
    b <- rep(0, p)
    b[1] <- log(max(mean(y_prev) * mean(intensity_data, na.rm = TRUE), 1)); b }

  par0 <- c(b0, log(k0), log(rho0))
  if (is.null(fix_alpha_W)) par0 <- c(par0, qlogis(aW0))
  if (is.null(fix_gamma_W)) par0 <- c(par0, log(gW0))

  if (messages) message("Optimising simplified model for initial values (STH)...")
  fit <- nlminb(par0, nll, control = list(eval.max = 2000, iter.max = 1000,
                                          trace = ifelse(messages, 1, 0)))
  pe <- if (fit$convergence != 0) {
    warning("STH initial value optimisation did not converge (code=",
            fit$convergence, "). Using starting values.")
    par0
  } else {
    if (messages) message(sprintf("  Converged: nll = %.2f", fit$objective))
    fit$par
  }

  beta_e <- pe[1:p]; k_e <- exp(pe[p + 1]); rho_e <- exp(pe[p + 2])
  idx <- p + 3
  aW_e <- if (is.null(fix_alpha_W)) { v <- plogis(pe[idx]); idx <- idx + 1; v } else fix_alpha_W
  gW_e <- if (is.null(fix_gamma_W))   exp(pe[idx])                               else fix_gamma_W

  s2_0  <- if (!is.null(start_pars$sigma2)) start_pars$sigma2 else 1.0
  phi_0 <- if (!is.null(start_pars$phi)) start_pars$phi else {
    d <- as.matrix(dist(coords)); median(d[upper.tri(d)]) / 3 }

  if (messages)
    message(sprintf("  k=%.3f  rho=%.3f  alpha_W=%.3f  gamma_W=%.3f  sigma2=%.3f  phi=%.3f",
                    k_e, rho_e, aW_e, gW_e, s2_0, phi_0))

  list(beta = beta_e, k = k_e, rho = rho_e, alpha_W = aW_e, gamma_W = gW_e,
       sigma2 = s2_0, phi = phi_0)
}


##' @title Compute initial parameter values for the LF multi-diagnostic model
##' @keywords internal
dsgm_initial_value_lf <- function(y_counts, units_m, which_diag, D, coords,
                                  int_mat, survey_times_data, mda_times,
                                  gamma_sens, penalty,
                                  fix_k, fix_alpha_W, fix_gamma_W,
                                  use_mda, start_pars, messages) {

  n <- length(y_counts)
  p <- ncol(D)

  nll <- function(par) {
    beta <- par[1:p]
    k    <- if (!is.null(fix_k)) fix_k else exp(par[p + 1])
    rho  <- exp(par[p + 2])
    idx  <- p + 3

    alpha_W <- if (use_mda && is.null(fix_alpha_W)) { v <- plogis(par[idx]); idx <<- idx + 1; v } else fix_alpha_W
    gamma_W <- if (use_mda && is.null(fix_gamma_W))   exp(par[idx])                               else fix_gamma_W

    mu <- exp(D %*% beta)
    if (use_mda)
      mu <- mu * compute_mda_effect(survey_times_data, mda_times, int_mat,
                                    alpha_W, gamma_W, kappa = 1)

    prob <- numeric(n)
    mf   <- which_diag == 1
    cfa  <- which_diag == 0
    prob[mf]  <- pmax(pmin(1 - (k / (k + mu[mf]  * (1 - exp(-rho))))^k, 1 - 1e-10), 1e-10)
    prob[cfa] <- pmax(pmin(gamma_sens * (1 - (k / (k + mu[cfa]))^k),     1 - 1e-10), 1e-10)

    ll <- sum(y_counts * log(prob) + (units_m - y_counts) * log(1 - prob))
    if (!is.finite(ll)) 1e10 else -ll
  }

  obs_prev <- mean(y_counts / pmax(units_m, 1))
  b0  <- if (!is.null(start_pars$beta)) start_pars$beta else {
    b <- rep(0, p); b[1] <- log(max(obs_prev, 0.01)); b }
  k0   <- if (!is.null(fix_k)) fix_k else
    if (!is.null(start_pars$k)) start_pars$k else 0.5
  rho0 <- if (!is.null(start_pars$rho))     start_pars$rho     else 0.5
  aW0  <- if (!is.null(fix_alpha_W)) fix_alpha_W else
    if (!is.null(start_pars$alpha_W)) start_pars$alpha_W else 0.5
  gW0  <- if (!is.null(fix_gamma_W)) fix_gamma_W else
    if (!is.null(start_pars$gamma_W)) start_pars$gamma_W else 2.0

  par0 <- c(b0, log(k0), log(rho0))
  if (use_mda) {
    if (is.null(fix_alpha_W)) par0 <- c(par0, qlogis(aW0))
    if (is.null(fix_gamma_W)) par0 <- c(par0, log(gW0))
  }

  if (messages) message("Optimising simplified model for initial values (LF)...")
  fit <- nlminb(par0, nll, control = list(eval.max = 2000, iter.max = 1000,
                                          trace = ifelse(messages, 1, 0)))
  pe <- if (fit$convergence != 0) {
    warning("LF initial value optimisation did not converge (code=",
            fit$convergence, "). Using starting values.")
    par0
  } else {
    if (messages) message(sprintf("  Converged: nll = %.2f", fit$objective))
    fit$par
  }

  beta_e <- pe[1:p]
  k_e    <- if (!is.null(fix_k)) fix_k else exp(pe[p + 1])
  rho_e  <- exp(pe[p + 2])
  idx    <- p + 3
  aW_e <- if (use_mda && is.null(fix_alpha_W)) { v <- plogis(pe[idx]); idx <- idx + 1; v } else aW0
  gW_e <- if (use_mda && is.null(fix_gamma_W))   exp(pe[idx])                               else gW0

  s2_0   <- if (!is.null(start_pars$sigma2)) start_pars$sigma2 else 1.0
  tau2_0 <- if (!is.null(start_pars$tau2))   start_pars$tau2   else 0.1
  phi_0  <- if (!is.null(start_pars$phi)) start_pars$phi else {
    d <- as.matrix(dist(coords)); median(d[upper.tri(d)]) / 3 }

  if (messages)
    message(sprintf("  k=%.3f  rho=%.3f  sigma2=%.3f  phi=%.3f  tau2=%.3f",
                    k_e, rho_e, s2_0, phi_0, tau2_0))

  list(beta = beta_e, k = k_e, rho = rho_e, gamma_sens = gamma_sens,
       sigma2 = s2_0, phi = phi_0, tau2 = tau2_0, alpha_W = aW_e, gamma_W = gW_e)
}


# =============================================================================
# TMB fitting engines
# =============================================================================

##' @title Fit DSGM using TMB
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
                         model             = "sth",
                         y_counts          = NULL,
                         units_m           = NULL,
                         which_diag        = NULL,
                         gamma_sens        = 0.97,
                         fix_k             = NULL,
                         fix_tau2          = NULL,
                         use_mda           = TRUE,
                         intensity_family  = 0L,
                         use_hessian_refinement = TRUE,
                         messages          = TRUE) {

  model <- match.arg(model, c("sth", "lf_mdiag"))

  if (model == "lf_mdiag") {
    return(.dsgm_fit_tmb_lf_mdiag(
      y_counts = y_counts, units_m = units_m, which_diag = which_diag,
      D = D, coords = coords, ID_coords = ID_coords, cov_offset = cov_offset,
      gamma_sens = gamma_sens, fix_k = fix_k, fix_tau2 = fix_tau2,
      use_mda = use_mda, int_mat = int_mat,
      survey_times_data = survey_times_data, mda_times = mda_times,
      fix_alpha_W = fix_alpha_W, fix_gamma_W = fix_gamma_W,
      penalty = penalty, par0 = par0, S_samples_obj = S_samples_obj,
      messages = messages
    ))
  }

  n         <- length(y_prev)
  n_loc     <- nrow(coords)
  S_samples <- S_samples_obj$S_samples
  pos_idx   <- which(y_prev == 1) - 1L

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

  obj_d <- TMB::MakeADFun(
    data       = make_data(1L, numeric(nrow(S_samples))),
    parameters = parameters,
    map        = tmb_map,
    DLL        = "RiskMap",
    silent     = TRUE
  )
  obj_d$fn()
  log_denominator_vals <- obj_d$report()$log_f_vals

  if (messages) message("Building TMB objective for optimization...")
  obj <- TMB::MakeADFun(
    data       = make_data(0L, log_denominator_vals),
    parameters = parameters,
    map        = tmb_map,
    DLL        = "RiskMap",
    silent     = !messages
  )

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
      list(beta   = obj$env$parList()$beta,
           k      = par_est["k",      "Estimate"],
           rho    = par_est["rho",    "Estimate"],
           sigma2 = par_est["sigma2", "Estimate"],
           phi    = par_est["phi",    "Estimate"]),
      if (use_mda) list(alpha_W = par_est["alpha_W", "Estimate"],
                        gamma_W = par_est["gamma_W", "Estimate"])
    ),
    params_se = c(
      list(beta   = summary(sdr, "fixed")[grepl("beta", rownames(summary(sdr, "fixed"))), "Std. Error"],
           k      = par_est["k",      "Std. Error"],
           rho    = par_est["rho",    "Std. Error"],
           sigma2 = par_est["sigma2", "Std. Error"],
           phi    = par_est["phi",    "Std. Error"]),
      if (use_mda) list(alpha_W = par_est["alpha_W", "Std. Error"],
                        gamma_W = par_est["gamma_W", "Std. Error"])
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
    log_nu2       = log(max(if (!is.null(fix_tau2)) fix_tau2 else par0$tau2, 1e-10) / par0$sigma2),
    log_omega     = log(if (!is.null(fix_k)) fix_k else par0$k),
    log_alpha     = log(par0$rho),
    logit_alpha_W = if (use_mda) qlogis(par0$alpha_W) else 0.0,
    log_gamma_W   = if (use_mda) log(par0$gamma_W)    else 0.0
  )

  map_list <- list()
  if (!is.null(fix_tau2)) map_list$log_nu2    <- factor(NA)
  if (!is.null(fix_k))    map_list$log_omega  <- factor(NA)
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


# =============================================================================
# Main user-facing function
# =============================================================================

##' @title Fit a Doubly Stochastic Geostatistical Model (DSGM)
##' @param intensity_family Distribution for C | C > 0 (STH model only).
##'   \code{"gamma"} (default) = shifted Gamma; \code{"negbin"} = zero-truncated NegBin.
##' @export
dsgm <- function(formula,
                 data,
                 model            = c("sth", "lf_mdiag"),
                 intensity_family = c("gamma", "negbin"),
                 time             = NULL,
                 mda_times        = NULL,
                 int_mat          = NULL,
                 den              = NULL,
                 gamma_sens       = 0.97,
                 fix_k            = NULL,
                 use_mda          = NULL,
                 penalty          = NULL,
                 drop_W           = NULL,
                 decay_W          = NULL,
                 crs              = NULL,
                 convert_to_crs   = NULL,
                 scale_to_km      = TRUE,
                 par0             = NULL,
                 n_samples        = 1000,
                 n_warmup         = 1000,
                 n_chains         = 1,
                 adapt_delta      = 0.8,
                 max_treedepth    = 10,
                 return_samples   = TRUE,
                 backend          = NULL,
                 messages         = TRUE,
                 start_pars       = list(beta    = NULL,
                                         k       = NULL,
                                         rho     = NULL,
                                         sigma2  = NULL,
                                         phi     = NULL,
                                         tau2    = NULL,
                                         alpha_W = NULL,
                                         gamma_W = NULL)) {

  model            <- match.arg(model)
  intensity_family <- match.arg(intensity_family)
  intensity_family_int <- if (intensity_family == "gamma") 0L else 1L

  if (!inherits(formula, "formula"))
    stop("'formula' must be a formula object")
  if (!inherits(data, c("data.frame", "sf")))
    stop("'data' must be a data.frame or sf object")
  if (n_samples <= 0 || n_warmup < 0)
    stop("'n_samples' must be positive and 'n_warmup' non-negative")
  if (n_chains < 1) stop("'n_chains' must be at least 1")
  if (adapt_delta <= 0 || adapt_delta >= 1)
    stop("'adapt_delta' must be in (0, 1)")

  if (is.null(use_mda)) use_mda <- !is.null(mda_times)

  fix_alpha_W <- drop_W
  fix_gamma_W <- decay_W
  no_penalty  <- is.null(penalty)
  n           <- nrow(data)

  inter_f  <- interpret.formula(formula)
  gp_terms <- inter_f$gp.spec$term
  gp_dim   <- inter_f$gp.spec$dim

  if (!(length(gp_terms) == 1 && gp_terms[1] == "sf") && gp_dim != 2)
    stop("Specify a 2-D spatial GP: gp(x, y) or gp(sf)")

  gp_nugget <- inter_f$gp.spec$nugget
  fix_tau2  <- if (is.null(gp_nugget)) NULL else as.numeric(gp_nugget)

  mf         <- model.frame(inter_f$pf, data = data, na.action = na.fail)
  D          <- as.matrix(model.matrix(attr(mf, "terms"), data = data))
  p          <- ncol(D)
  cov_offset <- if (is.null(inter_f$offset)) rep(0, n) else data[[inter_f$offset]]

  den_name <- deparse(substitute(den))
  if (den_name == "NULL") {
    den_vals <- rep(1L, n)
  } else {
    if (!den_name %in% names(data))
      stop(sprintf("'den' column '%s' not found in 'data'", den_name))
    den_vals <- as.integer(data[[den_name]])
    if (any(den_vals < 1, na.rm = TRUE)) stop("'den' values must all be >= 1")
    if (any(is.na(den_vals)))            stop("Missing values detected in 'den' column")
  }

  if (inherits(data, "sf")) {
    if (!is.null(convert_to_crs)) {
      data <- sf::st_transform(data, convert_to_crs); crs <- convert_to_crs
    } else {
      crs <- sf::st_crs(data)$input
    }
    coords_all <- sf::st_coordinates(data)
  } else {
    cn <- gp_terms[1:2]
    if (!all(cn %in% names(data)))
      stop("Coordinate columns specified in gp() not found in 'data'")
    coords_all <- as.matrix(data[, cn])
  }
  if (scale_to_km) {
    coords_all <- coords_all / 1000
    if (messages) message("Coordinates scaled to kilometres")
  }
  coords_u  <- unique(coords_all)
  n_loc     <- nrow(coords_u)
  ID_coords <- apply(coords_all, 1, function(x)
    which(apply(coords_u, 1, function(y) all(abs(x - y) < 1e-10)))[1])
  if (messages) message(sprintf("Identified %d unique spatial locations", n_loc))

  if (use_mda) {
    time_name <- deparse(substitute(time))
    if (time_name == "NULL" || !time_name %in% names(data))
      stop("'time' column not found in 'data' (required when use_mda = TRUE)")
    survey_times_data <- data[[time_name]]
    if (any(is.na(survey_times_data))) stop("Missing values in 'time' variable")
    if (is.null(mda_times)) stop("'mda_times' required when use_mda = TRUE")
    if (is.null(int_mat))   stop("'int_mat' required when use_mda = TRUE")
    int_mat <- as.matrix(int_mat)
    if (nrow(int_mat) != n)
      stop(sprintf("'int_mat' must have %d rows", n))
    if (ncol(int_mat) != length(mda_times))
      stop(sprintf("'int_mat' must have %d columns", length(mda_times)))
  } else {
    survey_times_data <- rep(0, n)
    mda_times         <- 0
    int_mat           <- matrix(0, nrow = n, ncol = 1)
  }

  # ===========================================================================
  # STH branch
  # ===========================================================================
  if (model == "sth") {

    egg_counts <- as.numeric(model.response(mf))
    if (any(egg_counts < 0, na.rm = TRUE)) stop("Egg counts cannot be negative")
    if (any(is.na(egg_counts)))            stop("Missing values in egg count response")

    y_prev         <- as.integer(egg_counts > 0)
    intensity_data <- egg_counts[egg_counts > 0]
    n_pos          <- sum(y_prev)
    if (n_pos == 0) stop("No positive egg counts; model cannot be fitted")

    if (messages)
      message(sprintf("STH data: %d observations, %d (%.1f%%) egg-positive [intensity: %s]",
                      n, n_pos, 100 * n_pos / n, intensity_family))

    if (is.null(par0)) {
      if (messages) message("\n=== Computing initial parameter values (STH) ===")
      par0 <- dsgm_initial_value(
        y_prev = y_prev, intensity_data = intensity_data, D = D,
        coords = coords_u, ID_coords = ID_coords,
        int_mat = int_mat, survey_times_data = survey_times_data,
        mda_times = mda_times, penalty = penalty,
        fix_alpha_W = fix_alpha_W, fix_gamma_W = fix_gamma_W,
        start_pars = start_pars, messages = messages)
    }
    req  <- c("beta", "k", "rho", "sigma2", "phi")
    if (is.null(fix_alpha_W)) req <- c(req, "alpha_W")
    if (is.null(fix_gamma_W)) req <- c(req, "gamma_W")
    miss <- setdiff(req, names(par0))
    if (length(miss) > 0)
      stop("Missing initial parameters: ", paste(miss, collapse = ", "))

    if (messages) {
      message("\n=== Sampling spatial process (STH, Stan) ===")
      message(sprintf("  n_samples=%d  n_warmup=%d  n_chains=%d  adapt_delta=%.2f",
                      n_samples, n_warmup, n_chains, adapt_delta))
    }
    sp <- sample_spatial_process_stan(
      y_prev            = y_prev,
      intensity_data    = intensity_data,
      D                 = D,
      coords            = coords_u,
      ID_coords         = ID_coords,
      int_mat           = int_mat,
      survey_times_data = survey_times_data,
      mda_times         = mda_times,
      par               = par0,
      n_samples         = n_samples,
      n_warmup          = n_warmup,
      n_chains          = n_chains,
      n_cores           = 1,
      adapt_delta       = adapt_delta,
      max_treedepth     = max_treedepth,
      intensity_family  = intensity_family_int,
      backend           = backend,
      messages          = messages)

    if (messages) message("\n=== Fitting via TMB-MCML (STH) ===")
    fit <- dsgm_fit_tmb(
      model             = "sth",
      y_prev            = y_prev,
      intensity_data    = intensity_data,
      D                 = D,
      coords            = coords_u,
      ID_coords         = ID_coords,
      int_mat           = int_mat,
      survey_times_data = survey_times_data,
      mda_times         = mda_times,
      par0              = par0,
      cov_offset        = cov_offset,
      fix_alpha_W       = fix_alpha_W,
      fix_gamma_W       = fix_gamma_W,
      fix_tau2          = fix_tau2,
      use_mda           = use_mda,
      penalty           = penalty,
      S_samples_obj     = sp,
      intensity_family  = intensity_family_int,
      use_hessian_refinement = TRUE,
      messages          = messages)

    res <- list(
      family            = "intprev",
      intensity_family  = intensity_family,
      prevalence_data   = y_prev,
      intensity_data    = intensity_data,
      egg_counts        = egg_counts,
      n_positive        = n_pos,
      D                 = D,
      coords            = coords_u,
      mda_times         = mda_times,
      survey_times_data = survey_times_data,
      int_mat           = int_mat,
      ID_coords         = ID_coords,
      fix_alpha_W       = fix_alpha_W,
      fix_gamma_W       = fix_gamma_W,
      use_mda           = use_mda,
      formula           = formula,
      crs               = crs,
      scale_to_km       = scale_to_km,
      data_sf           = data,
      kappa             = 0.5,
      cov_offset        = cov_offset,
      penalty           = if (!no_penalty) penalty else NULL,
      model_params      = fit$params,
      params_se         = fit$params_se,
      tmb_sdr           = fit$tmb_sdr,
      tmb_obj           = fit$tmb_obj,
      convergence       = fit$convergence,
      log_likelihood    = fit$log_likelihood,
      n_locations       = n_loc,
      n_observations    = n,
      n_covariates      = p,
      n_mda_rounds      = if (use_mda) length(mda_times) else 0L,
      stan_settings     = list(n_samples = n_samples, n_warmup = n_warmup,
                               n_chains = n_chains, adapt_delta = adapt_delta,
                               max_treedepth = max_treedepth),
      call              = match.call()
    )
    if (return_samples) res$spatial_samples <- sp
    class(res) <- "RiskMap"
    return(res)
  }

  # ===========================================================================
  # LF multi-diagnostic branch
  # ===========================================================================
  if (model == "lf_mdiag") {

    if ("diagnostic" %in% names(data)) {
      diag_raw <- data[["diagnostic"]]
      if (!all(diag_raw %in% c("par", "ser")))
        stop("Column 'diagnostic' must contain only 'par' or 'ser'.")
      which_diag <- ifelse(diag_raw == "par", 1L, 0L)
    } else {
      if (messages)
        message("No 'diagnostic' column found: assuming all observations are parasitological (MF).")
      which_diag <- rep(1L, n)
    }
    if (gamma_sens <= 0 || gamma_sens > 1)
      stop("'gamma_sens' must be in (0, 1]")

    y_counts <- as.integer(model.response(mf))
    if (any(y_counts < 0, na.rm = TRUE)) stop("Response (positive counts) cannot be negative")
    if (any(is.na(y_counts)))            stop("Missing values in response variable")

    units_m <- den_vals

    if (messages)
      message(sprintf("LF data: %d observations (%d MF, %d serological)",
                      n, sum(which_diag == 1), sum(which_diag == 0)))

    if (is.null(par0)) {
      if (messages) message("\n=== Computing initial parameter values (LF) ===")
      par0 <- dsgm_initial_value_lf(
        y_counts = y_counts, units_m = units_m, which_diag = which_diag, D = D,
        coords = coords_u,
        int_mat = int_mat, survey_times_data = survey_times_data,
        mda_times = mda_times, gamma_sens = gamma_sens, penalty = penalty,
        fix_k = fix_k, fix_alpha_W = fix_alpha_W, fix_gamma_W = fix_gamma_W,
        use_mda = use_mda, start_pars = start_pars, messages = messages)
    }
    par0$gamma_sens <- gamma_sens
    if (!is.null(fix_tau2)) {
      par0$tau2 <- fix_tau2
    } else if (is.null(par0$tau2)) {
      par0$tau2 <- if (!is.null(start_pars$tau2)) start_pars$tau2 else 0.1
    }

    mda_impact <- if (use_mda)
      compute_mda_effect(survey_times_data, mda_times, int_mat,
                         par0$alpha_W, par0$gamma_W, kappa = 1)
    else
      rep(1.0, n)

    if (messages) {
      message("\n=== Sampling spatial process (LF, Stan) ===")
      message(sprintf("  n_samples=%d  n_warmup=%d  n_chains=%d  adapt_delta=%.2f",
                      n_samples, n_warmup, n_chains, adapt_delta))
    }
    sp <- sample_spatial_process_stan_lf(
      y_counts = y_counts, units_m = units_m, which_diag = which_diag, D = D,
      coords = coords_u, ID_coords = ID_coords, par = par0,
      mda_impact = mda_impact,
      n_samples = n_samples, n_warmup = n_warmup,
      n_chains = n_chains, n_cores = 1,
      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
      backend = backend, messages = messages)

    if (messages) message("\n=== Fitting via TMB-MCML (LF) ===")
    fit <- dsgm_fit_tmb(
      model             = "lf_mdiag",
      y_counts          = y_counts,
      units_m           = units_m,
      which_diag        = which_diag,
      D                 = D,
      coords            = coords_u,
      ID_coords         = ID_coords,
      int_mat           = int_mat,
      survey_times_data = survey_times_data,
      mda_times         = mda_times,
      par0              = par0,
      cov_offset        = cov_offset,
      gamma_sens        = gamma_sens,
      fix_k             = fix_k,
      use_mda           = use_mda,
      fix_alpha_W       = fix_alpha_W,
      fix_gamma_W       = fix_gamma_W,
      fix_tau2          = fix_tau2,
      penalty           = penalty,
      S_samples_obj     = sp,
      use_hessian_refinement = TRUE,
      messages          = messages)

    res <- list(
      family            = "lf_mdiag",
      y_counts          = y_counts,
      units_m           = units_m,
      which_diag        = which_diag,
      D                 = D,
      coords            = coords_u,
      mda_times         = mda_times,
      survey_times_data = survey_times_data,
      int_mat           = int_mat,
      ID_coords         = ID_coords,
      fix_alpha_W       = fix_alpha_W,
      fix_gamma_W       = fix_gamma_W,
      fix_k             = fix_k,
      gamma_sens        = gamma_sens,
      use_mda           = use_mda,
      formula           = formula,
      crs               = crs,
      scale_to_km       = scale_to_km,
      data_sf           = data,
      kappa             = 0.5,
      cov_offset        = cov_offset,
      penalty           = if (!no_penalty) penalty else NULL,
      model_params      = fit$params,
      params_se         = fit$params_se,
      tmb_sdr           = fit$tmb_sdr,
      tmb_obj           = fit$tmb_obj,
      convergence       = fit$convergence,
      log_likelihood    = fit$log_likelihood,
      n_locations       = n_loc,
      n_observations    = n,
      n_covariates      = p,
      n_mda_rounds      = length(mda_times),
      stan_settings     = list(n_samples = n_samples, n_warmup = n_warmup,
                               n_chains = n_chains, adapt_delta = adapt_delta,
                               max_treedepth = max_treedepth),
      call              = match.call()
    )
    if (return_samples) res$spatial_samples <- sp
    class(res) <- "RiskMap"
    return(res)
  }
}
