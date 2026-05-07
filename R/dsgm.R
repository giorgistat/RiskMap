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
# dsgm_initial_value  —  now accepts vary_k and omega1_start
# =============================================================================

##' @title Compute initial parameter values for the STH model
##' @param vary_k Logical; if TRUE, include omega1 (slope of log(omega) on
##'   log(mu_W)) in the simplified NLL optimisation.
##' @param omega1_start Starting value for omega1 when vary_k = TRUE.
##' @keywords internal
dsgm_initial_value <- function(y_prev, intensity_data, D, coords, ID_coords,
                               int_mat, survey_times_data, mda_times,
                               penalty, fix_alpha_W, fix_gamma_W,
                               vary_k = FALSE, omega1_start = 0.5,
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
    idx     <- idx + as.integer(is.null(fix_gamma_W))
    # omega1 slope (only estimated when vary_k = TRUE)
    omega1  <- if (vary_k) par[idx] else 0.0

    mu_star <- exp(D %*% beta) *
      compute_mda_effect(survey_times_data, mda_times, int_mat,
                         alpha_W, gamma_W, kappa = 1)

    # Per-observation k when vary_k = TRUE
    k_vec <- if (vary_k) {
      pmax(exp(log(k) + omega1 * log(mu_star + 1e-10)), 1e-6)
    } else {
      rep(k, n)
    }

    c_rho <- 1 - exp(-rho)
    pr <- pmax(pmin(1 - (k_vec / (k_vec + mu_star * c_rho))^k_vec, 1 - 1e-10), 1e-10)

    ll <- sum(log(1 - pr[y_prev == 0])) + sum(log(pr[y_prev == 1]))

    pos    <- which(y_prev == 1)
    mu_pos <- mu_star[pos]; k_pos <- k_vec[pos]; pr_pos <- pr[pos]

    mu_C <- rho * mu_pos / pr_pos
    s2_C <- pmax((rho * mu_pos * (1 + rho)) / pr_pos +
                   (rho^2 * mu_pos^2 / pr_pos) * (1 / k_pos + 1 - 1 / pr_pos), 1e-10)
    kC   <- pmax((mu_C - 1)^2 / s2_C, 1e-10)
    tC   <- pmax(s2_C / (mu_C - 1), 1e-10)
    li   <- dgamma(intensity_data - 1, shape = kC, scale = tC, log = TRUE)
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
  o1_0 <- if (!is.null(start_pars$omega1))  start_pars$omega1  else omega1_start
  b0   <- if (!is.null(start_pars$beta)) start_pars$beta else {
    b <- rep(0, p)
    b[1] <- log(max(mean(y_prev) * mean(intensity_data, na.rm = TRUE), 1)); b }

  par0 <- c(b0, log(k0), log(rho0))
  if (is.null(fix_alpha_W)) par0 <- c(par0, qlogis(aW0))
  if (is.null(fix_gamma_W)) par0 <- c(par0, log(gW0))
  if (vary_k)               par0 <- c(par0, o1_0)

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
  gW_e <- if (is.null(fix_gamma_W))   { v <- exp(pe[idx]);   idx <- idx + 1; v } else fix_gamma_W
  o1_e <- if (vary_k) pe[idx] else 0.0

  s2_0  <- if (!is.null(start_pars$sigma2)) start_pars$sigma2 else 1.0
  phi_0 <- if (!is.null(start_pars$phi)) start_pars$phi else {
    d <- as.matrix(dist(coords)); median(d[upper.tri(d)]) / 3 }

  if (messages) {
    msg <- sprintf("  k=%.3f  rho=%.3f  alpha_W=%.3f  gamma_W=%.3f  sigma2=%.3f  phi=%.3f",
                   k_e, rho_e, aW_e, gW_e, s2_0, phi_0)
    if (vary_k) msg <- paste0(msg, sprintf("  omega1=%.3f", o1_e))
    message(msg)
  }

  out <- list(beta = beta_e, k = k_e, rho = rho_e, alpha_W = aW_e, gamma_W = gW_e,
              sigma2 = s2_0, phi = phi_0)
  if (vary_k) out$omega1 <- o1_e
  out
}

##' @title Fit a Doubly Stochastic Geostatistical Model (DSGM)
##'
##' @param vary_k Logical (STH model only). If \code{FALSE} (default), the
##'   NegBin aggregation parameter \eqn{\omega} is constant across all
##'   observations. If \code{TRUE}, it varies with mean worm burden according to
##'   \eqn{\log(\omega_i) = \omega_0 + \omega_1 \log(\mu_{W,i})}, adding one
##'   parameter \eqn{\omega_1}. Motivated by the empirical Coffeng relationship
##'   whereby more-burdened populations show less aggregation.
##' @param omega1_start Starting value for \eqn{\omega_1} when \code{vary_k = TRUE}.
##'   Default is 0.5, consistent with estimates from Kenya hookworm survey data.
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
                 vary_k           = FALSE,       # <-- NEW
                 omega1_start     = 0.5,         # <-- NEW
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
                                         gamma_W = NULL,
                                         omega1  = NULL)) {  # <-- omega1 added

  model            <- match.arg(model)
  intensity_family <- match.arg(intensity_family)
  intensity_family_int <- if (intensity_family == "gamma") 0L else 1L

  # vary_k only makes sense for the STH model
  if (vary_k && model != "sth") {
    warning("vary_k = TRUE is only supported for model = 'sth'. Ignoring.")
    vary_k <- FALSE
  }

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

    if (messages) {
      msg <- sprintf("STH data: %d observations, %d (%.1f%%) egg-positive [intensity: %s]",
                     n, n_pos, 100 * n_pos / n, intensity_family)
      if (vary_k) msg <- paste0(msg, "  [omega ~ mu_W]")
      message(msg)
    }

    if (is.null(par0)) {
      if (messages) {
        message("\n=== Computing initial parameter values (STH) ===")
      }

      par0 <- dsgm_initial_value(
        y_prev            = y_prev,
        intensity_data    = intensity_data,
        D                 = D,
        coords            = coords_u,
        ID_coords         = ID_coords,
        int_mat           = int_mat,
        survey_times_data = survey_times_data,
        mda_times         = mda_times,
        penalty           = penalty,
        fix_alpha_W       = fix_alpha_W,
        fix_gamma_W       = fix_gamma_W,
        vary_k            = vary_k,           # <-- pass through
        omega1_start      = omega1_start,     # <-- pass through
        start_pars        = start_pars,
        messages          = messages)
    }
    req  <- c("beta", "k", "rho", "sigma2", "phi")
    if (is.null(fix_alpha_W)) req <- c(req, "alpha_W")
    if (is.null(fix_gamma_W)) req <- c(req, "gamma_W")
    miss <- setdiff(req, names(par0))
    if (length(miss) > 0)
      stop("Missing initial parameters: ", paste(miss, collapse = ", "))

    # For the Stan sampling step, pass an effective scalar k to the Stan model.
    # When vary_k = TRUE, k in par0 is omega_0 (k at mu_W = 1). This is an
    # approximation for the proposal — the TMB step handles the exact likelihood.
    # No changes to sample_spatial_process_stan() are needed.

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
      vary_k            = vary_k,
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
    if(vary_k) {
      message("\n * Varying aggregation parameter as a function of worm burden")
    } else {
      message("\n * Constant aggregation parameter")
    }
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
      vary_k            = vary_k,             # <-- pass through
      omega1_start      = if (!is.null(par0$omega1)) par0$omega1 else omega1_start,
      use_hessian_refinement = TRUE,
      messages          = messages)

    res <- list(
      family            = "intprev",
      intensity_family  = intensity_family,
      vary_k            = vary_k,
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
  # LF multi-diagnostic branch  (UNCHANGED from original)
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
      vary_k            = FALSE,
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


##' @title Sample spatial process for the STH model
##' @param vary_k Logical; if TRUE, passes log_k and omega1 to Stan so that
##'   the per-observation aggregation parameter is computed as
##'   exp(log_k + omega1 * log(mu_W[i])). Must match the vary_k setting used
##'   in dsgm_fit_tmb().
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
                                        vary_k           = FALSE,
                                        backend          = NULL,
                                        messages         = TRUE) {

  n       <- length(y_prev)
  n_loc   <- nrow(coords)
  p       <- ncol(D)
  pos_idx <- which(y_prev == 1)
  n_pos   <- length(pos_idx)

  mda_impact <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                                   par$alpha_W, par$gamma_W, kappa = 1)

  # When vary_k = TRUE, pass log_k and omega1 separately so Stan can compute
  # k_i = exp(log_k + omega1 * log(mu_W[i])).
  # When vary_k = FALSE, omega1 = 0 so k_i = exp(log_k) = k (constant).
  log_k_val  <- log(par$k)
  omega1_val <- if (vary_k && !is.null(par$omega1)) par$omega1 else 0.0

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
    log_k            = log_k_val,    # replaces scalar k
    omega1           = omega1_val,   # = 0 when vary_k = FALSE
    rho              = par$rho,
    sigma2           = par$sigma2,
    phi              = par$phi,
    vary_k           = as.integer(vary_k),
    intensity_family = as.integer(intensity_family)
  )

  mod     <- get_stan_model(model = "sth", backend = backend, messages = messages)
  backend <- mod$backend

  if (messages)
    message(sprintf(
      "Sampling %d iter (%d warmup), %d chain(s) [sth, family=%s, vary_k=%s]...",
      n_samples + n_warmup, n_warmup, n_chains,
      ifelse(intensity_family == 0L, "shifted Gamma", "zero-trunc NegBin"),
      ifelse(vary_k, "TRUE", "FALSE")))

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

