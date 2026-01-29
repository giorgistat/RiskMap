##' @title Sample spatial process using STAN
##' @description Uses STAN to sample S(x) from posterior given fixed parameters
##' @param y_prev Binary infection indicators
##' @param intensity_data Egg counts for positives only
##' @param D Design matrix
##' @param coords Unique spatial coordinates (n_loc x 2)
##' @param ID_coords Location indices for each observation
##' @param int_mat Intervention matrix
##' @param survey_times_data Survey times
##' @param mda_times MDA times
##' @param par Parameters list with beta, k, rho, alpha_W, gamma_W, sigma2, phi
##' @param n_samples Number of MCMC samples to draw
##' @param n_warmup Number of warmup iterations
##' @param n_chains Number of MCMC chains
##' @param n_cores Number of cores for parallel chains
##' @param adapt_delta Target acceptance rate (default 0.8)
##' @param max_treedepth Maximum tree depth (default 10)
##' @param messages Print STAN messages
##' @return Matrix of S samples (n_samples x n_loc)
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
                                        n_samples = 1000,
                                        n_warmup = 1000,
                                        n_chains = 4,
                                        n_cores = 4,
                                        adapt_delta = 0.8,
                                        max_treedepth = 10,
                                        messages = TRUE) {

  if(!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required for this function")
  }

  # =============================================================================
  # PREPARE DATA FOR STAN
  # =============================================================================

  n <- length(y_prev)
  n_loc <- nrow(coords)
  p <- ncol(D)

  # Compute distance matrix
  dist_mat <- as.matrix(dist(coords))

  # Compute MDA impact
  mda_impact <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                                   par$alpha_W, par$gamma_W, kappa = 1)

  # Linear predictor (without spatial effect)
  eta_fixed <- as.numeric(D %*% par$beta)

  # Extract positives
  pos_idx <- which(y_prev == 1)
  n_pos <- length(pos_idx)

  # Create STAN data list
  stan_data <- list(
    n = n,
    n_loc = n_loc,
    n_pos = n_pos,
    p = p,

    # Response data
    y = y_prev,
    C_pos = intensity_data,
    pos_idx = pos_idx,

    # Location mapping
    ID_coords = ID_coords,

    # Distance matrix
    D_mat = dist_mat,

    # Fixed effects
    eta_fixed = eta_fixed,
    mda_impact = mda_impact,

    # Fixed parameters
    k = par$k,
    rho = par$rho,
    sigma2 = par$sigma2,
    phi = par$phi
  )

  # =============================================================================
  # STAN MODEL CODE
  # =============================================================================

  stan_code <- "
  functions {
    // Compute correlation matrix using exponential kernel
    matrix compute_corr_matrix(matrix D, real phi) {
      int n = rows(D);
      matrix[n, n] R;

      for (i in 1:n) {
        R[i, i] = 1.0;
        for (j in (i+1):n) {
          R[i, j] = exp(-D[i, j] / phi);
          R[j, i] = R[i, j];
        }
      }

      return R;
    }

    // Compute prevalence from mean worm burden
    real compute_prevalence(real mu_W, real k, real rho) {
      real term = k / (k + mu_W * (1 - exp(-rho)));
      real pr = 1 - pow(term, k);
      return pr;
    }

    // Compute conditional mean intensity
    real compute_mu_C(real mu_W, real pr, real rho) {
      return (rho * mu_W) / pr;
    }

    // Compute conditional variance
    real compute_sigma2_C(real mu_W, real pr, real k, real rho) {
      real var1 = (rho * mu_W * (1 + rho)) / pr;
      real var2 = (square(rho) * square(mu_W) / pr) * (1.0/k + 1 - 1.0/pr);
      return var1 + var2;
    }

    // Compute shifted Gamma parameters
    vector compute_gamma_params(real mu_W, real pr, real k, real rho) {
      vector[2] params;
      real mu_C = compute_mu_C(mu_W, pr, rho);
      real sigma2_C = compute_sigma2_C(mu_W, pr, k, rho);

      // kappa = (mu_C - 1)^2 / sigma2_C
      // theta = sigma2_C / (mu_C - 1)
      params[1] = square(mu_C - 1) / sigma2_C;  // kappa
      params[2] = sigma2_C / (mu_C - 1);        // theta

      return params;
    }
  }

  data {
    int<lower=1> n;              // Number of observations
    int<lower=1> n_loc;          // Number of unique locations
    int<lower=0> n_pos;          // Number of positive observations
    int<lower=1> p;              // Number of covariates

    // Response data
    int<lower=0, upper=1> y[n];  // Binary infection indicator
    vector<lower=1>[n_pos] C_pos; // Egg counts for positives
    int<lower=1, upper=n> pos_idx[n_pos]; // Indices of positives

    // Location mapping
    int<lower=1, upper=n_loc> ID_coords[n];

    // Distance matrix
    matrix[n_loc, n_loc] D_mat;

    // Fixed effects
    vector[n] eta_fixed;         // D * beta (without spatial effect)
    vector<lower=0, upper=1>[n] mda_impact; // MDA impact function

    // Fixed parameters
    real<lower=0> k;             // Aggregation parameter
    real<lower=0> rho;           // Detection rate
    real<lower=0> sigma2;        // Spatial variance
    real<lower=0> phi;           // Spatial range
  }

  parameters {
    vector[n_loc] S_raw;         // Standardized spatial process
  }

  transformed parameters {
    vector[n_loc] S;
    vector[n] mu_W_star;
    vector[n] mu_W;
    vector[n] pr0;

    // Correlation matrix
    matrix[n_loc, n_loc] R = compute_corr_matrix(D_mat, phi);
    matrix[n_loc, n_loc] L = cholesky_decompose(R);

    // Spatial process: S ~ N(0, sigma2 * R)
    S = sqrt(sigma2) * L * S_raw;

    // Mean worm burden
    for (i in 1:n) {
      mu_W_star[i] = exp(eta_fixed[i] + S[ID_coords[i]]);
      mu_W[i] = mu_W_star[i] * mda_impact[i];
      pr0[i] = compute_prevalence(mu_W[i], k, rho);

      // Clamp prevalence
      if (pr0[i] < 1e-10) pr0[i] = 1e-10;
      if (pr0[i] > 1 - 1e-10) pr0[i] = 1 - 1e-10;
    }
  }

  model {
    // Prior on standardized spatial process
    S_raw ~ normal(0, 1);

    // Likelihood for zeros
    for (i in 1:n) {
      if (y[i] == 0) {
        target += log(1 - pr0[i]);
      }
    }

    // Likelihood for positives
    for (i in 1:n_pos) {
      int idx = pos_idx[i];
      real mu_W_i = mu_W[idx];
      real pr0_i = pr0[idx];
      vector[2] gamma_params;

      // Prevalence contribution
      target += log(pr0_i);

      // Intensity contribution (shifted Gamma)
      gamma_params = compute_gamma_params(mu_W_i, pr0_i, k, rho);

      // C - 1 ~ Gamma(kappa, theta)
      if (gamma_params[1] > 0 && gamma_params[2] > 0) {
        target += gamma_lpdf(C_pos[i] - 1 | gamma_params[1], gamma_params[2]);
      }
    }
  }

  "

  # =============================================================================
  # COMPILE AND SAMPLE
  # =============================================================================

  if(messages) {
    message("Compiling STAN model...")
  }

  # Compile model
  stan_model <- rstan::stan_model(model_code = stan_code,
                                  model_name = "dsgm_spatial",
                                  verbose = FALSE)

  if(messages) {
    message(sprintf("Sampling %d iterations (%d warmup) across %d chains...",
                    n_samples + n_warmup, n_warmup, n_chains))
  }

  # Sample
  fit <- rstan::sampling(
    stan_model,
    data = stan_data,
    iter = n_samples + n_warmup,
    warmup = n_warmup,
    chains = n_chains,
    cores = n_cores,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    refresh = ifelse(messages, max(1, (n_samples + n_warmup) %/% 10), 0),
    show_messages = messages,
    verbose = messages
  )

  # =============================================================================
  # EXTRACT SAMPLES
  # =============================================================================

  # Extract S samples (spatial process at each location)
  S_samples <- rstan::extract(fit, pars = "S")$S

  # S_samples is now (n_samples_total x n_loc) where n_samples_total = n_samples * n_chains

  if(messages) {
    message(sprintf("Extracted %d samples of spatial process (%d locations)",
                    nrow(S_samples), ncol(S_samples)))

    # Diagnostics
    summary_S <- rstan::summary(fit, pars = "S")$summary
    n_divergent <- rstan::get_num_divergent(fit)
    n_max_treedepth <- rstan::get_num_max_treedepth(fit)

    message(sprintf("Divergent transitions: %d", sum(n_divergent)))
    message(sprintf("Max treedepth hits: %d", sum(n_max_treedepth)))

    if(sum(n_divergent) > 0) {
      warning("Divergent transitions detected. Consider increasing adapt_delta.")
    }
  }

  # =============================================================================
  # RETURN RESULTS
  # =============================================================================

  result <- list(
    S_samples = S_samples,           # Matrix: n_samples_total x n_loc
    stan_fit = fit,                  # Full STAN fit object
    n_samples = nrow(S_samples),
    n_loc = n_loc,
    coords = coords,
    par = par
  )

  class(result) <- "dsgm_spatial_samples"

  return(result)
}


##' @title Thin spatial process samples
##' @description Thin MCMC samples to reduce autocorrelation
##' @param spatial_samples Output from sample_spatial_process_stan
##' @param thin Thinning interval
##' @return Thinned samples
##' @keywords internal
thin_spatial_samples <- function(spatial_samples, thin = 10) {

  n_total <- nrow(spatial_samples$S_samples)
  keep_idx <- seq(1, n_total, by = thin)

  spatial_samples$S_samples <- spatial_samples$S_samples[keep_idx, , drop = FALSE]
  spatial_samples$n_samples <- nrow(spatial_samples$S_samples)

  return(spatial_samples)
}


##' @title Compute effective sample size for spatial process
##' @description Calculate ESS for each spatial location
##' @param spatial_samples Output from sample_spatial_process_stan
##' @return Vector of ESS values for each location
##' @keywords internal
compute_spatial_ess <- function(spatial_samples) {

  if(!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for ESS calculation")
  }

  n_loc <- spatial_samples$n_loc
  ess_vec <- numeric(n_loc)

  for(i in 1:n_loc) {
    mcmc_obj <- coda::as.mcmc(spatial_samples$S_samples[, i])
    ess_vec[i] <- coda::effectiveSize(mcmc_obj)
  }

  return(ess_vec)
}


##' @title Print method for spatial samples
##' @export
print.dsgm_spatial_samples <- function(x, ...) {
  cat("DSGM Spatial Process Samples\n")
  cat(sprintf("  Number of samples: %d\n", x$n_samples))
  cat(sprintf("  Number of locations: %d\n", x$n_loc))
  cat(sprintf("  Sample matrix dimensions: %d x %d\n",
              nrow(x$S_samples), ncol(x$S_samples)))

  if(requireNamespace("coda", quietly = TRUE)) {
    ess <- compute_spatial_ess(x)
    cat(sprintf("  Effective sample size: %.0f - %.0f (median: %.0f)\n",
                min(ess), max(ess), median(ess)))
  }

  invisible(x)
}

##' @title Initial values for DSGM model
##' @description Computes starting values for joint prevalence-intensity model
##'   by fitting simplified model without spatial random effects
##' @keywords internal
##' @title Initial values for DSGM model
##' @description Computes starting values for joint prevalence-intensity model
##'   by fitting simplified model without spatial random effects
##' @keywords internal
##' @title Initial values for DSGM model
##' @description Computes starting values for joint prevalence-intensity model
##'   by fitting simplified model without spatial random effects
##' @keywords internal
dsgm_initial_value <- function(y_prev, intensity_data, D, coords, ID_coords,
                               int_mat, survey_times_data, mda_times,
                               penalty, fix_alpha_W, fix_gamma_W, start_pars) {

  n <- length(y_prev)
  p <- ncol(D)
  n_loc <- nrow(coords)

  # =============================================================================
  # NEGATIVE LOG-LIKELIHOOD FOR SIMPLIFIED MODEL (NO SPATIAL EFFECTS)
  # =============================================================================

  nll_simplified <- function(par) {

    # Extract parameters (ONLY the fundamental parameters)
    beta <- par[1:p]                    # Covariate effects on log(mu_W*)
    k_log <- par[p + 1]                 # Aggregation parameter (log scale)
    rho_log <- par[p + 2]               # Detection rate (log scale)

    k <- exp(k_log)
    rho <- exp(rho_log)

    # MDA parameters
    idx <- p + 3
    if(is.null(fix_alpha_W)) {
      alpha_W_logit <- par[idx]
      alpha_W <- plogis(alpha_W_logit)
      idx <- idx + 1
    } else {
      alpha_W <- fix_alpha_W
    }

    if(is.null(fix_gamma_W)) {
      gamma_W_log <- par[idx]
      gamma_W <- exp(gamma_W_log)
      idx <- idx + 1
    } else {
      gamma_W <- fix_gamma_W
    }

    # Counterfactual mean worm burden (without MDA)
    mu_W_star <- exp(D %*% beta)

    # MDA impact function
    mda_imp <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                                   alpha_W, gamma_W, kappa = 1)

    # Adjusted mean worm burden
    mu_W <- mu_W_star * mda_imp

    # Compute prevalence using PGF of negative binomial
    # p(x,t) = 1 - [k / (k + mu_W * (1 - exp(-rho)))]^k
    term <- k / (k + mu_W * (1 - exp(-rho)))
    pr0 <- 1 - term^k

    # Clamp probabilities to avoid numerical issues
    pr0 <- pmax(pmin(pr0, 1 - 1e-10), 1e-10)

    # Initialize log-likelihood
    ll <- 0

    # Contribution from zeros (C_i = 0)
    zeros_idx <- which(y_prev == 0)
    if(length(zeros_idx) > 0) {
      ll <- ll + sum(log(1 - p[zeros_idx] + 1e-10))
    }

    # Contribution from positives (C_i > 0)
    positives_idx <- which(y_prev == 1)
    if(length(positives_idx) > 0) {
      # Prevalence contribution
      ll <- ll + sum(log(pr0[positives_idx] + 1e-10))

      # Extract mu_W and pr0 for positives
      mu_W_pos <- mu_W[positives_idx]
      p_pos <- pr0[positives_idx]

      # Conditional mean: mu_C = rho * mu_W / p
      mu_C <- (rho * mu_W_pos) / p_pos

      # Conditional variance (from tex file):
      # sigma^2_C = [rho * mu_W * (1 + rho)] / p +
      #             [rho^2 * mu_W^2 / p] * (1/k + 1 - 1/p)
      sigma2_C <- (rho * mu_W_pos * (1 + rho)) / p_pos +
        (rho^2 * mu_W_pos^2 / p_pos) * (1/k + 1 - 1/p_pos)

      # Ensure positive variance
      sigma2_C <- pmax(sigma2_C, 1e-10)

      # Gamma parameters (from moment matching):
      # kappa = (mu_C - 1)^2 / sigma^2_C
      # theta = sigma^2_C / (mu_C - 1)
      kappa_C <- ((mu_C - 1)^2) / sigma2_C
      theta_C <- sigma2_C / (mu_C - 1)

      # Ensure valid parameters
      kappa_C <- pmax(kappa_C, 1e-10)
      theta_C <- pmax(theta_C, 1e-10)

      # Extract observed intensities
      C_obs <- intensity_data  # Already filtered to positives

      # Gamma density for C - 1 (shifted)
      # C - 1 ~ Gamma(kappa_C, theta_C)
      ll_intensity <- dgamma(C_obs - 1, shape = kappa_C, scale = theta_C, log = TRUE)

      # Check for valid log-likelihood values
      ll_intensity[!is.finite(ll_intensity)] <- -1e10

      ll <- ll + sum(ll_intensity)
    }

    # Add penalties if specified
    if(!is.null(penalty)) {
      # Penalty on alpha_W (if not fixed)
      if(is.null(fix_alpha_W) && !is.null(penalty$alpha)) {
        ll <- ll - penalty$alpha(alpha_W)
      }

      # Penalty on gamma_W (if not fixed)
      if(is.null(fix_gamma_W) && !is.null(penalty$gamma)) {
        ll <- ll - penalty$gamma(gamma_W)
      }
    }

    # Return negative log-likelihood
    nll <- -ll

    # Check for numerical issues
    if(!is.finite(nll) || nll > 1e10) {
      return(1e10)
    }

    return(nll)
  }

  # =============================================================================
  # SET UP STARTING VALUES
  # =============================================================================

  # Observed statistics
  obs_prev <- mean(y_prev)
  obs_mean_epg <- mean(intensity_data, na.rm = TRUE)
  obs_var_epg <- var(intensity_data, na.rm = TRUE)
  if(is.na(obs_mean_epg)) obs_mean_epg <- 100
  if(is.na(obs_var_epg)) obs_var_epg <- obs_mean_epg^2

  # Starting values for beta (intercept-only model initially)
  if(is.null(start_pars$beta)) {
    # Back-calculate approximate mu_W from observed prevalence and intensity
    rho_approx <- 1  # Assume 1 egg per worm initially
    k_approx <- 1    # Initial aggregation

    # Rough approximation for mu_W
    mu_W_approx <- obs_prev * obs_mean_epg / rho_approx

    beta_init <- rep(0, p)
    beta_init[1] <- log(max(mu_W_approx, 1))  # Intercept
  } else {
    beta_init <- start_pars$beta
  }

  # Starting value for k (aggregation parameter)
  if(is.null(start_pars$k)) {
    # Estimate from overdispersion in intensity data
    # CV^2 for C - 1 gives rough estimate of 1/k
    mean_C_minus_1 <- mean(intensity_data - 1, na.rm = TRUE)
    var_C_minus_1 <- var(intensity_data - 1, na.rm = TRUE)

    if(var_C_minus_1 > 0 && mean_C_minus_1 > 0) {
      CV_sq <- var_C_minus_1 / (mean_C_minus_1^2)
      k_init <- max(1 / CV_sq, 0.1)  # Ensure reasonable value
    } else {
      k_init <- 0.5  # Default moderate aggregation
    }
  } else {
    k_init <- start_pars$k
  }

  # Starting value for rho (detection rate)
  if(is.null(start_pars$rho)) {
    rho_init <- 1  # 1 egg per worm
  } else {
    rho_init <- start_pars$rho
  }

  # Starting values for MDA parameters
  if(is.null(fix_alpha_W)) {
    if(is.null(start_pars$alpha_W)) {
      alpha_W_init <- 0.5  # 50% reduction
    } else {
      alpha_W_init <- start_pars$alpha_W
    }
  }

  if(is.null(fix_gamma_W)) {
    if(is.null(start_pars$gamma_W)) {
      gamma_W_init <- 2.0  # 2-year decay
    } else {
      gamma_W_init <- start_pars$gamma_W
    }
  }

  # =============================================================================
  # ASSEMBLE INITIAL PARAMETER VECTOR (ONLY FUNDAMENTAL PARAMETERS)
  # =============================================================================

  par_init <- c(
    beta_init,           # Covariate effects
    log(k_init),         # Log aggregation parameter
    log(rho_init)        # Log detection rate
  )

  # Add MDA parameters
  if(is.null(fix_alpha_W)) {
    par_init <- c(par_init, qlogis(alpha_W_init))
  }
  if(is.null(fix_gamma_W)) {
    par_init <- c(par_init, log(gamma_W_init))
  }

  # =============================================================================
  # OPTIMIZE
  # =============================================================================

  fit <- nlminb(
    start = par_init,
    objective = nll_simplified,
    control = list(eval.max = 2000, iter.max = 1000, trace = 0)
  )

  if(fit$convergence != 0) {
    warning("Initial value optimization did not converge (code = ",
            fit$convergence, "). Using starting values.")
    par_est <- par_init
  } else {
    par_est <- fit$par
  }

  # =============================================================================
  # EXTRACT ESTIMATED PARAMETERS
  # =============================================================================

  beta_est <- par_est[1:p]
  k_est <- exp(par_est[p + 1])
  rho_est <- exp(par_est[p + 2])

  idx <- p + 3
  if(is.null(fix_alpha_W)) {
    alpha_W_est <- plogis(par_est[idx])
    idx <- idx + 1
  } else {
    alpha_W_est <- fix_alpha_W
  }

  if(is.null(fix_gamma_W)) {
    gamma_W_est <- exp(par_est[idx])
    idx <- idx + 1
  } else {
    gamma_W_est <- fix_gamma_W
  }

  # =============================================================================
  # INITIALIZE SPATIAL PARAMETERS
  # =============================================================================

  # Compute residuals to initialize spatial variance
  mu_W_star_est <- exp(D %*% beta_est)
  mda_imp_est <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                                     alpha_W_est, gamma_W_est, kappa = 1)
  mu_W_est <- mu_W_star_est * mda_imp_est

  # Deviance residuals as proxy for spatial variation
  term_est <- k_est / (k_est + mu_W_est * (1 - exp(-rho_est)))
  p_est <- 1 - term_est^k_est
  p_est <- pmax(pmin(p_est, 1 - 1e-10), 1e-10)

  # Compute log-likelihood residuals
  resid <- numeric(n)
  resid[y_prev == 0] <- log(1 - p_est[y_prev == 0])
  resid[y_prev == 1] <- log(p_est[y_prev == 1])

  # Aggregate residuals by location
  resid_by_loc <- tapply(resid, ID_coords, mean, na.rm = TRUE)

  # Initialize spatial variance
  if(is.null(start_pars$sigma2)) {
    sigma2_init <- var(resid_by_loc, na.rm = TRUE)
    sigma2_init <- max(sigma2_init, 0.1)  # Ensure positive
  } else {
    sigma2_init <- start_pars$sigma2
  }

  # Initialize spatial range
  if(is.null(start_pars$phi)) {
    # Use median distance as initial range
    dist_mat <- as.matrix(dist(coords))
    median_dist <- median(dist_mat[upper.tri(dist_mat)])
    phi_init <- median_dist / 3  # Practical range ~ 3*phi
  } else {
    phi_init <- start_pars$phi
  }

  # =============================================================================
  # RETURN INITIAL PARAMETER LIST
  # =============================================================================

  par0 <- list(
    beta = beta_est,
    k = k_est,
    rho = rho_est,
    alpha_W = alpha_W_est,
    gamma_W = gamma_W_est,
    sigma2 = sigma2_init,
    phi = phi_init
  )

  return(par0)
}

##' @title Fitting of Doubly Stochastic Geostatistical Model (DSGM) for joint analysis of prevalence
##' and intensity of infection
##'
##' @description
##' The function fits a joint prevalence-intensity model for STH using Monte Carlo maximum likelihood.
##' The model incorporates a mechanistic worm burden framework where MDA effects decay over time,
##' affecting both prevalence and intensity through their impact on mean worm burden.
##' Uses STAN for sampling the spatial process conditional on fixed parameters.
##'
##' Spatial structure:
##' \itemize{
##'   \item \code{gp(x, y)} fits a spatial Gaussian process S(x).
##'   \item The spatial process affects mean worm burden: log(mu_W*) = D*beta + S(x)
##' }
##'
##' @param formula A model formula with egg counts (EPG) as response and covariates for mean worm burden.
##'   Example: \code{epg ~ elevation + rainfall + gp(x, y)} where epg is the egg count variable.
##'   Zero values indicate uninfected individuals.
##' @param data A \code{data.frame} or \code{sf} object containing the dataset.
##' @param time A variable in \code{data} giving the survey times of observations (required).
##' @param mda_times A vector specifying the mass drug administration (MDA) times.
##' @param int_mat Intervention matrix specifying MDA timing and coverage; dimension \code{n * n_mda}.
##' @param penalty Optional list specifying penalty functions for regularization of MDA parameters.
##' @param drop_W Optional fixed value for the immediate worm burden reduction parameter (alpha_W).
##' @param decay_W Optional fixed value for the worm burden decay parameter (gamma_W).
##' @param crs Optional coordinate reference system (CRS) for spatial data.
##' @param convert_to_crs CRS to which spatial data should be converted.
##' @param scale_to_km Logical; whether to scale distances to kilometers (default: \code{TRUE}).
##' @param par0 Optional list of initial parameter values.
##' @param n_samples Number of MCMC samples to draw from spatial process (default: 1000).
##' @param n_warmup Number of warmup/burnin iterations for STAN (default: 1000).
##' @param n_chains Number of MCMC chains for STAN (default: 1).
##' @param adapt_delta Target acceptance rate for STAN sampler (default: 0.8).
##' @param max_treedepth Maximum tree depth for STAN sampler (default: 10).
##' @param return_samples Logical; whether to return spatial process samples (default: \code{TRUE}).
##' @param messages Logical; whether to print messages (default: \code{TRUE}).
##' @param start_pars List of starting values for parameters including:
##'   \itemize{
##'     \item beta: Covariate effects on log mean worm burden
##'     \item k: Negative binomial aggregation parameter
##'     \item rho: Egg detection rate per worm
##'     \item sigma2: Variance of the spatial process
##'     \item phi: Scale parameter of the spatial correlation
##'     \item alpha_W: Immediate worm burden reduction
##'     \item gamma_W: Worm burden decay rate
##'   }
##'
##' @return A list containing:
##' \itemize{
##'   \item \code{prevalence_data}: Binary infection indicators (C > 0).
##'   \item \code{intensity_data}: Egg counts among positive individuals.
##'   \item \code{n_positive}: Number of egg-positive individuals per location.
##'   \item \code{D}: Covariate matrix.
##'   \item \code{coords}: Unique spatial coordinates.
##'   \item \code{mda_times}: MDA time points.
##'   \item \code{survey_times_data}: Survey times.
##'   \item \code{int_mat}: Intervention matrix.
##'   \item \code{ID_coords}: Spatial location indices.
##'   \item \code{fix_alpha_W}, \code{fix_gamma_W}: Fixed MDA parameters if specified.
##'   \item \code{formula}: Model formula.
##'   \item \code{model_params}: Estimated parameters for worm burden, detection, and spatial processes.
##'   \item \code{spatial_samples}: STAN spatial process samples if \code{return_samples = TRUE}.
##'   \item \code{family}: "zero_inflated_gamma" indicating joint model structure.
##' }
##'
##' @seealso \code{\link{dast}}
##' @author [Your name]
##' @export
dsgm <- function(formula,
                 data,
                 time,
                 mda_times,
                 int_mat,
                 penalty = NULL,
                 drop_W = NULL,
                 decay_W = NULL,
                 crs = NULL,
                 convert_to_crs = NULL,
                 scale_to_km = TRUE,
                 par0 = NULL,
                 n_samples = 1000,
                 n_warmup = 1000,
                 n_chains = 1,
                 adapt_delta = 0.8,
                 max_treedepth = 10,
                 return_samples = TRUE,
                 messages = TRUE,
                 start_pars = list(beta = NULL,
                                   k = NULL,
                                   rho = NULL,
                                   sigma2 = NULL,
                                   phi = NULL,
                                   alpha_W = NULL,
                                   gamma_W = NULL)) {

  # =============================================================================
  # ARGUMENT VALIDATION
  # =============================================================================

  # Check formula
  if(!inherits(formula, what = "formula", which = FALSE)) {
    stop("'formula' must be a 'formula' object with egg count as response and covariates for mean worm burden")
  }

  # Check data
  if(!inherits(data, c("data.frame", "sf"))) {
    stop("'data' must be a data.frame or sf object")
  }

  # Check STAN parameters
  if(n_samples <= 0 || n_warmup < 0) {
    stop("'n_samples' must be positive and 'n_warmup' must be non-negative")
  }
  if(n_chains < 1) {
    stop("'n_chains' must be at least 1")
  }
  if(adapt_delta <= 0 || adapt_delta >= 1) {
    stop("'adapt_delta' must be in (0, 1)")
  }

  # =============================================================================
  # EXTRACT AND VALIDATE TIME VARIABLE
  # =============================================================================

  time_name <- deparse(substitute(time))
  if (time_name == "NULL") {
    stop("You must supply a time column via the 'time' argument")
  }
  if (!time_name %in% names(data)) {
    stop("'time' column not found in 'data'")
  }
  survey_times_data <- data[[time_name]]

  if(any(is.na(survey_times_data))) {
    stop("Missing values detected in 'time' variable")
  }

  # =============================================================================
  # PROCESS FORMULA AND EXTRACT RESPONSE (EGG COUNT) AND COVARIATES
  # =============================================================================

  inter_f <- interpret.formula(formula)

  # Extract GP specification
  gp_terms <- inter_f$gp.spec$term
  gp_dim   <- inter_f$gp.spec$dim

  # Currently only support spatial GP (no spatio-temporal for joint model)
  if(length(gp_terms) == 1 && gp_terms[1] == "sf") {
    # Using sf geometry
  } else if(gp_dim == 2) {
    # Using gp(x, y)
  } else {
    stop("For joint model, specify gp(x, y) or gp(sf). Spatio-temporal GP not yet supported for dsgm")
  }

  # Extract model frame
  mf <- model.frame(inter_f$pf, data = data, na.action = na.fail)

  # Extract response variable (egg counts)
  egg_counts <- as.numeric(model.response(mf))
  n <- length(egg_counts)

  # Validate egg counts
  if(!is.numeric(egg_counts)) {
    stop("Response variable (egg count) must be numeric")
  }
  if(any(egg_counts < 0, na.rm = TRUE)) {
    stop("Response variable (egg count) cannot contain negative values")
  }
  if(any(is.na(egg_counts))) {
    stop("Missing values detected in response variable (egg count)")
  }

  # Get response variable name for messages
  egg_count_name <- as.character(formula[[2]])

  # Create binary infection indicator
  y_prev <- as.integer(egg_counts > 0)

  # Extract intensity data (only positive counts)
  intensity_data <- egg_counts[egg_counts > 0]
  n_positive_total <- sum(y_prev)

  if(n_positive_total == 0) {
    stop("No positive egg counts detected; model cannot be fitted")
  }

  if(messages) {
    message(sprintf("Response variable: %s", egg_count_name))
    message(sprintf("Data summary: %d observations, %d (%.1f%%) egg-positive",
                    length(egg_counts), n_positive_total,
                    100 * n_positive_total / length(egg_counts)))
  }

  # Extract covariate matrix (for mean worm burden)
  D <- as.matrix(model.matrix(attr(mf, "terms"), data = data))
  p <- ncol(D)

  # Handle offset
  if(is.null(inter_f$offset)) {
    cov_offset <- rep(0, n)
  } else {
    cov_offset <- data[[inter_f$offset]]
  }

  # =============================================================================
  # VALIDATE MDA DATA
  # =============================================================================

  if(missing(mda_times)) {
    stop("'mda_times' must be specified")
  }
  if(missing(int_mat)) {
    stop("'int_mat' (intervention matrix) must be specified")
  }

  if(!is.numeric(mda_times)) {
    stop("'mda_times' must be numeric")
  }

  if(!is.matrix(int_mat) && !is.data.frame(int_mat)) {
    stop("'int_mat' must be a matrix or data.frame")
  }

  int_mat <- as.matrix(int_mat)

  if(nrow(int_mat) != nrow(data)) {
    stop(sprintf("'int_mat' must have %d rows to match 'data'", nrow(data)))
  }
  if(ncol(int_mat) != length(mda_times)) {
    stop(sprintf("'int_mat' must have %d columns to match length of 'mda_times'",
                 length(mda_times)))
  }

  if(any(int_mat < 0 | int_mat > 1, na.rm = TRUE)) {
    warning("'int_mat' values outside [0,1] detected; these will be treated as coverage proportions")
  }

  # =============================================================================
  # PROCESS SPATIAL COORDINATES
  # =============================================================================

  if(inherits(data, "sf")) {
    # Convert to specified CRS if needed
    if(!is.null(convert_to_crs)) {
      data <- sf::st_transform(data, convert_to_crs)
      crs <- convert_to_crs
    } else {
      crs <- sf::st_crs(data)$input
    }

    coords_all <- sf::st_coordinates(data)
  } else {
    # Extract coordinates from gp() specification
    if(gp_terms[1] != "sf") {
      coord_names <- gp_terms[1:2]
      if(!all(coord_names %in% names(data))) {
        stop("Coordinate columns specified in gp() not found in 'data'")
      }
      coords_all <- as.matrix(data[, coord_names])
    } else {
      stop("If using gp(sf), data must be an sf object")
    }
  }

  # Scale to kilometers if requested
  if(scale_to_km) {
    coords_all <- coords_all / 1000
    if(messages) {
      message("Coordinates scaled to kilometers")
    }
  }

  # Get unique spatial locations
  coords_unique <- unique(coords_all)
  n_loc <- nrow(coords_unique)

  # Create location indices
  ID_coords <- apply(coords_all, 1, function(x) {
    which(apply(coords_unique, 1, function(y) all(abs(x - y) < 1e-10)))[1]
  })

  if(messages) {
    message(sprintf("Identified %d unique spatial locations", n_loc))
  }

  # =============================================================================
  # HANDLE FIXED PARAMETERS
  # =============================================================================

  fix_alpha_W <- NULL
  fix_gamma_W <- NULL

  if(!is.null(drop_W)) {
    if(drop_W < 0 || drop_W > 1) {
      stop("'drop_W' must be in [0, 1]")
    }
    fix_alpha_W <- drop_W
    if(messages) {
      message(sprintf("MDA worm burden reduction (alpha_W) fixed at %.3f", fix_alpha_W))
    }
  }

  if(!is.null(decay_W)) {
    if(decay_W <= 0) {
      stop("'decay_W' must be positive")
    }
    fix_gamma_W <- decay_W
    if(messages) {
      message(sprintf("MDA decay rate (gamma_W) fixed at %.3f years", fix_gamma_W))
    }
  }

  # =============================================================================
  # HANDLE PENALTY FUNCTIONS
  # =============================================================================

  no_penalty <- TRUE
  if(!is.null(penalty)) {
    if(!is.list(penalty)) {
      stop("'penalty' must be a list")
    }
    no_penalty <- FALSE
  }

  # =============================================================================
  # SET UP INITIAL PARAMETERS
  # =============================================================================

  if(is.null(par0)) {
    if(messages) {
      message("\n=== STEP 1: Computing initial parameter values ===")
    }

    par0 <- dsgm_initial_value(
      y_prev = y_prev,
      intensity_data = intensity_data,
      D = D,
      coords = coords_unique,
      ID_coords = ID_coords,
      int_mat = int_mat,
      survey_times_data = survey_times_data,
      mda_times = mda_times,
      penalty = penalty,
      fix_alpha_W = fix_alpha_W,
      fix_gamma_W = fix_gamma_W,
      start_pars = start_pars
    )

    if(messages) {
      message("\nInitial parameter values computed successfully")
    }
  }

  # Validate par0 structure
  required_pars <- c("beta", "k", "rho", "sigma2", "phi")
  if(is.null(fix_alpha_W)) required_pars <- c(required_pars, "alpha_W")
  if(is.null(fix_gamma_W)) required_pars <- c(required_pars, "gamma_W")

  missing_pars <- setdiff(required_pars, names(par0))
  if(length(missing_pars) > 0) {
    stop(sprintf("Missing initial parameters: %s", paste(missing_pars, collapse = ", ")))
  }

  # =============================================================================
  # SAMPLE SPATIAL PROCESS USING STAN
  # =============================================================================

  if(messages) {
    message("\n=== STEP 2: Sampling spatial process using STAN ===")
    message(sprintf("STAN settings:"))
    message(sprintf("  - Samples: %d", n_samples))
    message(sprintf("  - Warmup: %d", n_warmup))
    message(sprintf("  - Chains: %d", n_chains))
    message(sprintf("  - adapt_delta: %.2f", adapt_delta))
    message(sprintf("  - max_treedepth: %d", max_treedepth))
  }

  spatial_samples <- sample_spatial_process_stan(
    y_prev = y_prev,
    intensity_data = intensity_data,
    D = D,
    coords = coords_unique,
    ID_coords = ID_coords,
    int_mat = int_mat,
    survey_times_data = survey_times_data,
    mda_times = mda_times,
    par = par0,
    n_samples = n_samples,
    n_warmup = n_warmup,
    n_chains = n_chains,
    n_cores = 1,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    messages = messages
  )

  dsgm_res <- dsg_fit()

  if(messages) {
    message("\nSpatial process sampling complete")
    message(sprintf("Obtained %d samples for %d locations",
                    spatial_samples$n_samples, spatial_samples$n_loc))
  }

  # =============================================================================
  # FIT THE MODEL (placeholder - to be implemented)
  # =============================================================================

  if(messages) {
    message("\n=== STEP 3: Fitting joint prevalence-intensity model via MCML ===")
    message(sprintf("  - %d observations at %d locations", n, n_loc))
    message(sprintf("  - %d MDA rounds", length(mda_times)))
    message(sprintf("  - %d covariates", p))
  }

  # This would call the MCML optimization using the spatial samples
  # fit_result <- dsgm_fit(...)
  # For now, return a placeholder
  fit_result <- list(
    params = par0,
    convergence = 0,
    log_likelihood = NA
  )

  if(messages) {
    message("Model fitting not yet implemented (placeholder)")
  }

  # =============================================================================
  # PREPARE RETURN OBJECT
  # =============================================================================

  res <- list()
  res$prevalence_data <- y_prev
  res$intensity_data <- intensity_data
  res$egg_counts <- egg_counts
  res$n_positive <- sum(y_prev)
  res$D <- D
  res$coords <- coords_unique
  res$coords_all <- coords_all
  res$mda_times <- mda_times
  res$survey_times_data <- survey_times_data
  res$int_mat <- int_mat
  res$ID_coords <- ID_coords
  res$fix_alpha_W <- fix_alpha_W
  res$fix_gamma_W <- fix_gamma_W
  res$formula <- formula
  res$crs <- crs
  res$scale_to_km <- scale_to_km
  res$data_sf <- data
  res$family <- "zero_inflated_gamma"
  res$cov_offset <- cov_offset
  res$call <- match.call()

  if(!no_penalty) {
    res$penalty <- penalty
  } else {
    res$penalty <- NULL
  }

  # Add fitted results
  res$model_params <- fit_result$params
  res$convergence <- fit_result$convergence
  res$log_likelihood <- fit_result$log_likelihood

  # Add spatial samples
  if(return_samples) {
    res$spatial_samples <- spatial_samples
  }

  # Add model structure info
  res$n_locations <- n_loc
  res$n_observations <- n
  res$n_covariates <- p
  res$n_mda_rounds <- length(mda_times)

  # Add STAN settings
  res$stan_settings <- list(
    n_samples = n_samples,
    n_warmup = n_warmup,
    n_chains = n_chains,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )

  class(res) <- c("dsgm", "RiskMap")

  if(messages) {
    message("\n=== Model fitting complete! ===")
  }

  return(res)
}
