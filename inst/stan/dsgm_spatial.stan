// =============================================================================
// dsgm_spatial.stan
//
// Stan sampler for the doubly stochastic geostatistical model for joint
// prevalence-intensity modelling of soil-transmitted helminths (STH).
// Used to draw MCMC samples of the spatial process S(x) at fixed parameter
// values theta_0, which are then passed to the TMB MCML engine (dsgm_tmb.cpp).
//
// Model:
//   W_j(x_i) | S(x_i) ~ NegBin(mu(x_i), omega_i)
//     log mu(x_i) = eta_fixed[i] + S[ID_coords[i]]
//     mu(x_i)     = exp(eta_fixed[i] + S[ID_coords[i]]) * mda_impact[i]
//
//   Aggregation parameter (controlled by vary_k):
//     vary_k = 0:  omega_i = exp(log_k)                          [constant]
//     vary_k = 1:  log(omega_i) = log_k + omega1 * log(mu_W[i]) [varies with burden]
//
//   Y_ij | W_j(x_i) ~ Poisson(rho * W_j(x_i))
//
// Intensity likelihood for C | C > 0 (controlled by intensity_family):
//   0 = Shifted Gamma:           C - 1 ~ Gamma(kappa_C, rate_C)
//   1 = Zero-truncated NegBin:   C | C > 0 ~ NegBin2(mu_C, phi_C) / (1 - p0)
//
// Spatial process:
//   S = sqrt(sigma2) * L * S_raw,   S_raw ~ N(0, I)
//   R[i,j] = exp(-D_mat[i,j] / phi),  L = cholesky(R)
//   L is computed once in transformed_data since phi is fixed.
//
// All model parameters are FIXED at their theta_0 values and passed as data.
// Only S_raw is sampled.
// =============================================================================

data {

  int<lower=1> n;        // Total number of observations
  int<lower=1> n_loc;    // Number of unique spatial locations
  int<lower=0> n_pos;    // Number of egg-positive observations
  int<lower=1> p;        // Number of fixed-effect covariates (unused; eta pre-computed)

  // Binary response
  array[n] int<lower=0, upper=1> y;

  // EPG for positive observations
  vector<lower=1>[n_pos]          C_pos;      // real, used by Gamma branch
  array[n_pos] int<lower=1>       C_pos_int;  // integer, used by NegBin branch

  // Indices of positive observations (1-indexed, matching y)
  array[n_pos] int<lower=1, upper=n> pos_idx;

  // Location mapping (1-indexed)
  array[n] int<lower=1, upper=n_loc> ID_coords;

  // Pairwise distance matrix between unique locations
  matrix[n_loc, n_loc] D_mat;

  // Pre-computed fixed-effects linear predictor: D * beta + offset
  vector[n] eta_fixed;

  // Pre-computed MDA decay factor per observation (1.0 = no MDA)
  vector<lower=0>[n] mda_impact;

  // Fixed parameters at theta_0
  // log_k  = log of the baseline aggregation parameter (= log(k) when vary_k=0)
  // omega1 = slope of log(omega) on log(mu_W); set to 0.0 when vary_k=0
  real log_k;
  real omega1;
  real<lower=0> rho;     // Per-worm egg detection rate
  real<lower=0> sigma2;  // GP variance
  real<lower=0> phi;     // GP range

  // Aggregation parameterisation flag
  // 0 = constant omega (original behaviour)
  // 1 = omega_i = exp(log_k + omega1 * log(mu_W_i))
  int<lower=0, upper=1> vary_k;

  // Intensity likelihood family: 0 = shifted Gamma, 1 = zero-truncated NegBin
  int<lower=0, upper=1> intensity_family;

}

transformed data {

  matrix[n_loc, n_loc] R;
  matrix[n_loc, n_loc] L;
  real c_rho = 1.0 - exp(-rho);

  for (i in 1:n_loc) {
    R[i, i] = 1.0;
    for (j in (i + 1):n_loc) {
      R[i, j] = exp(-D_mat[i, j] / phi);
      R[j, i] = R[i, j];
    }
  }
  L = cholesky_decompose(R);

}

parameters {

  vector[n_loc] S_raw;

}

transformed parameters {

  vector[n_loc] S;
  vector[n]     mu_W;
  vector[n]     k_vec;   // per-observation aggregation parameter
  vector[n]     pr_pos;

  S = sqrt(sigma2) * (L * S_raw);

  for (i in 1:n) {
    mu_W[i] = exp(eta_fixed[i] + S[ID_coords[i]]) * mda_impact[i];

    // Compute per-observation omega
    if (vary_k == 1) {
      k_vec[i] = exp(log_k + omega1 * log(fmax(mu_W[i], 1e-10)));
    } else {
      k_vec[i] = exp(log_k);
    }

    real ratio = k_vec[i] / (k_vec[i] + mu_W[i] * c_rho);
    pr_pos[i]  = fmax(fmin(1.0 - pow(ratio, k_vec[i]), 1.0 - 1e-10), 1e-10);
  }

}

model {

  // Prior on standardised spatial process
  S_raw ~ std_normal();

  // ----- Prevalence likelihood -----
  for (i in 1:n) {
    if (y[i] == 1)
      target += log(pr_pos[i]);
    else
      target += log1m(pr_pos[i]);
  }

  // ----- Intensity likelihood (positives only) -----
  for (idx in 1:n_pos) {
    int  i   = pos_idx[idx];
    real k_i = k_vec[i];

    // Conditional moments of C | C > 0
    real mu_C     = (rho * mu_W[i]) / pr_pos[i];
    real sigma2_C = fmax(
      (rho * mu_W[i] * (1.0 + rho)) / pr_pos[i]
      + (rho^2 * mu_W[i]^2 / pr_pos[i]) * (1.0/k_i + 1.0 - 1.0/pr_pos[i]),
      1e-6);

    if (intensity_family == 0) {

      // ------------------------------------------------------------------
      // Shifted Gamma: C - 1 ~ Gamma(kappa_C, rate_C)
      // ------------------------------------------------------------------
      real mu_C1   = fmax(mu_C - 1.0, 0.1);
      real kappa_C = fmax(square(mu_C1) / sigma2_C, 1e-6);
      real rate_C  = fmax(mu_C1 / sigma2_C,          1e-6);
      target += gamma_lpdf(C_pos[idx] - 1.0 | kappa_C, rate_C);

    } else {

      // ------------------------------------------------------------------
      // Zero-truncated NegBin2: C | C > 0 ~ NegBin2(mu_C, phi_C) / (1 - p0)
      // ------------------------------------------------------------------
      real denom_nb = fmax(sigma2_C - mu_C, 1e-6);
      real phi_C    = fmax(square(mu_C) / denom_nb, 1e-4);
      real log_p0   = neg_binomial_2_lpmf(0 | mu_C, phi_C);
      target += neg_binomial_2_lpmf(C_pos_int[idx] | mu_C, phi_C)
              - log1m_exp(log_p0);

    }
  }

}
