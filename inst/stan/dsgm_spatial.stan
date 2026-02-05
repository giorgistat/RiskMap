functions {
  // Compute prevalence from mean worm burden
  real compute_prevalence(real mu_W, real k, real rho) {
    real term = k / (k + mu_W * (1.0 - exp( -rho )));
    real pr = 1.0 - pow(term, k);

    // Clamp prevalence
    if (pr < 1e-10) pr = 1e-10;
    if (pr > 1.0 - 1e-10) pr = 1.0 - 1e-10;

    return pr;
  }

  // Compute conditional mean intensity
  real compute_mu_C(real mu_W, real pr, real rho) {
    return (rho * mu_W) / pr;
  }

  // Compute conditional variance
  real compute_sigma2_C(real mu_W, real pr, real k, real rho) {
    real var1 = (rho * mu_W * (1.0+ rho)) / pr;
    real var2 = (square(rho) * square(mu_W) / pr) * (1.0/k + 1.0 - 1.0/pr);
    real result = var1 + var2;

    // Safeguard against numerical issues
    if (result < 1e-10) result = 1e-10;

    return result;
  }

  // Compute shifted Gamma parameters
  // CRITICAL FIX: Returns [shape, RATE] not [shape, scale]
  vector compute_gamma_params(real mu_W, real pr, real k, real rho) {
    vector[2] params;
    real mu_C = compute_mu_C(mu_W, pr, rho);
    real sigma2_C = compute_sigma2_C(mu_W, pr, k, rho);

    // Ensure mu_C > 1 and sigma2_C > 0
    if (mu_C <= 1.0) mu_C = 1.0 + 1e-10;
    if (sigma2_C < 1e-10) sigma2_C = 1e-10;

    // kappa = (mu_C - 1)^2 / sigma2_C (shape)
    // lambda = (mu_C - 1) / sigma2_C (RATE, not scale!)
    // Note: STAN's gamma_lpdf uses (shape, rate) parameterization
    params[1] = square(mu_C - 1.0) / sigma2_C;  // kappa (shape)
    params[2] = (mu_C - 1.0) / sigma2_C;         // lambda (rate = 1/theta)

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

transformed data {
  // PERFORMANCE FIX: Compute correlation matrix ONCE (not every iteration!)
  // Since phi is fixed (data), R and L are constants
  matrix[n_loc, n_loc] R;
  matrix[n_loc, n_loc] L;

  // Build correlation matrix using exponential kernel
  for (i in 1:n_loc) {
    R[i, i] = 1.0;
    for (j in (i+1):n_loc) {
      R[i, j] = exp(-D_mat[i, j] / phi);
      R[j, i] = R[i, j];
    }
  }

  // Cholesky decomposition
  L = cholesky_decompose(R);
}

parameters {
  vector[n_loc] S_raw;         // Standardized spatial process
}

transformed parameters {
  vector[n_loc] S;
  vector[n] mu_W_star;
  vector[n] mu_W;
  vector[n] pr0;

  // Spatial process: S ~ N(0, sigma2 * R)
  // Using pre-computed L from transformed data
  S = sqrt(sigma2) * L * S_raw;

  // Mean worm burden
  for (i in 1:n) {
    mu_W_star[i] = exp(eta_fixed[i] + S[ID_coords[i]]);
    mu_W[i] = mu_W_star[i] * mda_impact[i];
    pr0[i] = compute_prevalence(mu_W[i], k, rho);
  }
}

model {
  // Prior on standardized spatial process
  S_raw ~ normal(0, 1);

  // Likelihood for zeros
  for (i in 1:n) {
    if (y[i] == 0) {
      target += log(1.0 - pr0[i]);
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

    // C - 1 ~ Gamma(kappa, lambda) where lambda is RATE
    if (gamma_params[1] > 0 && gamma_params[2] > 0) {
      target += gamma_lpdf(C_pos[i] - 1.0 | gamma_params[1], gamma_params[2]);
    }
  }
}

