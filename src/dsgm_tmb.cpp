#include <TMB.hpp>

// Helper function to get distance from compressed vector
template<class Type>
Type get_distance(int i, int j, const vector<Type>& dist_vec, int n_loc) {
  if(i == j) return Type(0);
  if(i > j) {
    int temp = i;
    i = j;
    j = temp;
  }
  // Convert (i,j) to index in upper triangle stored as vector
  int idx = i * n_loc - (i * (i + 1)) / 2 + (j - i - 1);
  return dist_vec(idx);
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  // =============================================================================
  // DATA
  // =============================================================================

  DATA_VECTOR(y_prev);
  DATA_VECTOR(intensity_data);
  DATA_IVECTOR(pos_idx);
  DATA_MATRIX(D);
  DATA_VECTOR(cov_offset);
  DATA_MATRIX(S_samples);
  DATA_IVECTOR(ID_coords);

  // OPTIMIZATION 2: Compressed distance matrix
  DATA_VECTOR(dist_vec);
  DATA_INTEGER(n_loc);

  DATA_VECTOR(survey_times);
  DATA_VECTOR(mda_times);

  // OPTIMIZATION 1: Sparse MDA representation
  DATA_IVECTOR(mda_i);          // Row indices (0-indexed)
  DATA_IVECTOR(mda_j);          // Column indices (0-indexed)
  DATA_VECTOR(mda_coverage);    // Coverage values
  DATA_INTEGER(n_mda_pairs);    // Number of non-zero entries

  // Penalty configuration
  DATA_INTEGER(use_alpha_penalty);
  DATA_INTEGER(alpha_penalty_type);
  DATA_SCALAR(alpha_param1);
  DATA_SCALAR(alpha_param2);

  DATA_INTEGER(use_gamma_penalty);
  DATA_INTEGER(gamma_penalty_type);
  DATA_SCALAR(gamma_param1);
  DATA_SCALAR(gamma_param2);

  // OPTION to compute denominator only (for importance sampling)
  DATA_INTEGER(compute_denominator_only);
  DATA_VECTOR(log_denominator_vals);

  // =============================================================================
  // PARAMETERS (θ - current parameters being optimized)
  // =============================================================================

  PARAMETER_VECTOR(beta);
  PARAMETER(log_k);
  PARAMETER(log_rho);
  PARAMETER(logit_alpha);
  PARAMETER(log_gamma);
  PARAMETER(log_sigma2);
  PARAMETER(log_phi);

  // Transform to natural scale
  const Type k = exp(log_k);
  const Type rho = exp(log_rho);
  const Type alpha_W = invlogit(logit_alpha);
  const Type gamma_W = exp(log_gamma);
  const Type sigma2 = exp(log_sigma2);
  const Type phi = exp(log_phi);

  // Dimensions
  int n = y_prev.size();
  int n_samples = S_samples.rows();
  int n_pos = pos_idx.size();

  // =============================================================================
  // OPTIMIZATION 1: EFFICIENT MDA EFFECT COMPUTATION
  // =============================================================================

  vector<Type> mda_effect(n);
  mda_effect.setOnes();

  // Only loop over non-zero entries in int_mat
  for(int idx = 0; idx < n_mda_pairs; idx++) {
    int i = mda_i(idx);
    int m = mda_j(idx);
    Type coverage = mda_coverage(idx);

    Type time_since = survey_times(i) - mda_times(m);
    if(time_since > Type(0)) {
      Type decay = exp(-time_since / gamma_W);
      mda_effect(i) *= (Type(1.0) - alpha_W * coverage * decay);
    }
  }

  // =============================================================================
  // OPTIMIZATION 2: BUILD CORRELATION MATRIX FROM COMPRESSED DISTANCES
  // =============================================================================

  matrix<Type> R(n_loc, n_loc);
  for(int i = 0; i < n_loc; i++) {
    R(i, i) = Type(1.0);
    for(int j = i + 1; j < n_loc; j++) {
      Type dist_ij = get_distance(i, j, dist_vec, n_loc);
      R(i, j) = exp(-dist_ij / phi);
      R(j, i) = R(i, j);
    }
  }

  // Cholesky decomposition: R = L * L'
  Eigen::LLT<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> llt(R);
  matrix<Type> L = llt.matrixL();

  // Log-determinant via Cholesky
  Type log_det_R = Type(0.0);
  for(int i = 0; i < n_loc; i++) {
    log_det_R += log(L(i, i));
  }
  log_det_R *= Type(2.0);

  // =============================================================================
  // MONTE CARLO INTEGRATION WITH IMPORTANCE SAMPLING
  // =============================================================================

  vector<Type> log_f_vals(n_samples);

  vector<Type> mu = D * beta + cov_offset;
  Type k_term = Type(1.0) - exp(-rho);

  for(int s = 0; s < n_samples; s++) {
    vector<Type> S = S_samples.row(s);

    // =========================================================================
    // COMPUTE NUMERATOR: log f(y | S, θ) + log f(S | θ)
    // =========================================================================

    vector<Type> eta(n);
    for(int i = 0; i < n; i++) {
      eta(i) = mu(i) + S(ID_coords(i));
    }
    vector<Type> mu_W_star = exp(eta);
    vector<Type> mu_W = mu_W_star * mda_effect;

    // Prevalence via PGF
    vector<Type> denom = k + mu_W * k_term;
    vector<Type> ratio = k / denom;
    vector<Type> pr = Type(1.0) - pow(ratio, k);

    // Clamp prevalence
    for(int i = 0; i < n; i++) {
      pr(i) = CppAD::CondExpLt(pr(i), Type(1e-10), Type(1e-10), pr(i));
      pr(i) = CppAD::CondExpGt(pr(i), Type(1.0 - 1e-10), Type(1.0 - 1e-10), pr(i));
    }

    Type ll_theta = Type(0.0);

    // Zeros contribution
    for(int i = 0; i < n; i++) {
      if(y_prev(i) == 0) {
        ll_theta += log(Type(1.0) - pr(i));
      }
    }

    // Positives prevalence contribution
    for(int j = 0; j < n_pos; j++) {
      int i = pos_idx(j);
      ll_theta += log(pr(i));
    }

    // =======================================================================
    // OPTIMIZATION 3: VECTORIZED POSITIVE INTENSITY COMPUTATION
    // =======================================================================

    // Pre-compute vectors for all positives
    vector<Type> mu_C_vec(n_pos);
    vector<Type> sigma2_C_vec(n_pos);
    vector<Type> kappa_C_vec(n_pos);
    vector<Type> theta_C_vec(n_pos);

    for(int j = 0; j < n_pos; j++) {
      int i = pos_idx(j);
      Type mu_W_i = mu_W(i);
      Type pr_i = pr(i);

      // Conditional mean
      mu_C_vec(j) = (rho * mu_W_i) / pr_i;

      // Conditional variance
      Type term1 = (rho * mu_W_i * (Type(1.0) + rho)) / pr_i;
      Type rho_sq_mu_sq = rho * rho * mu_W_i * mu_W_i;
      Type term2 = (rho_sq_mu_sq / pr_i) * (Type(1.0)/k + Type(1.0) - Type(1.0)/pr_i);
      sigma2_C_vec(j) = term1 + term2;

      // Clamp variance
      sigma2_C_vec(j) = CppAD::CondExpLt(sigma2_C_vec(j), Type(1e-10),
                   Type(1e-10), sigma2_C_vec(j));

      // Gamma parameters
      Type mu_C_minus_1 = mu_C_vec(j) - Type(1.0);
      kappa_C_vec(j) = mu_C_minus_1 * mu_C_minus_1 / sigma2_C_vec(j);
      theta_C_vec(j) = sigma2_C_vec(j) / mu_C_minus_1;

      // Clamp Gamma parameters
      kappa_C_vec(j) = CppAD::CondExpLt(kappa_C_vec(j), Type(1e-10),
                  Type(1e-10), kappa_C_vec(j));
      theta_C_vec(j) = CppAD::CondExpLt(theta_C_vec(j), Type(1e-10),
                  Type(1e-10), theta_C_vec(j));
    }

    // Add intensity contributions
    for(int j = 0; j < n_pos; j++) {
      Type C_minus_1 = intensity_data(j) - Type(1.0);
      ll_theta += dgamma(C_minus_1, kappa_C_vec(j), theta_C_vec(j), true);
    }

    // Spatial prior using Cholesky solve
    Eigen::Matrix<Type, Eigen::Dynamic, 1> S_eig = S;
    Eigen::Matrix<Type, Eigen::Dynamic, 1> y_eig =
      L.template triangularView<Eigen::Lower>().solve(S_eig);

    vector<Type> y = y_eig;
    Type quad_form = (y * y).sum();

    Type log_prior_theta = -Type(0.5) * (Type(n_loc) * log(sigma2) + log_det_R +
      quad_form / sigma2);

    Type log_numerator = ll_theta + log_prior_theta;

    // =========================================================================
    // STORE RESULT
    // =========================================================================

    if(compute_denominator_only) {
      log_f_vals(s) = log_numerator;
    } else {
      log_f_vals(s) = log_numerator - log_denominator_vals(s);
    }
  }


  // =========================================================================
  // IF COMPUTING DENOMINATOR ONLY, REPORT AND EXIT
  // =========================================================================

  if(compute_denominator_only) {
    REPORT(log_f_vals);
    return Type(0);
  }

  // =============================================================================
  // MC LOG-LIKELIHOOD WITH IMPORTANCE SAMPLING
  // =============================================================================

  Type max_log_f = log_f_vals.maxCoeff();
  vector<Type> exp_vals = exp(log_f_vals - max_log_f);
  Type mean_exp = exp_vals.mean();
  Type mc_log_lik = log(mean_exp) + max_log_f;

  // =============================================================================
  // FLEXIBLE PENALTIES
  // =============================================================================

  Type penalty = Type(0.0);

  // Alpha penalty
  if(use_alpha_penalty) {
    if(alpha_penalty_type == 1) {
      // Beta(shape1, shape2) prior
      penalty -= (alpha_param1 - Type(1.0)) * log(alpha_W);
      penalty -= (alpha_param2 - Type(1.0)) * log(Type(1.0) - alpha_W);
    }
  }

  // Gamma penalty
  if(use_gamma_penalty) {
    if(gamma_penalty_type == 1) {
      // Gamma(shape, rate) prior
      Type shape = gamma_param1;
      Type rate = gamma_param2;
      penalty -= (shape - Type(1.0)) * log(gamma_W);
      penalty += rate * gamma_W;
    } else if(gamma_penalty_type == 2) {
      // Normal(mean, sd) prior
      Type mean = gamma_param1;
      Type sd = gamma_param2;
      Type diff = gamma_W - mean;
      penalty += Type(0.5) * diff * diff / (sd * sd);
    }
  }

  Type nll = -(mc_log_lik - penalty);

  // Report parameters on natural scale
  ADREPORT(alpha_W);
  ADREPORT(gamma_W);
  ADREPORT(k);
  ADREPORT(rho);
  ADREPORT(sigma2);
  ADREPORT(phi);

  return nll;
}
