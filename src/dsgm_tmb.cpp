#include <TMB.hpp>

// =============================================================================
// dsgm_tmb.cpp
//
// TMB template for the doubly stochastic geostatistical model for joint
// prevalence-intensity modelling of soil-transmitted helminths (STH).
//
// Hierarchical model:
//   W_j(x_i) | S(x_i) ~ NegBin(mu(x_i), k)
//     log mu*(x_i) = D_i * beta + S(x_i)
//     mu(x_i)      = mu*(x_i) * MDA_effect(x_i, t)
//
//   Y_ij | W_j(x_i) ~ Poisson(rho * W_j(x_i))
//
// Intensity likelihood for C | C > 0 (controlled by intensity_family):
//   0 = Shifted Gamma:           C - 1 ~ Gamma(kappa_C, theta_C)
//   1 = Zero-truncated NegBin2:  C | C > 0 ~ NegBin2(mu_C, phi_C) / (1 - p0)
//
// MDA effect (exponential decay):
//   MDA_effect(i) = prod_{m: t_i > u_m} [1 - alpha_W * cov_m * exp(-(t_i-u_m)/gamma_W)]
//
// Inference: MCML with importance sampling.
//   Step 1: compute_denominator_only = 1  -> report log_f_vals at theta_0
//   Step 2: compute_denominator_only = 0  -> MCML objective
//
// Penalty types for gamma_W:
//   1 = Gamma(shape, rate) on gamma_W          gradient ~ O(gamma_W)
//   2 = Normal(mean, sd)   on gamma_W          gradient ~ O(gamma_W - mean)
//   3 = Log-Normal: Normal(mu, sd) on log_gamma gradient ~ O(1)  [preferred]
// =============================================================================

template<class Type>
Type get_distance_sth(int i, int j, const vector<Type>& dist_vec, int n_loc) {
  if (i == j) return Type(0);
  if (i > j)  { int tmp = i; i = j; j = tmp; }
  int idx = i * n_loc - (i * (i + 1)) / 2 + (j - i - 1);
  return dist_vec(idx);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ===========================================================================
  // DATA
  // ===========================================================================

  DATA_VECTOR(y_prev);          // Binary infection status (0/1), length n
  DATA_VECTOR(intensity_data);  // EPG for egg-positive individuals, length n_pos
  DATA_IVECTOR(pos_idx);        // 0-indexed positions of positives in y_prev

  DATA_MATRIX(D);               // Covariate matrix, n x p
  DATA_VECTOR(cov_offset);      // Linear predictor offset, length n

  DATA_MATRIX(S_samples);       // Stan MCMC samples of S, n_samples x n_loc
  DATA_IVECTOR(ID_coords);      // 0-indexed location ID per obs, length n

  DATA_VECTOR(dist_vec);        // Compressed upper-triangle distance vector
  DATA_INTEGER(n_loc);          // Number of unique spatial locations

  // MDA block
  DATA_VECTOR(survey_times);    // Survey time per obs, length n
  DATA_VECTOR(mda_times);       // MDA round times, length n_mda
  DATA_IVECTOR(mda_i);          // Row indices (0-indexed) of non-zero int_mat entries
  DATA_IVECTOR(mda_j);          // Col indices (0-indexed) of non-zero int_mat entries
  DATA_VECTOR(mda_coverage);    // Coverage values for non-zero entries
  DATA_INTEGER(n_mda_pairs);    // Number of non-zero entries

  // Alpha penalty
  DATA_INTEGER(use_alpha_penalty);
  DATA_INTEGER(alpha_penalty_type);  // 1 = Beta(a,b) on alpha_W; 2 = Normal on logit
  DATA_SCALAR(alpha_param1);
  DATA_SCALAR(alpha_param2);

  // Gamma_W penalty
  DATA_INTEGER(use_gamma_penalty);
  DATA_INTEGER(gamma_penalty_type);  // 1 = Gamma(shape,rate); 2 = Normal; 3 = Log-Normal
  DATA_SCALAR(gamma_param1);
  DATA_SCALAR(gamma_param2);

  // Rho penalty
  DATA_INTEGER(use_rho_penalty);
  DATA_INTEGER(rho_penalty_type);  // 1 = Gamma(shape,rate); 2 = Normal; 3 = Log-Normal on log_rho
  DATA_SCALAR(rho_param1);
  DATA_SCALAR(rho_param2);

  // Importance sampling
  DATA_INTEGER(compute_denominator_only);
  DATA_VECTOR(log_denominator_vals);  // length n_samples

  // Intensity likelihood family: 0 = shifted Gamma, 1 = zero-truncated NegBin
  DATA_INTEGER(intensity_family);

  // ===========================================================================
  // PARAMETERS
  // ===========================================================================

  PARAMETER_VECTOR(beta);       // Fixed-effect coefficients
  PARAMETER(log_k);             // log NB aggregation parameter
  PARAMETER(log_rho);           // log per-worm egg detection rate
  PARAMETER(logit_alpha);       // logit immediate worm burden reduction in (0,1)
  PARAMETER(log_gamma);         // log worm burden recovery rate  <- raw TMB parameter
  PARAMETER(log_sigma2);        // log GP variance
  PARAMETER(log_phi);           // log GP range

  // Natural-scale transforms
  const Type k       = exp(log_k);
  const Type rho     = exp(log_rho);
  const Type alpha_W = Type(1.0) / (Type(1.0) + exp(-logit_alpha));
  const Type gamma_W = exp(log_gamma);
  const Type sigma2  = exp(log_sigma2);
  const Type phi     = exp(log_phi);

  int n         = y_prev.size();
  int n_pos     = intensity_data.size();
  int n_samples = S_samples.rows();

  // ===========================================================================
  // MDA EFFECT
  // ===========================================================================

  vector<Type> mda_effect(n);
  mda_effect.setOnes();

  for (int idx = 0; idx < n_mda_pairs; idx++) {
    int  i        = mda_i(idx);
    int  m        = mda_j(idx);
    Type coverage = mda_coverage(idx);
    Type dt       = survey_times(i) - mda_times(m);
    if (dt > Type(0))
      mda_effect(i) *= (Type(1.0) - alpha_W * coverage * exp(-dt / gamma_W));
  }

  // ===========================================================================
  // CORRELATION MATRIX AND CHOLESKY (outside MC loop)
  // ===========================================================================

  matrix<Type> R(n_loc, n_loc);
  for (int i = 0; i < n_loc; i++) {
    R(i, i) = Type(1.0);
    for (int j = i + 1; j < n_loc; j++) {
      Type r  = exp(-get_distance_sth(i, j, dist_vec, n_loc) / phi);
      R(i, j) = r;
      R(j, i) = r;
    }
  }

  Eigen::LLT<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> llt(R);
  matrix<Type> L = llt.matrixL();

  Type log_det_R = Type(0.0);
  for (int i = 0; i < n_loc; i++) log_det_R += log(L(i, i));
  log_det_R *= Type(2.0);

  // ===========================================================================
  // FIXED-EFFECTS LINEAR PREDICTOR (outside MC loop)
  // ===========================================================================

  vector<Type> mu_fixed = D * beta + cov_offset;
  const Type   c_rho    = Type(1.0) - exp(-rho);

  // ===========================================================================
  // MONTE CARLO LOOP
  // ===========================================================================

  vector<Type> log_f_vals(n_samples);

  for (int s = 0; s < n_samples; s++) {

    vector<Type> S_s = S_samples.row(s);

    // Mean worm burden with MDA
    vector<Type> mu_W(n);
    for (int i = 0; i < n; i++)
      mu_W(i) = exp(mu_fixed(i) + S_s(ID_coords(i))) * mda_effect(i);

    // ----- Prevalence likelihood -----
    Type ll = Type(0.0);

    for (int i = 0; i < n; i++) {
      Type ratio = k / (k + mu_W(i) * c_rho);
      Type pr    = Type(1.0) - pow(ratio, k);
      pr = CppAD::CondExpLt(pr, Type(1e-10),               Type(1e-10),               pr);
      pr = CppAD::CondExpGt(pr, Type(1.0) - Type(1e-10),   Type(1.0) - Type(1e-10),   pr);

      if (y_prev(i) > Type(0.5))
        ll += log(pr);
      else
        ll += log(Type(1.0) - pr);
    }

    // ----- Intensity likelihood (positives only) -----
    for (int idx = 0; idx < n_pos; idx++) {
      int  i     = pos_idx(idx);
      Type mu_i  = mu_W(i);
      Type ratio = k / (k + mu_i * c_rho);
      Type pr_i  = Type(1.0) - pow(ratio, k);
      pr_i = CppAD::CondExpLt(pr_i, Type(1e-10), Type(1e-10), pr_i);

      // Conditional moments of C | C > 0
      Type mu_C     = (rho * mu_i) / pr_i;
      Type sigma2_C = (rho * mu_i * (Type(1.0) + rho)) / pr_i
      + (rho * rho * mu_i * mu_i / pr_i)
        * (Type(1.0)/k + Type(1.0) - Type(1.0)/pr_i);
        sigma2_C = CppAD::CondExpLt(sigma2_C, Type(1e-6), Type(1e-6), sigma2_C);

        if (intensity_family == 0) {

          // ------------------------------------------------------------------
          // Shifted Gamma: C - 1 ~ Gamma(kappa_C, theta_C)
          // mu_C1 floored at 0.1 to prevent shape/rate collapsing to zero
          // when mu_C -> 1 at low worm burden
          // ------------------------------------------------------------------
          Type mu_C1 = mu_C - Type(1.0);
          mu_C1 = CppAD::CondExpLt(mu_C1, Type(0.1), Type(0.1), mu_C1);

          Type kappa_C = (mu_C1 * mu_C1) / sigma2_C;
          Type theta_C = sigma2_C / mu_C1;
          kappa_C = CppAD::CondExpLt(kappa_C, Type(1e-6), Type(1e-6), kappa_C);
          theta_C = CppAD::CondExpLt(theta_C, Type(1e-6), Type(1e-6), theta_C);

          Type y_shift = intensity_data(idx) - Type(1.0);
          y_shift = CppAD::CondExpLt(y_shift, Type(0.0), Type(0.0), y_shift);

          ll += (kappa_C - Type(1.0)) * log(y_shift + Type(1e-10))
            - y_shift / theta_C
          - kappa_C * log(theta_C)
          - lgamma(kappa_C);

        } else {

          // ------------------------------------------------------------------
          // Zero-truncated NegBin2: C | C > 0 ~ NegBin2(mu_C, phi_C) / (1 - p0)
          // phi_C = mu_C^2 / (sigma2_C - mu_C)
          // p0    = P(C = 0) under the untruncated NegBin2
          // log NegBin2(c | mu, phi) = lgamma(c+phi) - lgamma(phi) - lgamma(c+1)
          //                          + phi*log(phi/(phi+mu)) + c*log(mu/(phi+mu))
          // ------------------------------------------------------------------
          Type denom_nb = sigma2_C - mu_C;
          denom_nb = CppAD::CondExpLt(denom_nb, Type(1e-6), Type(1e-6), denom_nb);
          Type phi_C = mu_C * mu_C / denom_nb;
          phi_C = CppAD::CondExpLt(phi_C, Type(1e-4), Type(1e-4), phi_C);

          Type c        = intensity_data(idx);
          Type log_r    = log(phi_C) - log(phi_C + mu_C);  // log(phi/(phi+mu))
          Type log_1mr  = log(mu_C)  - log(phi_C + mu_C);  // log(mu/(phi+mu))

          // log NegBin2(c | mu_C, phi_C)
          Type log_nb = lgamma(c + phi_C) - lgamma(phi_C) - lgamma(c + Type(1.0))
            + phi_C * log_r
          + c     * log_1mr;

          // log p0 = NegBin2(0 | mu_C, phi_C) = phi_C * log(phi_C/(phi_C+mu_C))
          Type log_p0   = phi_C * log_r;

          // log(1 - p0): use log1m_exp for numerical stability when p0 is close to 1
          Type log1m_p0 = log(Type(1.0) - exp(log_p0));
          log1m_p0 = CppAD::CondExpLt(log1m_p0, Type(-30.0), Type(-30.0), log1m_p0);

          ll += log_nb - log1m_p0;
        }
    }

    // ----- GP prior -----
    Eigen::Matrix<Type, Eigen::Dynamic, 1> S_eig = S_s;
    Eigen::Matrix<Type, Eigen::Dynamic, 1> z =
      L.template triangularView<Eigen::Lower>().solve(S_eig);
    vector<Type> z_v = z;
    Type quad = (z_v * z_v).sum();

    Type log_prior = -Type(0.5) * (Type(n_loc) * log(sigma2) + log_det_R +
      quad / sigma2);

    Type log_num  = ll + log_prior;
    log_f_vals(s) = compute_denominator_only ?
    log_num : log_num - log_denominator_vals(s);
  }

  // ===========================================================================
  // DENOMINATOR PASS
  // ===========================================================================

  if (compute_denominator_only) {
    REPORT(log_f_vals);
    return Type(0);
  }

  // ===========================================================================
  // MC LOG-LIKELIHOOD (log-sum-exp stable)
  // ===========================================================================

  Type max_lf     = log_f_vals.maxCoeff();
  vector<Type> ef = exp(log_f_vals - max_lf);
  Type mc_loglik  = log(ef.mean()) + max_lf;

  // ===========================================================================
  // PENALTIES
  // ===========================================================================

  Type penalty = Type(0.0);

  // --- Alpha_W ---
  if (use_alpha_penalty == 1) {
    if (alpha_penalty_type == 1) {
      // Beta(a, b) on alpha_W
      penalty -= (alpha_param1 - Type(1.0)) * log(alpha_W);
      penalty -= (alpha_param2 - Type(1.0)) * log(Type(1.0) - alpha_W);
    } else if (alpha_penalty_type == 2) {
      // Normal(mean, sd) on logit(alpha_W)
      Type d = logit_alpha - alpha_param1;
      penalty += Type(0.5) * d * d / (alpha_param2 * alpha_param2);
    }
  }

  // --- Gamma_W ---
  if (use_gamma_penalty == 1) {
    if (gamma_penalty_type == 1) {
      // Gamma(shape, rate) on gamma_W — gradient O(gamma_W), can cause IS collapse
      penalty -= (gamma_param1 - Type(1.0)) * log(gamma_W);
      penalty += gamma_param2 * gamma_W;
    } else if (gamma_penalty_type == 2) {
      // Normal(mean, sd) on gamma_W directly
      Type d = gamma_W - gamma_param1;
      penalty += Type(0.5) * d * d / (gamma_param2 * gamma_param2);
    } else if (gamma_penalty_type == 3) {
      // Log-Normal: Normal(mu, sd) on log(gamma_W)
      // log_gamma is the raw TMB parameter so gradient is O(1) — preferred
      Type d = log_gamma - gamma_param1;
      penalty += Type(0.5) * d * d / (gamma_param2 * gamma_param2);
    }
  }

  // --- Rho ---
  if (use_rho_penalty == 1) {
    if (rho_penalty_type == 1) {
      // Gamma(shape, rate) on rho
      penalty -= (rho_param1 - Type(1.0)) * log(rho);
      penalty += rho_param2 * rho;
    } else if (rho_penalty_type == 2) {
      // Normal(mean, sd) on rho
      Type d = rho - rho_param1;
      penalty += Type(0.5) * d * d / (rho_param2 * rho_param2);
    } else if (rho_penalty_type == 3) {
      // Log-Normal: Normal(mu, sd) on log_rho — preferred, gradient O(1)
      Type d = log_rho - rho_param1;
      penalty += Type(0.5) * d * d / (rho_param2 * rho_param2);
    }
  }

  Type nll = -(mc_loglik - penalty);

  // ===========================================================================
  // REPORT ON NATURAL SCALE
  // ===========================================================================

  ADREPORT(k);
  ADREPORT(rho);
  ADREPORT(alpha_W);
  ADREPORT(gamma_W);
  ADREPORT(sigma2);
  ADREPORT(phi);

  return nll;
}
