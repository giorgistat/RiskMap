#include <TMB.hpp>

// =============================================================================
// dsgm_tmb.cpp
//
// TMB template for the doubly stochastic geostatistical model for joint
// prevalence-intensity modelling of soil-transmitted helminths (STH).
//
// Hierarchical model:
//   W_j(x_i) | S(x_i) ~ NegBin(mu(x_i), omega_i)
//     log mu*(x_i) = D_i * beta + S(x_i)
//     mu(x_i)      = mu*(x_i) * MDA_effect(x_i, t)
//
//   Aggregation parameter (vary_k controls which form is used):
//     vary_k = 0:  omega_i = exp(log_k)                      [constant, default]
//     vary_k = 1:  log(omega_i) = log_k + omega1*log(mu_W_i) [varies with burden]
//
//   omega1 is a free, unpenalised regression slope — well-identified from data.
//
//   Y_ij | W_j(x_i) ~ Poisson(rho * W_j(x_i))
//
// Intensity likelihood for C | C > 0 (controlled by intensity_family):
//   0 = Shifted Gamma:           C - 1 ~ Gamma(kappa_C, theta_C)
//   1 = Shifted NegBin2:         C - 1 ~ NB2(mu_C - 1, phi_C)
//
// MDA effect (exponential decay):
//   MDA_effect(i) = prod_{m: t_i > u_m} [1 - alpha_W * cov_m * exp(-(t_i-u_m)/gamma_W)]
//
// Inference: MCML with importance sampling.
//   Step 1: compute_denominator_only = 1  -> report log_f_vals at theta_0
//   Step 2: compute_denominator_only = 0  -> MCML objective
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

  DATA_VECTOR(y_prev);
  DATA_VECTOR(intensity_data);
  DATA_IVECTOR(pos_idx);

  DATA_MATRIX(D);
  DATA_VECTOR(cov_offset);

  DATA_MATRIX(S_samples);
  DATA_IVECTOR(ID_coords);

  DATA_VECTOR(dist_vec);
  DATA_INTEGER(n_loc);

  // MDA block
  DATA_VECTOR(survey_times);
  DATA_VECTOR(mda_times);
  DATA_IVECTOR(mda_i);
  DATA_IVECTOR(mda_j);
  DATA_VECTOR(mda_coverage);
  DATA_INTEGER(n_mda_pairs);

  // Alpha penalty
  DATA_INTEGER(use_alpha_penalty);
  DATA_INTEGER(alpha_penalty_type);  // 1 = Beta(a,b); 2 = Normal on logit
  DATA_SCALAR(alpha_param1);
  DATA_SCALAR(alpha_param2);

  // Gamma_W penalty
  DATA_INTEGER(use_gamma_penalty);
  DATA_INTEGER(gamma_penalty_type);  // 1 = Gamma; 2 = Normal; 3 = Log-Normal [preferred]
  DATA_SCALAR(gamma_param1);
  DATA_SCALAR(gamma_param2);

  // Rho penalty
  DATA_INTEGER(use_rho_penalty);
  DATA_INTEGER(rho_penalty_type);    // 1 = Gamma; 2 = Normal; 3 = Log-Normal [preferred]
  DATA_SCALAR(rho_param1);
  DATA_SCALAR(rho_param2);

  // Importance sampling
  DATA_INTEGER(compute_denominator_only);
  DATA_VECTOR(log_denominator_vals);

  // Intensity likelihood family: 0 = shifted Gamma, 1 = shifted NegBin2
  DATA_INTEGER(intensity_family);

  // Aggregation flag: 0 = constant omega, 1 = omega varies with mu_W
  DATA_INTEGER(vary_k);

  // ===========================================================================
  // PARAMETERS
  // ===========================================================================

  PARAMETER_VECTOR(beta);
  PARAMETER(log_k);       // log baseline aggregation (omega_0)
  PARAMETER(log_rho);
  PARAMETER(logit_alpha);
  PARAMETER(log_gamma);
  PARAMETER(log_sigma2);
  PARAMETER(log_phi);
  PARAMETER(omega1);      // slope of log(omega) on log(mu_W); fixed at 0 when vary_k=0

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
  // CORRELATION MATRIX AND CHOLESKY
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
  // FIXED-EFFECTS LINEAR PREDICTOR
  // ===========================================================================

  vector<Type> mu_fixed = D * beta + cov_offset;
  const Type   c_rho    = Type(1.0) - exp(-rho);

  // ===========================================================================
  // MONTE CARLO LOOP
  // ===========================================================================

  vector<Type> log_f_vals(n_samples);

  for (int s = 0; s < n_samples; s++) {

    vector<Type> S_s = S_samples.row(s);

    vector<Type> mu_W(n);
    for (int i = 0; i < n; i++)
      mu_W(i) = exp(mu_fixed(i) + S_s(ID_coords(i))) * mda_effect(i);

    // ----- Prevalence likelihood -----
    Type ll = Type(0.0);

    for (int i = 0; i < n; i++) {
      Type k_i = (vary_k == 1) ?
      exp(log_k + omega1 * log(mu_W(i) + Type(1e-10))) : k;

      Type ratio = k_i / (k_i + mu_W(i) * c_rho);
      Type pr    = Type(1.0) - pow(ratio, k_i);
      pr = CppAD::CondExpLt(pr, Type(1e-10),             Type(1e-10),             pr);
      pr = CppAD::CondExpGt(pr, Type(1.0) - Type(1e-10), Type(1.0) - Type(1e-10), pr);

      if (y_prev(i) > Type(0.5))
        ll += log(pr);
      else
        ll += log(Type(1.0) - pr);
    }

    // ----- Intensity likelihood (positives only) -----
    for (int idx = 0; idx < n_pos; idx++) {
      int  i    = pos_idx(idx);
      Type mu_i = mu_W(i);
      Type k_i  = (vary_k == 1) ?
      exp(log_k + omega1 * log(mu_i + Type(1e-10))) : k;

      Type ratio = k_i / (k_i + mu_i * c_rho);
      Type pr_i  = Type(1.0) - pow(ratio, k_i);
      pr_i = CppAD::CondExpLt(pr_i, Type(1e-10), Type(1e-10), pr_i);

      Type mu_C     = (rho * mu_i) / pr_i;
      Type sigma2_C = (rho * mu_i * (Type(1.0) + rho)) / pr_i
      + (rho * rho * mu_i * mu_i / pr_i)
        * (Type(1.0)/k_i + Type(1.0) - Type(1.0)/pr_i);
        sigma2_C = CppAD::CondExpLt(sigma2_C, Type(1e-6), Type(1e-6), sigma2_C);

        if (intensity_family == 0) {

          // Shifted Gamma: C - 1 ~ Gamma(kappa_C, theta_C)
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

          // Shifted NegBin2: C - 1 ~ NB2(mu_C1, phi_C)
          // Moments of C - 1: mean = mu_C - 1, variance = sigma2_C
          Type mu_C1 = mu_C - Type(1.0);
          mu_C1 = CppAD::CondExpLt(mu_C1, Type(0.1), Type(0.1), mu_C1);

          Type denom_nb = sigma2_C - mu_C1;
          denom_nb = CppAD::CondExpLt(denom_nb, Type(1e-6), Type(1e-6), denom_nb);
          Type phi_C = mu_C1 * mu_C1 / denom_nb;
          phi_C = CppAD::CondExpLt(phi_C, Type(1e-4), Type(1e-4), phi_C);

          Type c_shift = intensity_data(idx) - Type(1.0);   // C - 1 >= 0
          Type log_r   = log(phi_C)  - log(phi_C + mu_C1);
          Type log_1mr = log(mu_C1)  - log(phi_C + mu_C1);

          ll += lgamma(c_shift + phi_C) - lgamma(phi_C) - lgamma(c_shift + Type(1.0))
            + phi_C * log_r + c_shift * log_1mr;
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
  // PENALTIES  (alpha_W, gamma_W, rho only — omega1 is free)
  // ===========================================================================

  Type penalty = Type(0.0);

  if (use_alpha_penalty == 1) {
    if (alpha_penalty_type == 1) {
      penalty -= (alpha_param1 - Type(1.0)) * log(alpha_W);
      penalty -= (alpha_param2 - Type(1.0)) * log(Type(1.0) - alpha_W);
    } else if (alpha_penalty_type == 2) {
      Type d = logit_alpha - alpha_param1;
      penalty += Type(0.5) * d * d / (alpha_param2 * alpha_param2);
    }
  }

  if (use_gamma_penalty == 1) {
    if (gamma_penalty_type == 1) {
      penalty -= (gamma_param1 - Type(1.0)) * log(gamma_W);
      penalty += gamma_param2 * gamma_W;
    } else if (gamma_penalty_type == 2) {
      Type d = gamma_W - gamma_param1;
      penalty += Type(0.5) * d * d / (gamma_param2 * gamma_param2);
    } else if (gamma_penalty_type == 3) {
      Type d = log_gamma - gamma_param1;
      penalty += Type(0.5) * d * d / (gamma_param2 * gamma_param2);
    }
  }

  if (use_rho_penalty == 1) {
    if (rho_penalty_type == 1) {
      penalty -= (rho_param1 - Type(1.0)) * log(rho);
      penalty += rho_param2 * rho;
    } else if (rho_penalty_type == 2) {
      Type d = rho - rho_param1;
      penalty += Type(0.5) * d * d / (rho_param2 * rho_param2);
    } else if (rho_penalty_type == 3) {
      Type d = log_rho - rho_param1;
      penalty += Type(0.5) * d * d / (rho_param2 * rho_param2);
    }
  }

  Type nll = -(mc_loglik - penalty);

  ADREPORT(k);
  ADREPORT(rho);
  ADREPORT(alpha_W);
  ADREPORT(gamma_W);
  ADREPORT(sigma2);
  ADREPORT(phi);
  ADREPORT(omega1);

  return nll;
}
