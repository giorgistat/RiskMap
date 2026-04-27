#include <TMB.hpp>

// =============================================================================
// dsgm_mdiag.cpp
//
// TMB template for the doubly stochastic geostatistical model with
// multiple diagnostics (lf_mdiag).
//
// Two Binomial outcomes per observation sharing a latent NB worm burden:
//
//   Parasitological (is_mf == 1):
//     P(Y > 0) = 1 - [omega / (omega + mu*(1-exp(-alpha)))]^omega
//
//   Serological (is_mf == 0):
//     P(Y = 1) = gamma_sens * {1 - [omega / (omega + mu)]^omega}
//
// Optional MDA effect (use_mda == 1):
//   mu_W(i) = mu_W_star(i) * phi(i)
//   phi(i)  = prod_{m: t_i > u_m} [1 - alpha_W * cov_m * exp(-(t_i-u_m)/gamma_W)]
//
// Penalties (active when use_mda == 1):
//   alpha_W : Beta(a,b) or Normal on logit scale
//   gamma_W : Gamma(shape,rate) or Normal
//
// Inference: MCML with importance sampling.
//   Step 1: compute_denominator_only = 1  -> report log_f_vals at theta_0
//   Step 2: compute_denominator_only = 0  -> MCML objective
// =============================================================================

template<class Type>
Type get_distance(int i, int j, const vector<Type>& dist_vec, int n_loc) {
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

  DATA_VECTOR(y);           // Binomial counts (cases), length n
  DATA_VECTOR(units_m);     // Binomial denominators, length n
  DATA_IVECTOR(is_mf);      // 1 = Parasitological (MF), 0 = Serological, length n
  DATA_MATRIX(D);           // Covariate matrix, n x p
  DATA_VECTOR(cov_offset);  // Linear predictor offset, length n

  DATA_MATRIX(S_samples);   // MCMC samples of S, n_samples x n_loc
  DATA_IVECTOR(ID_coords);  // 0-indexed location ID for each obs, length n

  DATA_VECTOR(dist_vec);    // Compressed upper-triangle distance vector
  DATA_INTEGER(n_loc);      // Number of unique spatial locations

  DATA_SCALAR(gamma_sens);  // Fixed sensitivity of serological test

  DATA_INTEGER(fix_omega);        // 0 = estimate omega, 1 = fix to fixed_omega_val
  DATA_SCALAR(fixed_omega_val);

  // --- MDA block (only used when use_mda == 1) --------------------------------
  DATA_INTEGER(use_mda);
  DATA_VECTOR(survey_times);      // survey time per obs, length n
  DATA_VECTOR(mda_times);         // MDA round times, length n_mda
  DATA_IVECTOR(mda_i);            // row indices (0-indexed) of non-zero int_mat
  DATA_IVECTOR(mda_j);            // col indices (0-indexed) of non-zero int_mat
  DATA_VECTOR(mda_coverage);      // coverage values for non-zero entries
  DATA_INTEGER(n_mda_pairs);      // number of non-zero entries

  // Penalties on MDA parameters
  DATA_INTEGER(use_alpha_W_penalty);
  DATA_INTEGER(alpha_W_penalty_type);  // 1 = Beta on alpha_W; 2 = Normal on logit_alpha_W
  DATA_SCALAR(alpha_W_param1);
  DATA_SCALAR(alpha_W_param2);

  DATA_INTEGER(use_gamma_W_penalty);
  DATA_INTEGER(gamma_W_penalty_type);  // 1 = Gamma on gamma_W; 2 = Normal on gamma_W
  DATA_SCALAR(gamma_W_param1);
  DATA_SCALAR(gamma_W_param2);

  // Importance sampling
  DATA_INTEGER(compute_denominator_only);
  DATA_VECTOR(log_denominator_vals);  // length n_samples

  // ===========================================================================
  // PARAMETERS
  // ===========================================================================

  PARAMETER_VECTOR(beta);      // fixed-effect coefficients
  PARAMETER(log_sigma2);
  PARAMETER(log_phi);
  PARAMETER(log_nu2);          // log(tau2 / sigma2)
  PARAMETER(log_omega);        // ignored when fix_omega == 1
  PARAMETER(log_alpha);        // log MF detection rate per worm

  // MDA parameters (declared always; inactive when use_mda == 0)
  PARAMETER(logit_alpha_W);    // logit immediate worm burden drop in (0,1)
  PARAMETER(log_gamma_W);      // log recovery rate

  // Natural-scale transforms
  const Type sigma2  = exp(log_sigma2);
  const Type phi     = exp(log_phi);
  const Type nu2     = exp(log_nu2);
  const Type omega   = (fix_omega == 1) ? fixed_omega_val : exp(log_omega);
  const Type alpha   = exp(log_alpha);
  const Type alpha_W = Type(1.0) / (Type(1.0) + exp(-logit_alpha_W));
  const Type gamma_W = exp(log_gamma_W);

  int n         = y.size();
  int n_samples = S_samples.rows();

  // ===========================================================================
  // MDA EFFECT  phi(i)
  // ===========================================================================

  vector<Type> mda_effect(n);
  mda_effect.setOnes();

  if (use_mda == 1) {
    for (int idx = 0; idx < n_mda_pairs; idx++) {
      int  i        = mda_i(idx);
      int  m        = mda_j(idx);
      Type coverage = mda_coverage(idx);
      Type dt       = survey_times(i) - mda_times(m);
      if (dt > Type(0))
        mda_effect(i) *= (Type(1.0) - alpha_W * coverage * exp(-dt / gamma_W));
    }
  }

  // ===========================================================================
  // CORRELATION MATRIX R AND CHOLESKY
  // ===========================================================================

  matrix<Type> R(n_loc, n_loc);
  for (int i = 0; i < n_loc; i++) {
    R(i, i) = Type(1.0) + nu2;
    for (int j = i + 1; j < n_loc; j++) {
      Type r  = exp(-get_distance(i, j, dist_vec, n_loc) / phi);
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
  const Type   c_alpha  = Type(1.0) - exp(-alpha);

  // ===========================================================================
  // MONTE CARLO LOOP
  // ===========================================================================

  vector<Type> log_f_vals(n_samples);

  for (int s = 0; s < n_samples; s++) {

    vector<Type> S_s = S_samples.row(s);

    // Worm burden: exp(fixed + spatial) * MDA decay
    vector<Type> mu_W(n);
    for (int i = 0; i < n; i++)
      mu_W(i) = exp(mu_fixed(i) + S_s(ID_coords(i))) * mda_effect(i);

    // Binomial log-likelihood
    Type ll = Type(0.0);
    for (int i = 0; i < n; i++) {
      Type p;
      if (is_mf(i) == 1) {
        // Parasitological
        p = Type(1.0) - pow(omega / (omega + mu_W(i) * c_alpha), omega);
      } else {
        // Serological
        p = gamma_sens * (Type(1.0) - pow(omega / (omega + mu_W(i)), omega));
      }
      p = CppAD::CondExpLt(p, Type(1e-10), Type(1e-10), p);
      p = CppAD::CondExpGt(p, Type(1.0) - Type(1e-10), Type(1.0) - Type(1e-10), p);
      ll += y(i) * log(p) + (units_m(i) - y(i)) * log(Type(1.0) - p);
    }

    // GP prior via Cholesky solve: S^T R^{-1} S = z^T z where L z = S
    Eigen::Matrix<Type, Eigen::Dynamic, 1> S_eig = S_s;
    Eigen::Matrix<Type, Eigen::Dynamic, 1> z =
      L.template triangularView<Eigen::Lower>().solve(S_eig);
    vector<Type> z_v = z;
    Type quad = (z_v * z_v).sum();

    Type log_prior = -Type(0.5) * (Type(n_loc) * log(sigma2) + log_det_R +
                                   quad / sigma2);
    Type log_num   = ll + log_prior;

    log_f_vals(s) = compute_denominator_only ?
                    log_num : log_num - log_denominator_vals(s);
  }

  // ===========================================================================
  // DENOMINATOR PASS: report and exit
  // ===========================================================================

  if (compute_denominator_only) {
    REPORT(log_f_vals);
    return Type(0);
  }

  // ===========================================================================
  // MC LOG-LIKELIHOOD (log-sum-exp stable)
  // ===========================================================================

  Type max_lf     = log_f_vals.maxCoeff();
  vector<Type> ef  = exp(log_f_vals - max_lf);
  Type mc_loglik   = log(ef.mean()) + max_lf;

  // ===========================================================================
  // MDA PENALTIES
  // ===========================================================================

  Type penalty = Type(0.0);

  if (use_mda == 1) {

    if (use_alpha_W_penalty == 1) {
      if (alpha_W_penalty_type == 1) {
        // Beta(a, b) on alpha_W
        penalty -= (alpha_W_param1 - Type(1.0)) * log(alpha_W);
        penalty -= (alpha_W_param2 - Type(1.0)) * log(Type(1.0) - alpha_W);
      } else if (alpha_W_penalty_type == 2) {
        // Normal(mean, sd) on logit(alpha_W)
        Type d = logit_alpha_W - alpha_W_param1;
        penalty += Type(0.5) * d * d / (alpha_W_param2 * alpha_W_param2);
      }
    }

    if (use_gamma_W_penalty == 1) {
      if (gamma_W_penalty_type == 1) {
        // Gamma(shape, rate) on gamma_W
        penalty -= (gamma_W_param1 - Type(1.0)) * log(gamma_W);
        penalty += gamma_W_param2 * gamma_W;
      } else if (gamma_W_penalty_type == 2) {
        Type d = gamma_W - gamma_W_param1;
        penalty += Type(0.5) * d * d / (gamma_W_param2 * gamma_W_param2);
      } else if (gamma_W_penalty_type == 3) {
        Type d = log_gamma_W - gamma_W_param1;
        penalty += Type(0.5) * d * d / (gamma_W_param2 * gamma_W_param2);
      }
    }
  }

  Type nll = -(mc_loglik - penalty);

  // ===========================================================================
  // REPORT ON NATURAL SCALE
  // ===========================================================================

  ADREPORT(sigma2);
  ADREPORT(phi);
  ADREPORT(omega);
  ADREPORT(alpha);
  Type tau2 = nu2 * sigma2;
  ADREPORT(tau2);

  if (use_mda == 1) {
    ADREPORT(alpha_W);
    ADREPORT(gamma_W);
  }

  return nll;
}
