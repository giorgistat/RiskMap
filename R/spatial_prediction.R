##' @importFrom stats setNames
##' @importFrom utils head

##' @title Prediction of the random effects components and covariates effects over a spatial grid
##' @description Computes predictions over a spatial grid using a fitted model from
##'   \code{\link{glgpm}} or \code{\link{dsgm}}.
##' @param object A RiskMap object.
##' @param grid_pred An \code{sfc} or \code{sf} of POINT geometries, or a list thereof for joint predictions.
##' @param predictors Optional data frame of predictor variables at prediction locations.
##' @param re_predictors Optional data frame of random effect predictors.
##' @param pred_cov_offset Optional numeric vector of covariate offsets at prediction locations.
##' @param control_sim Control parameters from \code{\link{set_control_sim}}.
##' @param type \code{"marginal"} or \code{"joint"}.
##' @param messages Logical; display progress messages. Default \code{TRUE}.
##' @return An object of class \code{"RiskMap.pred.re"}.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @importFrom Matrix solve
##' @export
pred_over_grid <- function(object,
                           grid_pred = NULL,
                           predictors = NULL,
                           re_predictors = NULL,
                           pred_cov_offset = NULL,
                           control_sim = set_control_sim(),
                           type = "marginal",
                           messages = TRUE) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # ---------------------------------------------------------------------------
  # FIX 1: resolve par_hat
  # dsgm objects store parameters in $model_params; glgpm objects use coef()
  # ---------------------------------------------------------------------------
  is_dsgm <- !is.null(object$family) &&
    object$family %in% c("intprev", "lf_mdiag")

  par_hat <- if (is_dsgm) object$model_params else coef(object)

  # ---------------------------------------------------------------------------
  # grid_pred validation
  # ---------------------------------------------------------------------------
  list_mode <- is.list(grid_pred) && !is.null(grid_pred) &&
    !any(class(grid_pred) %in% c("sf", "sfc"))

  if (list_mode) {
    if (type != "joint")
      stop("When 'grid_pred' is a list, 'type' must be 'joint'.")
    if (length(grid_pred) == 0L)
      stop("'grid_pred' is a list but has length 0.")
    grid_pred <- lapply(grid_pred, function(g) {
      if (inherits(g, "sf")) sf::st_geometry(g) else g
    })
    ok_geom <- vapply(grid_pred, function(g) {
      inherits(g, "sfc") && all(sf::st_geometry_type(g) == "POINT")
    }, logical(1))
    if (!all(ok_geom))
      stop("Each element of 'grid_pred' must be an 'sf' or 'sfc' object with POINT geometries.")

  } else if (!is.null(grid_pred)) {
    if (inherits(grid_pred, "sf")) grid_pred <- sf::st_geometry(grid_pred)
    if (!inherits(grid_pred, "sfc") || !all(sf::st_geometry_type(grid_pred) == "POINT"))
      stop("'grid_pred' must be an 'sf' or 'sfc' object with POINT geometries.")
  }

  obs_loc <- is.null(grid_pred)

  if (!inherits(control_sim, what = "mcmc.RiskMap", which = FALSE))
    stop("the argument passed to 'control_sim' must be an output from set_control_sim")

  if (obs_loc) {
    predictors <- as.data.frame(sf::st_drop_geometry(object$data_sf))
    grid_pred  <- sf::st_as_sfc(object$data_sf)
    list_mode  <- FALSE
  } else {
    if (list_mode) {
      grid_pred <- lapply(grid_pred, sf::st_transform, crs = object$crs)
    } else {
      grid_pred <- sf::st_transform(grid_pred, crs = object$crs)
    }
  }

  if (list_mode) {
    grp    <- lapply(grid_pred, sf::st_coordinates)
    n_pred <- vapply(grp, nrow, integer(1))
  } else {
    grp    <- sf::st_coordinates(grid_pred)
    n_pred <- nrow(grp)
  }

  object$D <- as.matrix(object$D)
  p <- ncol(object$D)

  if (!all(length(object$cov_offset == 0))) {
    if (is.null(pred_cov_offset))
      stop("The covariate offset must be specified at each prediction location.")
    if (!inherits(pred_cov_offset, "numeric"))
      stop("'pred_cov_offset' must be a numeric vector")
    if (length(pred_cov_offset) != n_pred)
      stop("The length of 'pred_cov_offset' does not match the number of prediction locations")
  } else {
    pred_cov_offset <- 0
  }

  if (!type %in% c("marginal", "joint"))
    stop("the argument 'type' must be set to 'marginal' or 'joint'")

  inter_f      <- interpret.formula(object$formula)
  inter_lt_f   <- inter_f
  inter_lt_f$pf <- update(inter_lt_f$pf, NULL ~.)

  intercept_only <- if (p == 1) prod(object$D[, 1] == 1) == 1 else FALSE

  n_re <- length(object$re)
  if (n_re > 0 && type == "marginal" && !is.null(re_predictors))
    stop("Random effect predictions require type = 'joint'")

  # ---------------------------------------------------------------------------
  # Build mu_pred (fixed-effects linear predictor at prediction locations)
  # ---------------------------------------------------------------------------
  if (!is.null(predictors)) {

    .build_D_pred <- function(pred_i, n_i, idx = NULL) {
      if (!is.data.frame(pred_i))
        stop(if (is.null(idx)) "'predictors' must be a data.frame"
             else sprintf("'predictors[[%d]]' must be a data.frame", idx))
      if (nrow(pred_i) != n_i)
        stop(if (is.null(idx)) "Rows of 'predictors' do not match 'grid_pred'"
             else sprintf("Rows of 'predictors[[%d]]' do not match 'grid_pred[[%d]]'", idx, idx))
      mf <- model.frame(inter_lt_f$pf, data = pred_i, na.action = na.fail)
      D  <- as.matrix(model.matrix(attr(mf, "terms"), data = pred_i))
      if (ncol(D) != p)
        stop(if (is.null(idx)) "Predictors do not match the model formula"
             else sprintf("Predictors in group %d do not match the model formula", idx))
      D
    }

    if (list_mode) {
      D_pred  <- mapply(.build_D_pred, predictors, n_pred, seq_along(predictors), SIMPLIFY = FALSE)
      mu_pred <- lapply(D_pred, function(D) as.numeric(D %*% par_hat$beta))
    } else {
      D_pred  <- .build_D_pred(predictors, n_pred)
      mu_pred <- as.numeric(D_pred %*% par_hat$beta)
    }

  } else if (intercept_only) {
    mu_pred <- if (list_mode) lapply(n_pred, function(n) rep(par_hat$beta, n)) else par_hat$beta
  } else {
    mu_pred <- 0
  }

  # ---------------------------------------------------------------------------
  # Random effects (unchanged from original)
  # ---------------------------------------------------------------------------
  if (n_re > 0) {
    ID_g         <- as.matrix(cbind(object$ID_coords, object$ID_re))
    re_unique    <- object$re
    n_dim_re_tot <- sapply(seq_len(n_re + 1), function(i) length(unique(ID_g[, i])))

    if (!is.null(re_predictors)) {
      if (any(is.na(re_predictors))) {
        warning("Missing values in 're_predictors'; removing affected locations")
        ind_c       <- complete.cases(re_predictors)
        re_predictors <- re_predictors[ind_c, , drop = FALSE]
        grid_pred   <- if (list_mode) lapply(grid_pred, `[`, ind_c) else grid_pred[ind_c]
        grp         <- if (list_mode) lapply(grid_pred, sf::st_coordinates) else sf::st_coordinates(grid_pred)
        n_pred      <- if (list_mode) vapply(grp, nrow, integer(1)) else nrow(grp)
      }
      if (!is.data.frame(re_predictors)) stop("'re_predictors' must be a data.frame")
      if (nrow(re_predictors) != n_pred)  stop("'re_predictors' rows do not match 'grid_pred'")
      if (ncol(re_predictors) != n_re)    stop("'re_predictors' columns do not match number of random effects")

      n_dim_re <- sapply(seq_len(n_re), function(i) length(unique(object$ID_re[, i])))
      D_re_pred <- lapply(seq_len(n_re), function(i) {
        mat <- matrix(0, nrow = n_pred, ncol = n_dim_re[i])
        for (j in seq_along(object$re[[i]]))
          mat[re_predictors[, i] == object$re[[i]][j], j] <- 1
        mat
      })
    } else {
      D_re_pred <- NULL
    }
  } else {
    D_re_pred <- NULL
  }

  # ---------------------------------------------------------------------------
  # Spatial quantities
  # ---------------------------------------------------------------------------
  out <- list(mu_pred = mu_pred, grid_pred = grid_pred, par_hat = par_hat)

  if (object$scale_to_km) {
    grp <- if (list_mode) lapply(grp, function(g) g / 1000) else grp / 1000
  }

  if (object$family != "gaussian" && !obs_loc) {
    if (list_mode) {
      U_pred <- lapply(grp, function(g) {
        t(sapply(seq_len(nrow(g)), function(i)
          sqrt((object$coords[, 1] - g[i, 1])^2 + (object$coords[, 2] - g[i, 2])^2)))
      })
    } else {
      U_pred <- t(sapply(seq_len(n_pred), function(i)
        sqrt((object$coords[, 1] - grp[i, 1])^2 + (object$coords[, 2] - grp[i, 2])^2)))
    }
  } else if (object$family == "gaussian" && !obs_loc) {
    if (list_mode) {
      U_pred <- lapply(grp, function(g) {
        t(sapply(seq_len(nrow(g)), function(i)
          sqrt((object$coords[object$ID_coords, 1] - g[i, 1])^2 +
                 (object$coords[object$ID_coords, 2] - g[i, 2])^2)))
      })
    } else {
      U_pred <- t(sapply(seq_len(n_pred), function(i)
        sqrt((object$coords[object$ID_coords, 1] - grp[i, 1])^2 +
               (object$coords[object$ID_coords, 2] - grp[i, 2])^2)))
    }
  }

  U <- dist(object$coords)
  if (!obs_loc) {
    C <- if (list_mode)
      lapply(U_pred, function(u) par_hat$sigma2 * matern_cor(u, phi = par_hat$phi, kappa = object$kappa))
    else
      par_hat$sigma2 * matern_cor(U_pred, phi = par_hat$phi, kappa = object$kappa)
  }

  mu <- as.numeric(object$D %*% par_hat$beta)

  n_samples <- if (control_sim$linear_model) {
                 control_sim$n_sim
               } else {
                 (control_sim$n_sim - control_sim$burnin) / control_sim$thin
               }

  R <- matern_cor(U, phi = par_hat$phi, kappa = object$kappa, return_sym_matrix = TRUE)

  # ---------------------------------------------------------------------------
  # FIX 2: nu2 / nugget
  # STH: tau2 fixed at 0 -> near-zero for numerical stability
  # LF:  tau2 may be estimated -> use it
  # glgpm: existing behaviour
  # ---------------------------------------------------------------------------
  if (is_dsgm) {
    tau2_val <- par_hat$tau2 %||% 0
    nu2 <- if (object$family == "intprev" || tau2_val == 0) 1e-10
    else tau2_val / par_hat$sigma2
  } else {
    nu2 <- if (!is.null(object$fix_tau2)) object$fix_tau2 / par_hat$sigma2
    else par_hat$tau2 / par_hat$sigma2
    if (nu2 == 0) nu2 <- 1e-10
  }
  diag(R) <- diag(R) + nu2

  # diff.y only needed for non-dsgm non-Gaussian and Gaussian paths
  diff.y <- if (!is_dsgm && object$family != "gaussian") {
               NULL
            }  else if (!is_dsgm) {
              object$y - mu
            }  else NULL

  dast_model <- !is.null(object$power_val)

  # ===========================================================================
  # NON-GAUSSIAN MODELS
  # ===========================================================================
  if (object$family != "gaussian") {

    Sigma     <- par_hat$sigma2 * R
    Sigma_inv <- solve(Sigma)

    if (!obs_loc) {
      A <- if (list_mode) lapply(C, function(Ci) Ci %*% Sigma_inv)
      else C %*% Sigma_inv
    }

    if (is_dsgm) {

      if (messages)
        message(sprintf("Sampling spatial process for DSGM (%s) using Stan...",
                        object$family))

      n_samples_stan <- control_sim$n_sim
      n_warmup_stan  <- control_sim$burnin
      n_chains_stan  <- control_sim$n_chains %||% 1
      n_cores_stan   <- control_sim$n_cores  %||% 1

      use_mda <- isTRUE(object$use_mda) ||
        (object$family == "intprev" &&
           !is.null(object$mda_times) && length(object$mda_times) > 0)

      if (object$family == "intprev") {
        par_stan <- list(beta    = par_hat$beta,
                         k       = par_hat$k,
                         rho     = par_hat$rho,
                         sigma2  = par_hat$sigma2,
                         phi     = par_hat$phi,
                         omega1  = par_hat$omega1 %||% 0.0)
        if (use_mda) {
          par_stan$alpha_W <- par_hat$alpha_W %||% object$fix_alpha_W
          par_stan$gamma_W <- par_hat$gamma_W %||% object$fix_gamma_W
        }
        S_obj <- sample_spatial_process_stan(
          y_prev            = object$prevalence_data,
          intensity_data    = object$intensity_data,
          vary_k            = object$vary_k,
          D                 = object$D,
          coords            = object$coords,
          ID_coords         = object$ID_coords,
          int_mat           = object$int_mat,
          survey_times_data = object$survey_times_data,
          mda_times         = object$mda_times,
          par               = par_stan,
          n_samples         = n_samples_stan,
          n_warmup          = n_warmup_stan,
          n_chains          = n_chains_stan,
          n_cores           = n_cores_stan,
          intensity_family  = if (identical(object$intensity_family, "negbin")) 1L else 0L,
          messages          = messages
        )

      } else {
        # lf_mdiag
        mda_impact <- if (use_mda)
          compute_mda_effect(object$survey_times_data, object$mda_times, object$int_mat,
                             par_hat$alpha_W %||% object$fix_alpha_W,
                             par_hat$gamma_W %||% object$fix_gamma_W, kappa = 1)
        else
          rep(1.0, nrow(object$D))

        S_obj <- sample_spatial_process_stan_lf(
          y_counts   = object$y_counts,
          units_m    = object$units_m,
          which_diag = object$which_diag,
          D          = object$D,
          coords     = object$coords,
          ID_coords  = object$ID_coords,
          par        = list(beta = par_hat$beta, k = par_hat$k, rho = par_hat$rho,
                            gamma_sens = object$gamma_sens,
                            sigma2 = par_hat$sigma2, phi = par_hat$phi),
          mda_impact = mda_impact,
          n_samples  = n_samples_stan,
          n_warmup   = n_warmup_stan,
          n_chains   = n_chains_stan,
          n_cores    = n_cores_stan,
          messages   = messages
        )
      }

      simulation           <- list()
      simulation$samples   <- list()
      simulation$samples$S <- S_obj$S_samples
      n_samples            <- nrow(S_obj$S_samples)

    } else if (dast_model) {
      alpha      <- par_hat$alpha %||% object$fix_alpha
      mda_effect <- compute_mda_effect(object$survey_times_data, object$mda_times,
                                       object$int_mat, alpha, par_hat$gamma,
                                       kappa = object$power_val)
      simulation <- Laplace_sampling_MCMC_dast(
        y = object$y, units_m = object$units_m, mu = mu, Sigma = Sigma,
        sigma2_re = par_hat$sigma2_re, mda_effect = mda_effect,
        ID_coords = object$ID_coords, ID_re = object$ID_re,
        control_mcmc = control_sim, messages = messages)
    } else {
      simulation <- Laplace_sampling_MCMC(
        y = object$y, units_m = object$units_m, mu = mu, Sigma = Sigma,
        sigma2_re = par_hat$sigma2_re, invlink = object$linkf,
        ID_coords = object$ID_coords, ID_re = object$ID_re,
        family = object$family, control_mcmc = control_sim, messages = messages)
    }

    if (obs_loc) {
      out$S_samples <- t(simulation$samples$S)
    } else {
      mu_cond_S <- if (list_mode)
        lapply(A, function(Ai) Ai %*% t(simulation$samples$S))
      else
        A %*% t(simulation$samples$S)

      if (type == "marginal") {
        sd_cond_S <- sqrt(par_hat$sigma2 - diag(A %*% t(C)))
        out$S_samples <- sapply(seq_len(n_samples), function(i)
          mu_cond_S[, i] + sd_cond_S * rnorm(n_pred))

      } else {
        if (list_mode) {
          out$S_samples <- lapply(seq_along(mu_cond_S), function(i) {
            Sp    <- par_hat$sigma2 * matern_cor(dist(grp[[i]]), phi = par_hat$phi,
                                                 kappa = object$kappa, return_sym_matrix = TRUE)
            Sc    <- Sp - A[[i]] %*% t(C[[i]])
            Scr   <- t(chol(Sc))
            sapply(seq_len(n_samples), function(j)
              mu_cond_S[[i]][, j] + Scr %*% rnorm(nrow(mu_cond_S[[i]])))
          })
        } else {
          Sp  <- par_hat$sigma2 * matern_cor(dist(grp), phi = par_hat$phi,
                                             kappa = object$kappa, return_sym_matrix = TRUE)
          Sc  <- Sp - A %*% t(C)
          Scr <- t(chol(Sc))
          out$S_samples <- sapply(seq_len(n_samples), function(i)
            mu_cond_S[, i] + Scr %*% rnorm(nrow(mu_cond_S)))
        }
      }
    }

  } else {
    # =========================================================================
    # GAUSSIAN MODEL (unchanged)
    # =========================================================================
    if (!is.null(object$fix_var_me) && object$fix_var_me > 0 || is.null(object$fix_var_me)) {
      m        <- length(object$y)
      s_unique <- unique(object$ID_coords)
      ID_g     <- as.matrix(cbind(object$ID_coords, object$ID_re))
      n_dim_re_tot <- sapply(seq_len(n_re + 1), function(i) length(unique(ID_g[, i])))
      C_g      <- matrix(0, nrow = m, ncol = sum(n_dim_re_tot))
      for (i in seq_len(m)) {
        ind_s_i <- which(s_unique == ID_g[i, 1])
        C_g[i, seq_len(n_dim_re_tot[1])][ind_s_i] <- 1
      }
      if (n_re > 0) {
        for (j in seq_len(n_re)) {
          select_col <- sum(n_dim_re_tot[seq_len(j)])
          for (i in seq_len(m)) {
            ind_re_j_i <- which(re_unique[[j]] == ID_g[i, j + 1])
            C_g[i, select_col + seq_len(n_dim_re_tot[j + 1])][ind_re_j_i] <- 1
          }
        }
      }
      C_g      <- Matrix(C_g, sparse = TRUE, doDiag = FALSE)
      C_g_m    <- forceSymmetric(Matrix::t(C_g) %*% C_g)
      Sigma_g  <- matrix(0, nrow = sum(n_dim_re_tot), ncol = sum(n_dim_re_tot))
      Sigma_g_inv <- matrix(0, nrow = sum(n_dim_re_tot), ncol = sum(n_dim_re_tot))
      Sigma_g[seq_len(n_dim_re_tot[1]), seq_len(n_dim_re_tot[1])] <- par_hat$sigma2 * R
      Sigma_g_inv[seq_len(n_dim_re_tot[1]), seq_len(n_dim_re_tot[1])] <- solve(R) / par_hat$sigma2
      if (n_re > 0) {
        for (j in seq_len(n_re)) {
          sc <- sum(n_dim_re_tot[seq_len(j)])
          diag(Sigma_g[sc + seq_len(n_dim_re_tot[j+1]), sc + seq_len(n_dim_re_tot[j+1])]) <- par_hat$sigma2_re[j]
          diag(Sigma_g_inv[sc + seq_len(n_dim_re_tot[j+1]), sc + seq_len(n_dim_re_tot[j+1])]) <- 1 / par_hat$sigma2_re[j]
        }
      }
      Sigma_star     <- Sigma_g_inv + C_g_m / par_hat$sigma2_me
      Sigma_star_inv <- forceSymmetric(Matrix::solve(Sigma_star))
      B    <- -C_g %*% Sigma_star_inv %*% Matrix::t(C_g) / (par_hat$sigma2_me^2)
      diag(B) <- Matrix::diag(B) + 1 / par_hat$sigma2_me
      A    <- C %*% B
    } else {
      Sigma     <- par_hat$sigma2 * R
      Sigma_inv <- solve(Sigma)
      A         <- C %*% Sigma_inv
    }

    mu_cond_S <- as.numeric(A %*% diff.y)

    if (type == "marginal") {
      sd_cond_S <- sqrt(par_hat$sigma2 - Matrix::diag(A %*% t(C)))
      out$S_samples <- sapply(seq_len(n_samples), function(i)
        mu_cond_S + sd_cond_S * rnorm(n_pred))
    } else {
      Sp  <- par_hat$sigma2 * matern_cor(dist(grp), phi = par_hat$phi,
                                         kappa = object$kappa, return_sym_matrix = TRUE)
      Sc  <- Sp - A %*% t(C)
      Scr <- t(chol(Sc))
      out$S_samples <- sapply(seq_len(n_samples), function(i)
        mu_cond_S + Scr %*% rnorm(n_pred))
    }
  }

  # ---------------------------------------------------------------------------
  # Random effect posterior samples (unchanged from original)
  # ---------------------------------------------------------------------------
  if (n_re > 0 && !is.null(re_predictors)) {
    out$re          <- list()
    out$re$D_pred   <- D_re_pred
    out$re$samples  <- list()
    re_names        <- colnames(object$ID_re)
    if (object$family == "gaussian") {
      Sigma_cond_inv <- solve(Sigma_cond)
      C_Z  <- C_g[, -(seq_len(n_dim_re_tot[1]))]
      add  <- 0
      for (i in seq_along(n_dim_re_tot[-1])) {
        C_Z[, add + seq_len(n_dim_re_tot[i+1])] <- par_hat$sigma2_re[i] *
          C_Z[, add + seq_len(n_dim_re_tot[i+1])]
        add <- n_dim_re_tot[i+1]
      }
      A_Z          <- Matrix::t(C_Z) %*% B %*% t(C) %*% Sigma_cond_inv
      Sigma_Z_cond <- diag(rep(par_hat$sigma2_re, n_dim_re_tot[-1])) -
        Matrix::t(C_Z) %*% B %*% C_Z -
        A_Z %*% C %*% Matrix::t(B) %*% C_Z
      Scr_Z        <- t(chol(Sigma_Z_cond))
      mu_Z_cond    <- sapply(seq_len(n_samples), function(i)
        as.matrix(A_Z %*% (out$S_samples[, i] - mu_cond_S)))
      add <- 0
      for (i in seq_len(n_re)) {
        for (j in seq_len(n_dim_re_tot[1 + i])) {
          add <- add + 1
          ind_ij   <- which(object$ID_re[[i]] == re_unique[[i]][j])
          C_re_ij  <- matrix(0, ncol = m)
          C_re_ij[, ind_ij] <- par_hat$sigma2_re[i]
          mu_Z_cond[add, ] <- as.numeric(C_re_ij %*% B %*% diff.y) + mu_Z_cond[add, ]
        }
      }
      re_samples <- sapply(seq_len(n_samples), function(i)
        as.numeric(mu_Z_cond[, i] + Scr_Z %*% rnorm(sum(n_dim_re_tot[-1]))))
    } else {
      re_samples <- matrix(0, nrow = sum(n_dim_re_tot[-1]), ncol = n_samples)
      add <- 0
      for (i in seq_len(n_re)) {
        re_samples[seq_len(n_dim_re_tot[i+1]), ] <- t(simulation$samples[[i+1]])
        add <- add + n_dim_re_tot[i+1]
      }
    }
    add <- 0
    for (i in seq_len(n_re)) {
      for (j in seq_len(n_dim_re_tot[1 + i])) {
        add <- add + 1
        out$re$samples[[re_names[[i]]]][[as.character(re_unique[[i]][[j]])]] <- re_samples[add, ]
      }
    }
  } else {
    out$re <- list(D_pred = NULL, samples = NULL)
  }

  if (dast_model || object$family %in% c("intprev", "lf_mdiag")) {
    out$mda_times <- object$mda_times
    if (is.null(par_hat$alpha)) out$fix_alpha <- object$fix_alpha
    out$power_val <- object$power_val
    out$use_mda   <- if (dast_model) {
      !is.null(object$mda_times) && length(object$mda_times) > 0
    } else {
      isTRUE(object$use_mda) ||
        (object$family == "intprev" &&
           !is.null(object$mda_times) && length(object$mda_times) > 0)
    }
  }

  out$obs_loc  <- obs_loc
  if (obs_loc) out$ID_coords <- object$ID_coords
  out$inter_f  <- inter_f
  out$family   <- object$family
  if (is_dsgm) out$vary_k <- isTRUE(object$vary_k)
  # Forward gamma_sens for lf_mdiag so pred_target_grid can read it
  if (is_dsgm && !is.null(object$gamma_sens))
    out$gamma_sens <- object$gamma_sens
  out$par_hat  <- par_hat
  out$cov_offset <- pred_cov_offset
  out$type     <- type
  class(out)   <- "RiskMap.pred.re"
  return(out)
}

##' @title Predictive Target Over a Regular Spatial Grid
##'
##' @description Computes predictions over a regular spatial grid using outputs from
##' \code{\link{pred_over_grid}}.
##'
##' For \strong{STH models} (\code{family = "intprev"}) the default predictive
##' targets are:
##' \describe{
##'   \item{\code{prevalence}}{Probability of observing at least one egg,
##'     \eqn{P(Y>0) = 1 - [k/(k + \mu_W(1-e^{-\rho}))]^k}.}
##'   \item{\code{worm_burden}}{Expected mean worm burden \eqn{\mu_W = e^{\eta}}.}
##'   \item{\code{intensity}}{Mean egg count \eqn{\rho\,\mu_W} (unconditional).}
##' }
##'
##' For \strong{LF models} (\code{family = "lf_mdiag"}) the default predictive
##' targets are:
##' \describe{
##'   \item{\code{mf_prevalence}}{Microfilarial (parasitological) prevalence,
##'     \eqn{P(\text{MF}>0) = 1 - [k/(k + \mu_W(1-e^{-\rho}))]^k}.}
##'   \item{\code{antigen_prevalence}}{Circulating filarial antigen (serological)
##'     prevalence, \eqn{\gamma_s(1 - [k/(k+\mu_W)]^k)}, where \eqn{\gamma_s}
##'     is the serological test sensitivity.}
##'   \item{\code{worm_burden}}{Expected mean worm burden \eqn{\mu_W = e^{\eta}}.}
##' }
##'
##' Custom targets can always be supplied via \code{f_target}.
##'
##' @param object Output from \code{\link{pred_over_grid}}, a
##'   \code{RiskMap.pred.re} object.
##' @param include_covariates Logical. Include covariate effects in the linear
##'   predictor. Default \code{TRUE}.
##' @param include_nugget Logical. Add a nugget draw to each spatial sample.
##'   Default \code{FALSE}.
##' @param include_cov_offset Logical. Include the covariate offset.
##'   Default \code{FALSE}.
##' @param include_mda_effect Logical. Apply the MDA decay to the linear
##'   predictor. Default \code{TRUE}.
##' @param mda_grid Optional. Matrix (or list of matrices for list-mode) of MDA
##'   coverage values at each prediction location and MDA round
##'   (\code{n_pred x n_mda_rounds}). Required when \code{include_mda_effect =
##'   TRUE} and the model was fitted with MDA.
##' @param time_pred Optional. Scalar or vector giving the prediction time(s).
##'   Required when the model includes MDA.
##' @param include_re Logical. Include unstructured random effects.
##'   Default \code{FALSE}.
##' @param f_target Optional named list of functions to apply on the linear
##'   predictor samples. Each function receives a matrix
##'   (\code{n_pred x n_samples}) and returns a matrix of the same dimensions.
##'   Overrides the model-specific defaults described above.
##' @param pd_summary Optional named list of summary functions applied
##'   row-wise to each target matrix (default: mean, median, sd, 2.5\% and
##'   97.5\% quantiles).
##'
##' @return An object of class \code{"RiskMap_pred_target_grid"}.
##' @seealso \code{\link{pred_over_grid}}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @importFrom Matrix solve
##' @export
pred_target_grid <- function(object,
                             include_covariates  = TRUE,
                             include_nugget      = FALSE,
                             include_cov_offset  = FALSE,
                             include_mda_effect  = TRUE,
                             mda_grid            = NULL,
                             time_pred           = NULL,
                             include_re          = FALSE,
                             f_target            = NULL,
                             pd_summary          = NULL) {

  if (!inherits(object, "RiskMap.pred.re"))
    stop("'object' must be an output of pred_over_grid()")

  # ---------------------------------------------------------------------------
  # Model-type flags
  # ---------------------------------------------------------------------------
  sth_model  <- isTRUE(object$family == "intprev")
  lf_model   <- isTRUE(object$family == "lf_mdiag")
  dsgm_model <- sth_model || lf_model          # any DSGM variant
  dast_model <- !dsgm_model && !is.null(object$par_hat$gamma)

  # ---------------------------------------------------------------------------
  # list-mode detection
  # ---------------------------------------------------------------------------
  list_mode <- is.list(object$grid_pred) &&
    !inherits(object$grid_pred, "sfc") &&
    !inherits(object$grid_pred, "sf")

  if (list_mode) {
    n_pred <- vapply(object$grid_pred,
                     function(g) nrow(sf::st_coordinates(g)), integer(1))
  } else {
    n_pred <- nrow(object$S_samples)
  }

  # ---------------------------------------------------------------------------
  # MDA checks: only required when the model was fitted with MDA
  # ---------------------------------------------------------------------------
  use_mda_obj <- isTRUE(object$use_mda) ||
    (sth_model && !is.null(object$mda_times) &&
       length(object$mda_times) > 0 && is.null(object$use_mda))

  needs_mda <- (dsgm_model || dast_model) && include_mda_effect && use_mda_obj

  if (needs_mda) {
    if (is.null(mda_grid))
      stop("'mda_grid' must be supplied when include_mda_effect = TRUE and the model includes MDA")
    if (is.null(time_pred))
      stop("'time_pred' must be supplied when include_mda_effect = TRUE and the model includes MDA")

    if (list_mode) {
      if (!is.list(mda_grid))
        stop("When 'object$grid_pred' is a list, 'mda_grid' must also be a list")
      if (length(mda_grid) != length(object$grid_pred))
        stop("Length of 'mda_grid' must equal length of 'object$grid_pred'")
      for (i in seq_along(mda_grid)) {
        if (!is.matrix(mda_grid[[i]]) && !is.data.frame(mda_grid[[i]]))
          stop(sprintf("'mda_grid[[%d]]' must be a matrix or data.frame", i))
        if (nrow(mda_grid[[i]]) != n_pred[i])
          stop(sprintf("Rows of 'mda_grid[[%d]]' must equal locations in 'object$grid_pred[[%d]]'", i, i))
      }
    }
  }

  if (!is.null(object$par_hat$tau2) == FALSE && include_nugget)
    stop("No nugget was estimated; cannot include it in the predictive target")

  # ---------------------------------------------------------------------------
  # Default f_target — model-specific
  # ---------------------------------------------------------------------------
  if (is.null(f_target)) {

    if (sth_model) {
      k_val       <- object$par_hat$k
      rho_val     <- object$par_hat$rho
      c_rho       <- 1 - exp(-rho_val)
      log_k_val   <- log(k_val)
      vary_k_flag <- isTRUE(object$vary_k) && !is.null(object$par_hat$omega1)
      omega1_val  <- if (vary_k_flag) object$par_hat$omega1 else 0.0

      k_fun <- function(mu_W) {
        if (vary_k_flag)
          exp(log_k_val + omega1_val * log(pmax(mu_W, 1e-10)))
        else
          array(k_val, dim = dim(mu_W))
      }

      f_target <- list(
        prevalence = function(lp) {
          mu_W <- exp(lp)
          k_i  <- k_fun(mu_W)
          pr   <- 1 - (k_i / (k_i + mu_W * c_rho))^k_i
          pmin(pmax(pr, 1e-10), 1 - 1e-10)
        },
        worm_burden = function(lp) {
          exp(lp)
        },
        intensity = function(lp) {
          rho_val * exp(lp)
        }
      )

      if (vary_k_flag) {
        f_target$aggregation <- function(lp) {
          k_fun(exp(lp))
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Default pd_summary
  # ---------------------------------------------------------------------------
  if (is.null(pd_summary)) {
    pd_summary <- list(
      mean   = mean,
      median = median,
      sd     = sd,
      lower  = function(x) quantile(x, 0.025),
      upper  = function(x) quantile(x, 0.975)
    )
  }

  n_f         <- length(f_target)
  n_summaries <- length(pd_summary)
  names_f     <- names(f_target)
  if (is.null(names_f)) names_f <- paste0("f_target_", seq_len(n_f))
  names_s     <- names(pd_summary)
  if (is.null(names_s)) names_s <- paste0("pd_summary_", seq_len(n_summaries))

  if (list_mode) {
    n_samples <- ncol(object$S_samples[[1]])
  } else {
    n_samples <- ncol(object$S_samples)
  }

  n_re <- length(object$re$samples)

  # ---------------------------------------------------------------------------
  # Covariate / offset checks
  # ---------------------------------------------------------------------------
  if (length(object$mu_pred) == 1 && object$mu_pred == 0 && include_covariates)
    stop("Covariates were not provided in pred_over_grid(); rerun with 'predictors'")

  if (n_re == 0 && include_re)
    stop("Random effect categories not provided; rerun pred_over_grid() with 're_predictors'")

  if (list_mode) {
    mu_target  <- if (include_covariates) object$mu_pred  else lapply(n_pred, function(n) rep(0, n))
    cov_offset <- if (include_cov_offset) object$cov_offset else lapply(n_pred, function(n) rep(0, n))
  } else {
    mu_target  <- if (include_covariates) object$mu_pred  else 0
    cov_offset <- if (include_cov_offset) object$cov_offset else 0
  }

  if (include_cov_offset && length(object$cov_offset) == 1)
    stop("No covariate offset was included in the model")

  # ---------------------------------------------------------------------------
  # Optional nugget draw
  # ---------------------------------------------------------------------------
  if (include_nugget) {
    tau2 <- object$par_hat$tau2
    if (list_mode) {
      object$S_samples <- lapply(seq_along(object$S_samples), function(i) {
        object$S_samples[[i]] +
          matrix(rnorm(n_samples * n_pred[i], sd = sqrt(tau2)), ncol = n_samples)
      })
    } else {
      object$S_samples <- object$S_samples +
        matrix(rnorm(n_samples * n_pred, sd = sqrt(tau2)), ncol = n_samples)
    }
  }

  # ---------------------------------------------------------------------------
  # Build linear predictor samples:  lp = mu + offset + S(x)
  # ---------------------------------------------------------------------------
  if (list_mode) {
    object$S_samples <- lapply(object$S_samples, function(x) {
      if (is.numeric(x) && is.vector(x)) matrix(x, nrow = 1) else x
    })

    lp_samples <- vector("list", length(object$grid_pred))
    for (i in seq_along(object$grid_pred)) {
      mu_i <- if (is.matrix(mu_target[[i]])) mu_target[[i]] else mu_target[[i]]
      lp_samples[[i]] <- sapply(seq_len(n_samples), function(j)
        mu_i + cov_offset[[i]] + object$S_samples[[i]][, j])
    }

  } else {
    ID_coords <- if (object$obs_loc) object$ID_coords else seq_len(n_pred)

    lp_samples <- if (is.matrix(mu_target)) {
      sapply(seq_len(n_samples), function(i)
        mu_target[, i] + cov_offset + object$S_samples[ID_coords, i])
    } else {
      sapply(seq_len(n_samples), function(i)
        mu_target + cov_offset + object$S_samples[ID_coords, i])
    }
    # sapply drops dimensions when n_pred == 1; restore matrix shape
    if (!is.matrix(lp_samples))
      lp_samples <- matrix(lp_samples, nrow = length(ID_coords))
  }

  # ---------------------------------------------------------------------------
  # Unstructured random effects
  # ---------------------------------------------------------------------------
  if (include_re) {
    n_dim_re <- sapply(seq_len(n_re), function(i) length(object$re$samples[[i]]))
    for (i in seq_len(n_re)) {
      for (j in seq_len(n_dim_re[i])) {
        re_samp <- object$re$samples[[i]][[j]]   # length n_samples
        if (list_mode) {
          for (g in seq_along(lp_samples))
            lp_samples[[g]] <- lp_samples[[g]] +
              outer(object$re$D_pred[[i]][, j], re_samp)
        } else {
          lp_samples <- lp_samples +
            outer(object$re$D_pred[[i]][, j], re_samp)
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # MDA decay applied on the log scale (DSGM: affects mu_W directly)
  # DAST: applied *after* transformation — handled in the summary loop below
  # ---------------------------------------------------------------------------
  if (dsgm_model && needs_mda) {
    alpha_mda <- object$par_hat$alpha_W %||% object$fix_alpha_W
    gamma_mda <- object$par_hat$gamma_W %||% object$fix_gamma_W

    if (list_mode) {
      for (i in seq_along(lp_samples)) {
        mda_vals <- compute_mda_effect(
          rep(time_pred, n_pred[i]),
          mda_times    = object$mda_times,
          intervention = mda_grid[[i]],
          alpha        = alpha_mda,
          gamma        = gamma_mda,
          kappa        = 1
        )
        lp_samples[[i]] <- lp_samples[[i]] + log(mda_vals)
      }
    } else {
      mda_vals <- compute_mda_effect(
        rep(time_pred, n_pred),
        mda_times    = object$mda_times,
        intervention = mda_grid,
        alpha        = alpha_mda,
        gamma        = gamma_mda,
        kappa        = 1
      )
      lp_samples <- lp_samples + log(mda_vals[ID_coords])
    }
  }

  # ---------------------------------------------------------------------------
  # Apply f_target transformations + summaries
  # ---------------------------------------------------------------------------
  out <- list()
  out$target  <- list()
  if (dsgm_model) out$samples <- list()

  if (list_mode) {
    group_names <- names(object$grid_pred) %||%
      paste0("group_", seq_along(object$grid_pred))

    for (i in seq_along(object$grid_pred)) {
      out$target[[group_names[i]]] <- list()
      if (dsgm_model) out$samples[[group_names[i]]] <- list()

      for (fi in seq_len(n_f)) {
        target_mat <- f_target[[fi]](lp_samples[[i]])

        # DAST: MDA applied multiplicatively after transformation
        if (dast_model && needs_mda) {
          mda_vals <- compute_mda_effect(
            rep(time_pred, n_pred[i]),
            mda_times    = object$mda_times,
            intervention = mda_grid[[i]],
            alpha        = object$par_hat$alpha %||% object$fix_alpha,
            gamma        = object$par_hat$gamma,
            kappa        = object$power_val
          )
          target_mat <- target_mat * mda_vals
        }

        out$target[[group_names[i]]][[names_f[fi]]] <- list()
        for (si in seq_len(n_summaries))
          out$target[[group_names[i]]][[names_f[fi]]][[names_s[si]]] <-
          apply(target_mat, 1, pd_summary[[si]])

        if (dsgm_model)
          out$samples[[group_names[i]]][[names_f[fi]]] <- target_mat
      }
    }

  } else {
    for (fi in seq_len(n_f)) {
      target_mat <- f_target[[fi]](lp_samples)
      if (!is.matrix(target_mat))
        target_mat <- matrix(target_mat, nrow = nrow(lp_samples))

      # DAST: MDA applied multiplicatively after transformation
      if (dast_model && needs_mda) {
        mda_vals <- compute_mda_effect(
          rep(time_pred, n_pred),
          mda_times    = object$mda_times,
          intervention = mda_grid,
          alpha        = object$par_hat$alpha %||% object$fix_alpha,
          gamma        = object$par_hat$gamma,
          kappa        = object$power_val
        )
        target_mat <- target_mat * mda_vals[ID_coords]
      }

      out$target[[names_f[fi]]] <- list()
      for (si in seq_len(n_summaries))
        out$target[[names_f[fi]]][[names_s[si]]] <-
        apply(target_mat, 1, pd_summary[[si]])

      if (dsgm_model)
        out$samples[[names_f[fi]]] <- target_mat
    }
  }

  # ---------------------------------------------------------------------------
  # Metadata
  # ---------------------------------------------------------------------------
  out$grid_pred  <- object$grid_pred
  out$f_target   <- names_f
  out$pd_summary <- names_s
  out$family     <- object$family
  out$lp_samples <- lp_samples

  if (dsgm_model || dast_model) {
    out$mda_effect_applied <- needs_mda
    out$time_pred          <- time_pred
    if (needs_mda) out$mda_grid <- mda_grid
  }

  class(out) <- "RiskMap_pred_target_grid"
  return(out)
}



##' Plot Method for RiskMap_pred_target_grid Objects
##'
##' Generates a plot of the predicted values or summaries over the regular spatial grid
##' from an object of class 'RiskMap_pred_target_grid'.
##'
##' @param x An object of class 'RiskMap_pred_target_grid'.
##' @param which_target Character string specifying which target prediction to plot.
##' @param which_summary Character string specifying which summary statistic to plot (e.g., "mean", "sd").
##' @param ... Additional arguments passed to the \code{\link[terra]{plot}} function of the \code{terra} package.
##' @return A \code{ggplot} object representing the specified prediction target or summary statistic over the spatial grid.
##' @details
##' This function requires the 'terra' package for spatial data manipulation and plotting.
##' It plots the values or summaries over a regular spatial grid, allowing for visual examination of spatial patterns.
##'
##' @seealso \code{\link{pred_target_grid}}
##'
##' @importFrom terra as.data.frame rast plot
##' @method plot RiskMap_pred_target_grid
##' @export
##'
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
plot.RiskMap_pred_target_grid <- function(x, which_target = "linear_target", which_summary = "mean", ...) {
  t_data.frame <-
    terra::as.data.frame(cbind(st_coordinates(x$grid_pred),
                               x$target[[which_target]][[which_summary]]),
                         xy = TRUE)
  raster_out <- terra::rast(t_data.frame, crs = st_crs(x$grid_pred)$input)

  terra::plot(raster_out, ...)
}

##' @title Predictive Targets over a Shapefile (grid-aggregated)
##'
##' @description
##' Computes predictive targets over polygon features using joint prediction
##' samples from \code{\link{pred_over_grid}}. Targets can incorporate
##' covariates, offsets, optional unstructured random effects, and (if fitted)
##' mass drug administration (MDA) effects from a DAST model.
##'
##' @param object Output from \code{\link{pred_over_grid}} (class \code{RiskMap.pred.re}),
##'   typically fitted with \code{type = "joint"} so that linear predictor samples are available.
##' @param shp An \pkg{sf} polygon object (preferred) or a \code{data.frame} with an
##'   attached geometry column, representing regions over which predictions are aggregated.
##' @param shp_target A function that aggregates grid-cell values within each polygon to a
##'   single regional value (default \code{mean}). Examples: \code{mean}, \code{sum},
##'   a custom weighted mean, etc.
##' @param weights Optional numeric vector of weights used inside \code{shp_target}.
##'   If supplied with \code{standardize_weights = TRUE}, weights are normalized within each region.
##' @param standardize_weights Logical; standardize \code{weights} within each region (\code{FALSE} by default).
##' @param col_names Name or column index in \code{shp} containing region identifiers to use in outputs.
##' @param include_covariates Logical; include fitted covariate effects in the linear predictor (default \code{TRUE}).
##' @param include_nugget Logical; include the nugget (unstructured measurement error) in the linear predictor (default \code{FALSE}).
##' @param include_cov_offset Logical; include any covariate offset term (default \code{FALSE}).
##' @param include_mda_effect Logical; include the MDA effect as defined by the fitted DAST model
##'   (default \code{TRUE}). Requires \code{time_pred} and, when applicable, \code{mda_grid}.
##' @param return_shp Logical; if \code{TRUE}, return the shapefile with appended summary columns
##'   defined by \code{pd_summary} (default \code{TRUE}).
##' @param time_pred Optional numeric scalar (or time index) at which to evaluate the predictive target
##!'   when MDA effects are included. If \code{NULL}, uses the default time implied by \code{object}.
##' @param mda_grid Optional structure describing MDA schedules aligned with prediction grid cells
##'   (e.g., a \code{data.frame}/matrix/list). Used only when \code{include_mda_effect = TRUE}.
##' @param include_re Logical; include unstructured random effects (RE) in the linear predictor (default \code{FALSE}).
##' @param f_target List of target functions applied to linear predictor samples (e.g.,
##'   \code{list(prev = plogis)} for prevalence on the probability scale). If \code{NULL},
##'   the identity is used.
##' @param pd_summary Named list of summary functions applied to each region's target samples
##'   (e.g., \code{list(mean = mean, sd = sd, q025 = function(x) quantile(x, 0.025), q975 = function(x) quantile(x, 0.975))}).
##'   Names are used as column suffixes in the outputs.
##' @param messages Logical; if \code{TRUE}, print progress messages while computing regional targets.
##' @param return_target_samples Logical; if \code{TRUE}, also return the raw target samples per region
##'   (default \code{FALSE}).
##'
##' @details
##' For each polygon in \code{shp}, grid-cell samples of the linear predictor are transformed with
##' \code{f_target}, optionally adjusted for covariates, offset, nugget, MDA effects and/or REs, and
##' then aggregated via \code{shp_target} (optionally weighted). The list \code{pd_summary} is applied
##' to each region's target samples to produce summary statistics.
##'
##' @return An object of class \code{RiskMap_pred_target_shp} with components:
##' \itemize{
##'   \item \code{target}: \code{data.frame} of region-level summaries (one row per region).
##'   \item \code{target_samples}: (optional) \code{list} with one element per region; each contains
##'         a \code{data.frame}/matrix of raw samples for each named target in \code{f_target},
##'         if \code{return_target_samples = TRUE}.
##'   \item \code{shp}: (optional) the input \code{sf} object with appended summary columns,
##'         included if \code{return_shp = TRUE}.
##'   \item \code{f_target}, \code{pd_summary}, \code{grid_pred}: inputs echoed for reproducibility.
##' }
##'
##' @seealso \code{\link{pred_over_grid}}, \code{\link{pred_target_grid}}
##'
##' @importFrom terra rast as.data.frame
##' @importFrom stats plogis
##' @export
pred_target_shp <- function(object, shp, shp_target = mean,
                            weights = NULL, standardize_weights = FALSE,
                            col_names = NULL,
                            include_covariates = TRUE,
                            include_nugget = FALSE,
                            include_cov_offset = FALSE,
                            include_mda_effect = TRUE,
                            return_shp = TRUE,
                            time_pred = NULL,
                            mda_grid = NULL,
                            include_re = FALSE,
                            f_target = NULL,
                            pd_summary = NULL,
                            messages = TRUE,
                            return_target_samples = FALSE) {

  if(!inherits(object, what = "RiskMap.pred.re", which = FALSE)) {
    stop("The object passed to 'object' must be an output of
         the function 'pred_S'")
  }

  if(object$type != "joint") {
    stop("To run predictions with a shape file, joint predictions must be used;
         rerun 'pred_over_grid' and set type='joint'")
  }

  dast_model <- !is.null(object$par_hat$gamma)

  if(dast_model) {
    if(include_mda_effect & is.null(mda_grid)) {
      stop("The MDA coverage must be specified for each point on the grid through the argument 'mda_grid'")
    }
    if(is.null(time_pred)) {
      stop("For a DAST model, the time of prediction must be specified through the argument 'time_pred'")
    }
  }

  list_mode <- is.list(object$grid_pred) & !(inherits(object$grid_pred,"sfc") |
                                               inherits(object$grid_pred,"sf"))

  if (list_mode) {
    ok_geom <- vapply(
      object$grid_pred,
      function(g) {
        inherits(g, "sf") || inherits(g, "sfc")
      },
      logical(1)
    )
    if (!all(ok_geom)) {
      stop("When 'object$grid_pred' is a list, each element must be an 'sf' or 'sfc' object.")
    }

    n_pred <- vapply(object$grid_pred, function(g) nrow(sf::st_coordinates(g)), integer(1))

    if (!is.null(weights)) {
      if (!is.list(weights)) {
        stop("When 'object$grid_pred' is a list, 'weights' must also be a list, ",
             "with one numeric vector per element of 'object$grid_pred'.")
      }
      if (length(weights) != length(object$grid_pred)) {
        stop("Length of 'weights' must match length of 'object$grid_pred'.")
      }
      for (i in seq_along(weights)) {
        if (!is.numeric(weights[[i]])) {
          stop(sprintf("'weights[[%d]]' must be numeric.", i))
        }
        if (length(weights[[i]]) != n_pred[i]) {
          stop(sprintf("Length of 'weights[[%d]]' (%d) must equal number of locations in 'object$grid_pred[[%d]]' (%d).",
                       i, length(weights[[i]]), i, n_pred[i]))
        }
        if (anyNA(weights[[i]])) {
          warning(sprintf("NA values found in 'weights[[%d]]'; they will be treated as 0.", i))
          weights[[i]][is.na(weights[[i]])] <- 0
        }
      }
      no_weights <- FALSE
    } else {
      weights <- lapply(n_pred, function(n) rep(1, n))
      no_weights <- TRUE
    }

    if (dast_model && include_mda_effect) {
      if (is.null(mda_grid)) {
        stop("With DAST and include_mda_effect = TRUE, 'mda_grid' must be provided (as a list matching 'object$grid_pred').")
      }
      if (!is.list(mda_grid)) {
        stop("When 'object$grid_pred' is a list, 'mda_grid' must also be a list (one element per group).")
      }
      if (length(mda_grid) != length(object$grid_pred)) {
        stop("Length of 'mda_grid' must match length of 'object$grid_pred'.")
      }
      for (i in seq_along(mda_grid)) {
        if (!is.matrix(mda_grid[[i]]) && !is.data.frame(mda_grid[[i]])) {
          stop(sprintf("'mda_grid[[%d]]' must be a matrix or data.frame.", i))
        }
        if (nrow(mda_grid[[i]]) != n_pred[i]) {
          stop(sprintf("Rows of 'mda_grid[[%d]]' (%d) must equal number of locations in 'object$grid_pred[[%d]]' (%d).",
                       i, nrow(mda_grid[[i]]), i, n_pred[i]))
        }
      }
    }
  } else {
    n_pred <- nrow(st_coordinates(object$grid_pred))
  }

  n_re <- length(object$re$samples)
  re_names <- names(object$re$samples)

  if(n_re == 0 && include_re) {
    stop("The categories of the randome effects variables have not been provided;
         re-run pred_over_grid and provide the covariates through the argument 're_predictors'")
  }

  if(!is.null(weights)) {
    no_weights <- FALSE
  } else {
    if(list_mode) {
      weights <- sapply(unlist(n_pred), function(i) rep(1, i))
    } else {
      weights <- rep(1, n_pred)
    }
    no_weights <- TRUE
  }

  if(is.null(object$par_hat$tau2) & include_nugget) {
    stop("The nugget cannot be included in the predictive target
             because it was not included when fitting the model")
  }

  if(is.null(f_target)) {
    f_target <- list(linear_target = function(x) x)
  }

  if(is.null(pd_summary)) {
    pd_summary <- list(mean = mean, sd = sd)
  }

  n_summaries <- length(pd_summary)
  n_f <- length(f_target)
  if(list_mode) {
    n_samples <- ncol(object$S_samples[[1]])
  } else {
    n_samples <- ncol(object$S_samples)
  }

  if(!list_mode & (length(object$mu_pred) == 1 && object$mu_pred == 0 && include_covariates)) {
    stop("Covariates have not been provided; re-run pred_over_grid
         and provide the covariates through the argument 'predictors'")
  }

  if(!include_covariates) {
    mu_target <- 0
  } else {
    if(is.null(object$mu_pred)) stop("the output obtained from 'pred_S' does not
                                     contain any covariates; if including covariates
                                     in the predictive target these shuold be included
                                     when running 'pred_S'")
    mu_target <- object$mu_pred
  }

  if(!include_cov_offset) {
    if(list_mode) {
      cov_offset <- sapply(unlist(n_pred), function(i) rep(0, i))
    } else {
      cov_offset <- 0
    }
  } else {
    if(length(object$cov_offset) == 1) {
      stop("No covariate offset was included in the model;
           set include_cov_offset = FALSE, or refit the model and include
           the covariate offset")
    }
    cov_offset <- object$cov_offset
  }

  if(include_nugget) {
    Z_sim <- matrix(rnorm(n_samples * n_pred, sd = sqrt(object$par_hat$tau2)),
                    ncol = n_samples)
    object$S_samples <- object$S_samples + Z_sim
  }

  out <- list()
  if (return_target_samples) out$target_samples <- list()

  if(list_mode) {
    object$S_samples <- lapply(object$S_samples, function(x) {
      if (is.numeric(x) && is.vector(x)) {
        matrix(x, nrow = 1)
      } else {
        x
      }
    })

    if(is.matrix(mu_target[[1]])) {
      out$lp_samples <- sapply(1:length(object$grid_pred), function(j)
        sapply(1:n_samples,
               function(i)
                 mu_target[[j]][, i] + cov_offset[[j]] +
                 object$S_samples[[j]][, i]))
    } else {
      out$lp_samples <- sapply(1:length(object$grid_pred), function(j)
        sapply(1:n_samples,
               function(i)
                 mu_target[[j]] + cov_offset[[j]] +
                 object$S_samples[[j]][, i]))
    }
  } else {
    if(is.matrix(mu_target)) {
      out$lp_samples <- sapply(1:n_samples, function(i)
        mu_target[, i] + cov_offset + object$S_samples[, i])
    } else {
      out$lp_samples <- sapply(1:n_samples, function(i)
        mu_target + cov_offset + object$S_samples[, i])
    }
  }

  if(include_re) {
    n_dim_re <- sapply(1:n_re, function(i) length(object$re$samples[[i]]))
    for(i in 1:n_re) {
      for(j in 1:n_dim_re[i]) {
        for(h in 1:n_samples) {
          out$lp_samples[, h] <- out$lp_samples[, h] +
            object$re$D_pred[[i]][, j] * object$re$samples[[i]][[j]][h]
        }
      }
    }
  }

  names_f <- names(f_target)
  names_s <- names(pd_summary)
  out$target <- list()

  n_reg <- nrow(shp)
  if(is.null(col_names)) {
    shp$region <- paste("reg", 1:n_reg, sep = "")
    col_names <- "region"
    names_reg <- shp$region
  } else {
    names_reg <- shp[[col_names]]
    if(n_reg != length(names_reg)) {
      stop("The names in the column identified by 'col_names' do not
         provide a unique set of names, but there are duplicates")
    }
  }

  if(list_mode) {
    shp <- st_transform(shp, crs = st_crs(object$grid_pred[[1]])$input)
  } else {
    shp <- st_transform(shp, crs = st_crs(object$grid_pred)$input)
  }

  if(!list_mode) {
    inter <- st_intersects(shp, object$grid_pred)
    if(any(is.na(weights))) {
      warning("Missing values found in 'weights' are set to 0 \n")
      weights[is.na(weights)] <- 0
    }
  } else {
    for(i in 1:length(object$grid_pred)) {
      if(any(is.na(weights[[i]]))) {
        warning("Missing values found in 'weights' are set to 0 \n")
        weights[[i]][is.na(weights[[i]])] <- 0
      }
    }
  }

  no_comp <- NULL
  for(h in 1:n_reg) {

    if(list_mode) {
      if(messages) message("Computing predictive target for: ", shp[[col_names]][h])
      if(standardize_weights & !no_weights) {
        weights_h <- weights[[h]] / sum(weights[[h]])
      } else {
        weights_h <- weights[[h]]
      }
      for(i in 1:n_f) {
        target_grid_samples_i <- as.matrix(f_target[[i]](out$lp_samples[[h]]))
        if(dast_model && include_mda_effect) {
          alpha <- object$par_hat$alpha
          if(is.null(alpha)) alpha <- object$fix_alpha
          gamma <- object$par_hat$gamma
          mda_effect_time_pred <- compute_mda_effect(
            rep(time_pred, n_pred[h]),
            mda_times = object$mda_times,
            intervention = mda_grid[[h]],
            alpha = alpha, gamma = gamma, kappa = object$power_val
          )
          target_grid_samples_i <- target_grid_samples_i * mda_effect_time_pred
        }

        target_samples_i <- apply(target_grid_samples_i, 2, function(x) shp_target(weights_h * x))

        if (return_target_samples) {
          if (is.null(out$target_samples[[ names_reg[h] ]])) out$target_samples[[ names_reg[h] ]] <- list()
          out$target_samples[[ names_reg[h] ]][[ names_f[i] ]] <- target_samples_i
        }

        out$target[[ paste(names_reg[h]) ]][[ paste(names_f[i]) ]] <- list()
        for(j in 1:n_summaries) {
          out$target[[ paste(names_reg[h]) ]][[ paste(names_f[i]) ]][[ paste(names_s[j]) ]] <-
            pd_summary[[j]](target_samples_i)
        }
      }
      if(messages) message(" \n")

    } else {
      if(messages) message("Computing predictive target for:", shp[[col_names]][h])
      if(length(inter[[h]]) == 0) {
        warning(paste("No points on the grid fall within", shp[[col_names]][h],
                      "and no predictions are carried out for this area"))
        no_comp <- c(no_comp, h)
      } else {
        ind_grid_h <- inter[[h]]
        if(standardize_weights & !no_weights) {
          weights_h <- weights[ind_grid_h] / sum(weights[ind_grid_h])
        } else {
          weights_h <- weights[ind_grid_h]
        }
        for(i in 1:n_f) {
          target_grid_samples_i <- as.matrix(f_target[[i]](out$lp_samples[ind_grid_h, ]))
          if(dast_model && include_mda_effect) {
            alpha <- object$par_hat$alpha
            if(is.null(alpha)) alpha <- object$fix_alpha
            gamma <- object$par_hat$gamma
            mda_effect_time_pred <- compute_mda_effect(
              rep(time_pred, length(ind_grid_h)),
              mda_times = object$mda_times,
              intervention = mda_grid[ind_grid_h, ],
              alpha = alpha, gamma = gamma, kappa = object$power_val
            )
            target_grid_samples_i <- target_grid_samples_i * mda_effect_time_pred
          }

          target_samples_i <- apply(target_grid_samples_i, 2, function(x) shp_target(weights_h * x))

          if (return_target_samples) {
            if (is.null(out$target_samples[[ names_reg[h] ]])) out$target_samples[[ names_reg[h] ]] <- list()
            out$target_samples[[ names_reg[h] ]][[ names_f[i] ]] <- target_samples_i
          }

          out$target[[ paste(names_reg[h]) ]][[ paste(names_f[i]) ]] <- list()
          for(j in 1:n_summaries) {
            out$target[[ paste(names_reg[h]) ]][[ paste(names_f[i]) ]][[ paste(names_s[j]) ]] <-
              pd_summary[[j]](target_samples_i)
          }
        }
      }
      if(messages) message(" \n")
    }
  }

  if(return_shp) {
    if(length(no_comp) > 0) {
      ind_reg <- (1:n_reg)[-no_comp]
    } else {
      ind_reg <- 1:n_reg
    }
    for(i in 1:n_f) {
      for(j in 1:n_summaries) {
        name_ij <- paste(names_f[i], "_", paste(names_s[j]), sep = "")
        shp[[name_ij]] <- rep(NA, n_reg)
        for(h in ind_reg) {
          which_reg <- which(shp[[col_names]] == names_reg[h])
          shp[which_reg, ][[name_ij]] <-
            out$target[[ paste(names_reg[h]) ]][[ paste(names_f[i]) ]][[ paste(names_s[j]) ]]
        }
      }
    }
  }

  out$shp <- shp
  out$f_target <- names(f_target)
  out$pd_summary <- names(pd_summary)
  out$grid_pred <- object$grid_pred
  class(out) <- "RiskMap_pred_target_shp"
  return(out)
}


##' Plot Method for RiskMap_pred_target_shp Objects
##'
##' @description
##' Generates a plot of predictive target values or summaries over a shapefile.
##'
##' @param x An object of class 'RiskMap_pred_target_shp' containing computed targets,
##' summaries, and associated spatial data.
##' @param which_target Character indicating the target type to plot (e.g., "linear_target").
##' @param which_summary Character indicating the summary type to plot (e.g., "mean", "sd").
##' @param ... Additional arguments passed to 'scale_fill_distiller' in 'ggplot2'.
##' @return A \code{ggplot} object showing the plot of the specified predictive target or summary.
##' @details
##' This function plots the predictive target values or summaries over a shapefile.
##' It requires the 'ggplot2' package for plotting and 'sf' objects for spatial data.
##'
##' @seealso
##' \code{\link{pred_target_shp}}, \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_sf}},
##' \code{\link[ggplot2]{aes}}, \code{\link[ggplot2]{scale_fill_distiller}}
##'
##' @importFrom ggplot2 ggplot geom_sf aes scale_fill_distiller
##' @method plot RiskMap_pred_target_shp
##' @export
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
plot.RiskMap_pred_target_shp <- function(x, which_target = "linear_target",
                                         which_summary = "mean", ...) {
  col_shp_name <- paste(which_target,"_",which_summary,sep="")

  out <- ggplot(x$shp) +
    geom_sf(aes(fill = x$shp[[col_shp_name]])) +
    scale_fill_distiller(...)
  return(out)
}

##' @title Update Predictors for a RiskMap Prediction Object
##'
##' @description
##' This function updates the predictors of a given RiskMap prediction object. It ensures that the new predictors match the original prediction grid and updates the relevant components of the object accordingly.
##'
##' @param object A `RiskMap.pred.re` object, which is the output of the \code{\link{pred_over_grid}} function.
##' @param predictors A data frame containing the new predictor values. The number of rows must match the prediction grid in the `object`.
##'
##' @details
##' The function performs several checks and updates:
##' \itemize{
##'   \item Ensures that `object` is of class `RiskMap.pred.re`.
##'   \item Ensures that the number of rows in `predictors` matches the prediction grid in `object`.
##'   \item Removes any rows with missing values in `predictors` and updates the corresponding components of the `object`.
##'   \item Updates the prediction locations, the predictive samples for the random effects, and the linear predictor.
##' }
##'
##' @return The updated `RiskMap.pred.re` object.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @export
update_predictors <- function(object, predictors) {
  if (!inherits(object, what = "RiskMap.pred.re", which = FALSE)) {
    stop("The object passed to 'object' must be an output of the function 'glgpm'")
  }

  list_mode <- is.list(object$grid_pred) && !is.null(object$grid_pred) &
    !any(class(object$grid_pred)=="sf" | class(object$grid_pred)=="sfc")
  par_hat <- object$par_hat
  p <- length(par_hat$beta)

  if (p == 1) {
    stop("No update of the predictors can be done for an intercept-only model")
  }

  inter_f <- object$inter_f
  inter_lt_f <- inter_f
  inter_lt_f$pf <- update(inter_lt_f$pf, NULL ~ .)

  if (list_mode) {
    if (!is.list(predictors)) stop("If object$grid_pred is a list, then 'predictors' must also be a list")
    if (length(predictors) != length(object$grid_pred)) stop("Length of 'predictors' list must match length of 'grid_pred' list")

    n_pred_list <- lapply(object$grid_pred, function(x) nrow(st_coordinates(x)))
    for (i in seq_along(predictors)) {
      if (!is.data.frame(predictors[[i]])) stop(sprintf("predictors[[%d]] must be a data.frame", i))
      if (nrow(predictors[[i]]) != n_pred_list[i]) stop(sprintf("predictors[[%d]] does not match the number of prediction locations", i))
    }

    if (!is.null(object$re_predictors)) {
      for (i in seq_along(predictors)) {
        comb_pred <- data.frame(predictors[[i]], object$re_predictors[[i]])
        ind_c <- complete.cases(comb_pred)

        predictors_aux <- predictors[[i]][ind_c, , drop = FALSE]
        predictors[[i]] <- predictors_aux

        re_predictors_aux <- object$re_predictors[[i]][ind_c, , drop = FALSE]
        object$re_predictors[[i]] <- re_predictors_aux

        n_re <- length(object$re$samples)
        for (r in 1:n_re) {
          object$re$D_pred[[r]][[i]] <- object$re$D_pred[[r]][[i]][ind_c, , drop = FALSE]
        }

        object$grid_pred[[i]] <- object$grid_pred[[i]][ind_c, ]
        if (is.matrix(object$S_samples[[i]])) {
          object$S_samples[[i]] <- object$S_samples[[i]][ind_c, , drop = FALSE]
        } else {
          object$S_samples[[i]] <- object$S_samples[[i]][ind_c]
        }
      }
    } else {
      for (i in seq_along(predictors)) {
        comb_pred <- predictors[[i]]
        ind_c <- complete.cases(comb_pred)
        predictors[[i]] <- predictors[[i]][ind_c, , drop = FALSE]

        object$grid_pred[[i]] <- object$grid_pred[[i]][ind_c, ]
        if (is.matrix(object$S_samples[[i]])) {
          object$S_samples[[i]] <- object$S_samples[[i]][ind_c, , drop = FALSE]
        } else {
          object$S_samples[[i]] <- object$S_samples[[i]][ind_c]
        }
      }
    }

    object$mu_pred <- vector("list", length(predictors))
    for (i in seq_along(predictors)) {
      mf_pred <- model.frame(inter_lt_f$pf, data = predictors[[i]], na.action = na.fail)
      D_pred <- as.matrix(model.matrix(attr(mf_pred, "terms"), data = predictors[[i]]))
      if (ncol(D_pred) != p) stop("Mismatch in number of predictors")
      object$mu_pred[[i]] <- as.numeric(D_pred %*% par_hat$beta)
    }
  } else {
    grp <- st_coordinates(object$grid_pred)
    n_pred <- nrow(grp)

    if (!is.data.frame(predictors)) stop("'predictors' must be an object of class 'data.frame'")
    if (nrow(predictors) != n_pred) stop("the values provided for 'predictors' do not match the prediction grid")

    if (any(is.na(predictors))) {
      warning("There are missing values in 'predictors'; these values have been removed alongside the corresponding prediction locations")
    }

    if (!is.null(object$re_predictors)) {
      comb_pred <- data.frame(predictors, object$re_predictors)
      ind_c <- complete.cases(comb_pred)
      predictors <- predictors[ind_c, , drop = FALSE]
      object$re_predictors <- object$re_predictors[ind_c, , drop = FALSE]

      n_re <- length(object$re$samples)
      for (i in 1:n_re) {
        object$re$D_pred[[i]] <- object$re$D_pred[[i]][ind_c, , drop = FALSE]
      }
    } else {
      comb_pred <- predictors
      ind_c <- complete.cases(comb_pred)
      predictors <- predictors[ind_c, , drop = FALSE]
    }

    object$grid_pred <- object$grid_pred[ind_c, ]
    grp <- st_coordinates(object$grid_pred)
    n_pred <- nrow(grp)
    object$S_samples <- object$S_samples[ind_c, , drop = FALSE]

    mf_pred <- model.frame(inter_lt_f$pf, data = predictors, na.action = na.fail)
    D_pred <- as.matrix(model.matrix(attr(mf_pred, "terms"), data = predictors))
    if (ncol(D_pred) != p) stop("Mismatch in number of predictors")
    object$mu_pred <- as.numeric(D_pred %*% par_hat$beta)
  }

  class(object) <- "RiskMap.pred.re"
  return(object)
}



##' @title Assess Predictive Performance via Spatial Cross-Validation
##'
##' @description
##' This function evaluates the predictive performance of spatial models fitted to `RiskMap` objects using cross-validation. It supports two classes of diagnostic tools:
##'
##' - **Scoring rules**, including the Continuous Ranked Probability Score (CRPS) and its scaled version (SCRPS), which quantify the sharpness and calibration of probabilistic forecasts;
##' - **Calibration diagnostics**, based on the Probability Integral Transform (PIT) for Gaussian outcomes and Aggregated nonparametric PIT (AnPIT) curves for discrete outcomes (e.g., Poisson or Binomial).
##'
##' Cross-validation can be performed using either spatial clustering or regularized subsampling with a minimum inter-point distance. For each fold or subset, models can be refitted or evaluated with fixed parameters, offering flexibility in model validation. The function also provides visualizations of the spatial distribution of test folds.
##'
##' @param object A list of `RiskMap` objects, each representing a model fitted with `glgpm`.
##' @param keep_par_fixed Logical; if `TRUE`, parameters are kept fixed across folds, otherwise the model is re-estimated for each fold.
##' @param iter Integer; number of times to repeat the cross-validation.
##' @param fold Integer; number of folds for cross-validation (required if `method = "cluster"`).
##' @param n_size Optional; the size of the test set, required if `method = "regularized"`.
##' @param control_sim Control settings for simulation, an output from `set_control_sim`.
##' @param method Character; either `"cluster"` or `"regularized"` for the cross-validation method. The `"cluster"` method uses
##' spatial clustering as implemented by the \code{spatial_clustering_cv} function from the `spatialEco` package, while the `"regularized"` method
##' selects a subsample of the dataset by imposing a minimum distance, set by the `min_dist` argument, for a randomly selected
##' subset of locations.
##' @param min_dist Optional; minimum distance for regularized subsampling (required if `method = "regularized"`).
##' @param plot_fold Logical; if `TRUE`, plots each fold's test set.
##' @param messages Logical; if `TRUE`, displays progress messages.
##' @param which_metric Character vector; one or more of `"CRPS"`, `"SCRPS"`, or `"AnPIT"`, to specify the predictive performance metrics to compute.
##' @param user_split A user-defined cross-validation split. Either:
##'   * a matrix with \code{nrow = n} (number of observations) and
##'     \code{ncol = iter} (number of iterations), where entries of \code{1}
##'     indicate membership in the test set for that iteration and \code{0}
##'     indicate training set; or
##'   * a list of length \code{iter}, where each element is either a vector of
##'     test indices, or a list with components \code{in_id} (training indices)
##'     and \code{out_id} (test indices).
##'   When supplied, \code{user_split} overrides the automatic clustering or
##'   regularized distance splitting defined by \code{method}.
##' @param ... Additional arguments passed to clustering or subsampling functions.
##'
##' @return A list of class `RiskMap.spatial.cv`, containing:
##' \describe{
##'   \item{test_set}{A list of test sets used for validation, each of class `'sf'`.}
##'   \item{model}{A named list, one per model, each containing:
##'     \describe{
##'       \item{score}{A list with CRPS and/or SCRPS scores for each fold if requested.}
##'       \item{PIT}{(if `family = "gaussian"` and `which_metric` includes `"AnPIT"`) A list of PIT values for test data.}
##'       \item{AnPIT}{(if `family` is discrete and `which_metric` includes `"AnPIT"`) A list of AnPIT curves for test data.}
##'     }
##'   }
##' }
##'
##' @seealso \code{\link{plot_AnPIT}}
##'
##' @references
##' Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of proper scoring rules. *Statistical Science*, 38(1), 140–159. \doi{10.1214/22-STS864}.
##'
##' @importFrom terra match
##' @importFrom ggplot2 ggplot geom_sf theme_minimal ggtitle
##' @importFrom gridExtra grid.arrange
##' @importFrom stats ecdf integrate rbinom rpois
##' @importFrom spatialEco subsample.distance
##' @importFrom spatialsample spatial_clustering_cv autoplot
##' @importFrom sf st_as_sfc
##' @export
##' @author Emanuele Giorgi
assess_pp <- function(object,
                      keep_par_fixed = TRUE,
                      iter = 1,
                      fold = NULL, n_size = NULL,
                      control_sim = set_control_sim(),
                      method,
                      min_dist = NULL,
                      plot_fold = TRUE,
                      messages = TRUE,
                      which_metric = c("AnPIT", "CRPS", "SCRPS"),
                      user_split = NULL,
                      ...) {

  ## ─────────────────────────── helpers ─────────────────────────── ##
  is_list_of_riskmap <- function(x) {
    is.list(x) && all(vapply(x, inherits, logical(1), what = "RiskMap"))
  }
  is_dast_fit <- function(fit) {
    as.character(fit$call[[1]])=="dast"
  }

  is_dsgm_fit <- function(fit) {
    as.character(fit$call[[1]])=="dsgm"
  }

  crps_gaussian <- function(y, mu, sigma) {
    if (sigma == 0) return(0)
    z <- (y - mu) / sigma
    2 * stats::dnorm(z) + z * (2 * stats::pnorm(z) - 1) - 1 / sqrt(pi)
  }
  crps_discrete <- function(y, Fk) {
    k <- seq_along(Fk) - 1
    sum((Fk - as.numeric(k >= y))^2)
  }
  exp_crps_discrete <- function(Fk, pk) {
    k <- seq_along(pk) - 1
    sum(pk * vapply(k, crps_discrete, numeric(1), Fk = Fk))
  }

  if(!is.null(user_split)) {
    iter <- ncol(user_split)
  }
  ## ────────────────────── sanity checks (unchanged) ─────────────────────── ##
  if (!is_list_of_riskmap(object))
    stop("`object` must be a list of fitted models of class 'RiskMap'.")

  if (!all(which_metric %in% c("CRPS", "SCRPS", "AnPIT")))
    stop("`which_metric` must only contain 'CRPS', 'SCRPS' or 'AnPIT'")

  if (is.null(user_split)) {
    if (!method %in% c("cluster", "regularized"))
      stop("`method` must be either 'cluster' or 'regularized' (unless `user_split` is supplied).")

    if (method == "regularized") {
      if (is.null(min_dist)) stop("for 'regularized', supply `min_dist`")
      if (is.null(n_size))   stop("for 'regularized', supply `n_size`")
    }
    if (method == "cluster" && is.null(fold))
      stop("for 'cluster', supply `fold`")
  }

  if (!inherits(control_sim, "mcmc.RiskMap"))
    stop("`control_sim` must come from `set_control_sim()`")

  get_CRPS  <- "CRPS"  %in% which_metric
  get_SCRPS <- "SCRPS" %in% which_metric
  get_AnPIT <- "AnPIT" %in% which_metric

  ## ───────────────────────────── data & splits ───────────────────────────── ##
  object1 <- object[[1]]
  data_sf <- object1$data_sf
  n_obs   <- nrow(data_sf)

  make_splits_from_user <- function(usr, n_iter_expected) {
    spl <- vector("list", n_iter_expected)
    if (is.matrix(usr)) {
      if (nrow(usr) != n_obs)
        stop("`user_split` matrix must have nrow == nrow(data).")
      if (ncol(usr) != n_iter_expected)
        stop("`user_split` matrix must have ncol == `iter`.")
      for (i in seq_len(n_iter_expected)) {
        out_id <- which(usr[, i] != 0 & !is.na(usr[, i]))
        in_id  <- setdiff(seq_len(n_obs), out_id)
        spl[[i]] <- list(in_id = in_id, out_id = out_id,
                         data = data_sf[in_id, ],
                         data_test = data_sf[out_id, ])
      }
    } else if (is.list(usr)) {
      if (length(usr) != n_iter_expected)
        stop("`user_split` list must have length == `iter`.")
      for (i in seq_len(n_iter_expected)) {
        ui <- usr[[i]]
        if (is.list(ui) && !is.null(ui$in_id) && !is.null(ui$out_id)) {
          in_id  <- ui$in_id
          out_id <- ui$out_id
        } else if (is.integer(ui) || is.double(ui)) {
          out_id <- as.integer(ui)
          in_id  <- setdiff(seq_len(n_obs), out_id)
        } else {
          stop("Each element of `user_split` must be a vector of test indices or a list(in_id=..., out_id=...).")
        }
        spl[[i]] <- list(in_id = in_id, out_id = out_id,
                         data = data_sf[in_id, ],
                         data_test = data_sf[out_id, ])
      }
    } else {
      stop("`user_split` must be a matrix (nrow=n, ncol=iter) or a list.")
    }
    list(splits = spl)
  }

  if (!is.null(user_split)) {
    data_split <- make_splits_from_user(user_split, iter)
    n_iter <- iter

    if (isTRUE(plot_fold)) {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("plot_fold = TRUE requires the 'ggplot2' package; skipping plots.", call. = FALSE)
      } else {
        if (n_iter == 1) {
          p <- ggplot2::ggplot(data_split$splits[[1]]$data_test) +
            ggplot2::geom_sf() +
            ggplot2::theme_minimal() +
            ggplot2::ggtitle("Test set")
          print(p)
        } else {
          plots <- lapply(seq_len(n_iter), function(i) {
            ggplot2::ggplot(data_split$splits[[i]]$data_test) +
              ggplot2::geom_sf() +
              ggplot2::theme_minimal() +
              ggplot2::ggtitle(paste("Test", i))
          })

          if (requireNamespace("gridExtra", quietly = TRUE)) {
            do.call(gridExtra::grid.arrange, c(plots, ncol = 2))
          } else {
            warning("Optional package 'gridExtra' not installed; printing plots sequentially.", call. = FALSE)
            for (p in plots) print(p)
          }
        }
      }
    }
  } else if (method == "cluster") {
    data_split <- spatial_clustering_cv(data = data_sf, v = fold, repeats = iter, ...)
    # ensure out_id present
    all_ids <- seq_len(nrow(data_sf))
    for (i in seq_along(data_split$splits)) {
      if (is.null(data_split$splits[[i]]$out_id) || all(is.na(data_split$splits[[i]]$out_id))) {
        in_id <- data_split$splits[[i]]$in_id
        data_split$splits[[i]]$out_id <- setdiff(all_ids, in_id)
      }
    }
    n_iter <- iter * fold
    if (plot_fold) print(autoplot(data_split))
  } else { # regularized
    data_split <- list(splits = vector("list", iter))
    for (i in seq_len(iter)) {
      locations_sf <- data_sf[!duplicated(sf::st_as_text(data_sf$geometry)), ]
      data_split$splits[[i]] <- list()
      data_split$splits[[i]]$data_test <- subsample.distance(locations_sf, size = n_size, d = min_dist * 1000, ...)
      test_geom <- sf::st_as_text(data_split$splits[[i]]$data_test$geometry)
      in_test   <- sf::st_as_text(data_sf$geometry) %in% test_geom
      data_split$splits[[i]]$out_id <- which(in_test)
      data_split$splits[[i]]$in_id  <- which(!in_test)
      data_split$splits[[i]]$data   <- data_sf[!in_test, ]
    }
    n_iter <- iter
    if (isTRUE(plot_fold)) {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("plot_fold = TRUE requires the 'ggplot2' package; skipping plots.", call. = FALSE)
      } else if (!requireNamespace("sf", quietly = TRUE)) {
        warning("plot_fold = TRUE with geom_sf() requires the 'sf' package; skipping plots.", call. = FALSE)
      } else {
        plots <- lapply(seq_len(n_iter), function(i) {
          ggplot2::ggplot(data_split$splits[[i]]$data_test) +
            ggplot2::geom_sf() +
            ggplot2::theme_minimal() +
            ggplot2::ggtitle(paste("Subset", i))
        })

        if (n_iter > 1 && requireNamespace("gridExtra", quietly = TRUE)) {
          do.call(gridExtra::grid.arrange, c(plots, ncol = 2))
        } else {
          # Either only one plot or gridExtra not available: print sequentially
          for (p in plots) print(p)
        }
      }
    }
  }

  ## ───────────────────────── initialise output ───────────────────────── ##
  n_models    <- length(object)
  model_names <- names(object)
  out <- list(test_set = vector("list", n_iter), model = list())

  ## ─────────────────────────── iterate over models ───────────────────────── ##
  for (h in seq_len(n_models)) {

    fit0      <- object[[h]]
    par_hat   <- coef(fit0)
    den_name  <- as.character(fit0$call$den)
    fam       <- fit0$family
    linkfun   <- switch(fam,
                        gaussian = identity,
                        binomial = plogis,
                        poisson  = exp)
    dast_flag <- is_dast_fit(fit0)
    dsgm_flag <- is_dsgm_fit(fit0)

    if (messages) {
      message(sprintf("\nModel '%s' (%s)", model_names[h], if (dast_flag) {
                                                             "DAST"
                                                           } else if(dsgm_flag) {
                                                             "DSGM"
                                                           } else {
                                                             "GLGM"
                                                           }
                                                          ))
    }

    ## containers for this model
    if (get_CRPS)   CRPS  <- vector("list", n_iter)
    if (get_SCRPS) { y_CRPS <- vector("list", n_iter); SCRPS <- vector("list", n_iter) }
    if (get_AnPIT) {
      if (fam == "gaussian") PIT <- vector("list", n_iter) else AnPIT <- vector("list", n_iter)
    }

    ## ─────────────────── iterate over CV splits (refit as needed) ─────────────────── ##
    for (i in seq_len(n_iter)) {

      in_id  <- data_split$splits[[i]]$in_id
      out_id <- data_split$splits[[i]]$out_id

      ## ----- refit or slice -----
      if (!keep_par_fixed) {
        message("\nRe-estimating model for subset ", i)
        crs_num <- if (is.numeric(fit0$crs)) {
          as.integer(fit0$crs)
        } else {
          as.integer(sub(".*:(\\d+)$", "\\1", fit0$crs))
        }

        fit0$crs <- crs_num
        if (!dast_flag & !dsgm_flag) {
          ## Original path (unchanged)
          refit_i <- eval(bquote(
            glgpm(.(
              formula      = fit0$formula,
              data         = data_sf[in_id, ],
              cov_offset   = .(fit0$cov_offset),
              family       = .(fam),
              crs          = .(fit0$crs),
              scale_to_km  = .(fit0$scale_to_km),
              control_mcmc = control_sim,
              fix_var_me   = .(fit0$fix_var_me),
              den          = .(as.name(den_name)),
              messages     = FALSE,
              start_pars   = par_hat
            ))
          ))
        } else if(dast_flag) {
          ## DAST refit
          time_sym <- fit0$call$time
          refit_i <- eval(bquote(
            dast(
              formula       = .(fit0$formula),
              data          = data_sf[in_id, ],
              den           = .(as.name(den_name)),
              time          = .(time_sym),
              mda_times     = .(fit0$mda_times),
              int_mat       = .(fit0$int_mat[in_id, , drop = FALSE]),
              penalty       = .(fit0$penalty),
              drop          = .(fit0$fix_alpha),
              power_val     = .(fit0$power_val),
              crs           = .(fit0$crs),
              scale_to_km   = .(fit0$scale_to_km),
              control_mcmc  = control_sim,
              messages      = FALSE
            )
          ))
        } else if (dsgm_flag) {
          ## DSGM refit
          time_sym <- fit0$call$time
          refit_i <- eval(bquote(
            dsgm(
              formula       = .(fit0$formula),
              data          = data_sf[in_id, ],
              den           = .(as.name(den_name)),
              time          = .(time_sym),
              intensity_family = .(fit0$intensity_family),
              gamma_sens = .(fit0$gamma_sens),
              mda_times     = .(fit0$mda_times),
              int_mat       = .(fit0$int_mat[in_id, , drop = FALSE]),
              penalty       = .(fit0$penalty),
              drop_W        = .(fit0$fix_alpha_W),
              decay_W       = .(fit0$fix_gamma_W),
              vary_k        = fit_dsgm$vary_k,
              power_val     = .(fit0$power_val),
              crs           = .(fit0$crs),
              scale_to_km   = .(fit0$scale_to_km),
              n_samples = fit0$stan_settings$n_samples,
              n_warmup = fit0$stan_settings$n_warmup,
              n_chains = fit0$stan_settings$n_chains,
              adapt_delta = fit0$stan_settings$adapt_delta,
              max_treedepth = fit0$stan_settings$max_treedepth,
              return_samples = FALSE,
              messages      = FALSE
            )
          ))
        }
      } else {
        ## quick slice without re-fitting
        refit_i <- fit0
        keep <- in_id
        refit_i$data_sf  <- refit_i$data_sf [keep, ]
        refit_i$units_m  <- refit_i$units_m[keep]
        keep_coord <- unique(refit_i$ID_coords[keep])
        refit_i$coords   <- refit_i$coords[keep_coord, , drop = FALSE]
        refit_i$y        <- refit_i$y      [keep]
        refit_i$D        <- refit_i$D      [keep, , drop = FALSE]
        if (!is.null(refit_i$cov_offset))
          refit_i$cov_offset <- refit_i$cov_offset[keep]
        if (!is.null(refit_i$ID_re)) {
          refit_i$ID_re <- refit_i$ID_re[keep, , drop = FALSE]
        }
        if (dast_flag | dsgm_flag) {
          refit_i$survey_times_data <- refit_i$survey_times_data[keep]
          refit_i$int_mat           <- refit_i$int_mat[keep, , drop = FALSE]
        }
        if(dsgm_flag) {
          if(refit_i$family=="intprev") {
            refit_i$prevalence_data <- refit_i$prevalence_data[keep]
            refit_i$egg_counts <- refit_i$egg_counts[keep]
            refit_i$intensity_data <- refit_i$egg_counts[refit_i$egg_counts>0]
          } else if(refit_i$family=="lf_mdiag") {
            refit_i$y_counts <- refit_i$y_counts[keep]
            refit_i$which_diag <- refit_i$which_diag[keep]
          }
        }
        ## recompute ID_coords mapping
        refit_i$ID_coords <- compute_ID_coords(refit_i$data_sf)$ID_coords
      }

      ## ----- held-out set and offsets -----
      data_test_i <- data_sf[out_id, ]
      #data_test_i <- data_test_i[complete.cases(sf::st_drop_geometry(data_test_i)), ]
      out$test_set[[i]] <- data_test_i

      pred_coff_i <- if (is.null(fit0$cov_offset)) rep(0, nrow(data_test_i)) else fit0$cov_offset[out_id]

      message("\nModel: ", model_names[h], "\nSpatial prediction for subset ", i)

      ## ----- prediction over test set -----
      pred_S <- pred_over_grid(
        object          = refit_i,
        grid_pred       = sf::st_as_sfc(data_test_i),
        control_sim     = control_sim,
        predictors      = data_test_i,
        pred_cov_offset = pred_coff_i,
        type            = "marginal",
        messages        = FALSE
      )
      if (dast_flag) {
        ## Build the DAST prediction inputs
        time_col  <- deparse(fit0$call$time)
        grid_pred_list <- list(
          geometry            = sf::st_as_sfc(data_test_i),
          survey_times_data   = data_test_i[[time_col]],
          int_mat             = fit0$int_mat[out_id, , drop = FALSE],
          mda_times           = fit0$mda_times
        )

        mda_eff_cov <- compute_mda_effect(data_test_i[[time_col]],
                                           mda_times = refit_i$mda_times,
                                           intervention = grid_pred_list$int_mat,
                                           alpha = coef.RiskMap(refit_i)$alpha,
                                           gamma = coef.RiskMap(refit_i)$gamma,
                                           kappa = refit_i$power_val)
        eta_samp <- t(pred_S$S_samples+pred_S$mu_pred)
        mu_samp  <- t((1/(1+exp(-eta_samp))))*mda_eff_cov
        eta_samp <- t(eta_samp)

      } else if(dsgm_flag) {
        ## Build the DSGM prediction inputs
        time_col  <- deparse(fit0$call$time)
        grid_pred_list <- list(
          geometry            = sf::st_as_sfc(data_test_i),
          survey_times_data   = data_test_i[[time_col]],
          int_mat             = fit0$int_mat[out_id, , drop = FALSE],
          mda_times           = fit0$mda_times
        )

        mda_eff_cov <- compute_mda_effect(data_test_i[[time_col]],
                                          mda_times = refit_i$mda_times,
                                          intervention = grid_pred_list$int_mat,
                                          alpha = coef.RiskMap(refit_i)$alpha_W,
                                          gamma = coef.RiskMap(refit_i)$gamma_W,
                                          kappa = 1)

        mu_W <- exp(pred_S$S_samples+pred_S$mu_pred)*mda_eff_cov
        rho_i <- coef(refit_i)$rho; k_i <- coef(refit_i)$k
        omega1_i <- coef(refit_i)$omega1

        if(refit_i$vary_k) {
          agg_W <- k_i*mu_W^omega1_i
        } else {
          agg_W <- k_i
        }

        prev_samples <- 1-(agg_W/(agg_W+mu_W*(1-exp(-rho_i))))^agg_W

        mu_C     = (rho_i * mu_W) / prev_samples
        sigma2_C =
          (rho_i * mu_W * (1 + rho_i)) / prev_samples +
          (rho_i^2 * mu_W^2 / prev_samples) * (1/agg_W + 1 - 1/prev_samples)

        if(refit_i$intensity_family=="negbin") {
          mu_C1    = mu_C - 1
          denom_nb = sigma2_C - mu_C1
          phi_C    = mu_C1^2 / denom_nb
        }

      } else {
        pred_lp  <- pred_target_grid(
          pred_S,
          include_nugget     = is.null(refit_i$fix_tau2) || refit_i$fix_tau2 != 0,
          include_cov_offset = !all(refit_i$cov_offset == 0)
        )
        eta_samp <- pred_lp$lp_samples
        if (fam == "gaussian") {
          sigma2_me <- if (is.null(refit_i$fix_var_me)) coef(refit_i)$sigma2_me else refit_i$fix_var_me
          eta_samp <- eta_samp + sqrt(sigma2_me) * stats::rnorm(length(eta_samp))
        }
        mu_samp <- linkfun(eta_samp)
      }

      if(!dsgm_flag) {
        n_pred <- nrow(eta_samp)
        n_draw <- ncol(eta_samp)
      } else {
        n_pred <- nrow(mu_C1)
        n_draw <- ncol(mu_C1)
      }

      if (get_CRPS)   CRPS [[i]] <- numeric(n_pred)
      if (get_SCRPS) { y_CRPS[[i]] <- numeric(n_pred); SCRPS[[i]] <- numeric(n_pred) }

      if (get_AnPIT) {
        if (fam == "gaussian") {
          PIT_i <- numeric(n_pred)
        } else {
          u_val   <- seq(0, 1, length.out = 1000)
          AnPIT_i <- matrix(NA_real_, nrow = length(u_val), ncol = n_pred)
          npit_fun <- function(y, u, Fk) {
            F_y1 <- if (y == 0) 0 else Fk[y]
            F_y  <- Fk[y + 1]
            ifelse(u <= F_y1, 0,
                   ifelse(u <= F_y, (u - F_y1) / (F_y - F_y1), 1))
          }
        }
      }

      if(dsgm_flag) {
        y_i <- fit0$egg_counts[out_id]
      } else {
        units_m_i <- fit0$units_m[out_id]
        y_i       <- fit0$y[out_id]
      }

      for (j in seq_len(n_pred)) {
        if (fam == "gaussian") {
          mu_j <- mean(mu_samp[j, ])
          sd_j <- stats::sd(mu_samp[j, ])
          if (get_CRPS)  CRPS[[i]][j] <- sd_j * crps_gaussian(y_i[j], mu_j, sd_j)
          if (get_SCRPS) {
            y_CRPS[[i]][j] <- sd_j / sqrt(pi)
            SCRPS [[i]][j] <- -0.5 * (1 + CRPS[[i]][j] / y_CRPS[[i]][j] + log(2 * abs(y_CRPS[[i]][j])))
          }
          if (get_AnPIT) PIT_i[j] <- stats::pnorm(y_i[j], mean = mu_j, sd = sd_j)
        } else {
          if (fam == "binomial") {
            y_samp  <- stats::rbinom(n_draw, size = units_m_i[j], prob = mu_samp[j, ])
            support <- 0:units_m_i[j]
          } else if(fam == "poisson") { # Poisson
            lambda  <- units_m_i[j] * mu_samp[j, ]
            y_samp  <- stats::rpois(n_draw, lambda)
            support <- 0:max(max(y_samp), y_i[j], stats::qpois(0.999, mean(lambda)))
          } else if (fam == "intprev") {
            if (refit_i$intensity_family == "negbin") {
              y_samp <- integer(n_draw)
              y_ind  <- rbinom(n_draw, size = 1, prob = prev_samples[j, ])
              pos    <- which(y_ind == 1)
              if (length(pos) > 0) {
                y_samp[pos] <- MASS::rnegbin(
                  length(pos),
                  mu    = mu_C1[j, pos],
                  theta = phi_C[j, pos]
                ) + 1
              }
              # y_samp[y_ind == 0] stays 0  -> structural zeros retained
              support <- 0:max(max(y_samp), y_i[j],
                               qnbinom(0.999,
                                       size = mean(phi_C[j, ]),
                                       mu   = mean(mu_C1[j, ])))
            }
          }
          pk <- tabulate(y_samp + 1, nbins = length(support)) / n_draw
          Fk <- cumsum(pk)
          if (get_CRPS)  CRPS[[i]][j] <- crps_discrete(y_i[j], Fk)
          if (get_SCRPS) {
            y_CRPS[[i]][j] <- exp_crps_discrete(Fk, pk)
            SCRPS [[i]][j] <- -0.5 * (1 + CRPS[[i]][j] / y_CRPS[[i]][j] + log(2 * abs(y_CRPS[[i]][j])))
          }
          if (get_AnPIT) {
            AnPIT_i[, j] <- vapply(u_val, npit_fun, numeric(1), y = y_i[j], Fk = Fk)
          }
        }
      }

      if (get_AnPIT) {
        if (fam == "gaussian") PIT[[i]] <- PIT_i else AnPIT[[i]] <- rowMeans(AnPIT_i)
      }

    } # end i loop

    ## ─────────────── finalise output for this model ─────────────── ##
    out$model[[model_names[h]]] <- list(score = list())
    if (get_CRPS)  out$model[[model_names[h]]]$score$CRPS  <- CRPS
    if (get_SCRPS) out$model[[model_names[h]]]$score$SCRPS <- SCRPS
    if (get_AnPIT) {
      if (fam == "gaussian") out$model[[model_names[h]]]$PIT <- PIT else out$model[[model_names[h]]]$AnPIT <- AnPIT
    }

  } # end h loop

  class(out) <- "RiskMap.spatial.cv"
  return(out)
}



##' Simulate surface data based on a spatial model
##'
##' This function simulates surface data based on a user-defined formula and other parameters. It allows for simulation of spatial data with various model families (Gaussian, Binomial, or Poisson). The simulation involves creating spatially correlated random fields and generating outcomes for data points in a given prediction grid.
##'
##' @param n_sim The number of simulations to run.
##' @param pred_grid A spatial object (either `sf` or `data.frame`) representing the prediction grid where the simulation will take place.
##' @param formula A formula object specifying the model to be fitted. It should include both fixed effects and random effects if applicable.
##' @param sampling_f A function that returns a sampled dataset (of class `sf` or `data.frame`) to simulate data from.
##' @param family A character string specifying the family of the model. Must be one of "gaussian", "binomial", or "poisson".
##' @param scale_to_km A logical indicating whether the coordinates should be scaled to kilometers. Defaults to `TRUE`.
##' @param control_mcmc A list of control parameters for MCMC (not used in this implementation but can be expanded later).
##' @param par0 A list containing initial parameter values for the simulation, including `beta`, `sigma2`, `phi`, `tau2`, and `sigma2_me`.
##' @param include_covariates A logical indicateing if the covariates (or the intercept if no covariates are used) should be included in the linear
##' predictor. By default \code{include_covariates = TRUE}
##' @param nugget_over_grid A logical indicating whether to include a nugget effect over the entire prediction grid.
##' @param fix_var_me A parameter to fix the variance of the random effects for the measurement error. Defaults to `NULL`.
##' @param messages A logical value indicating whether to print messages during the simulation. Defaults to `TRUE`.
##'
##' @return A list containing the simulated data (\code{data_sim}), the linear predictors (\code{lp_grid_sim}),
##' a logical value indicating if covariates have been included in the linear predictor (\code{include_covariates}),
##' a logical value indicating if the nugget has been included into the simulations of the linear predictor over the grid
##' (\code{nugget_over_grid}), a logical  indicating if a covariate offset has been included in the linear predictor (\code{include_cov_offset}),
##' the model parameters set for the simulation (\code{par0}) and the family used in the model (\code{family}).
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @export
surf_sim <- function(n_sim, pred_grid,
                     formula, sampling_f,
                     family,
                     scale_to_km = TRUE,
                     control_mcmc = set_control_sim(),
                     par0, nugget_over_grid = FALSE,
                     include_covariates = TRUE,
                     fix_var_me = NULL,
                     messages = TRUE) {


  if(!inherits(formula,
               what = "formula", which = FALSE)) {
    stop("'formula' must be a 'formula'
                                     object indicating the variables of the
                                     model to be fitted")
  }
  inter_f <- interpret.formula(formula)
  include_cov_offset <- !is.null(inter_f$offset)
  if(!inherits(pred_grid,
               what = c("sf", "data.frame"), which = FALSE)) {
    stop("'pred_grid' must be an 'sf'
          object indicating the variables of the
          model to be fitted")
  }

  if(!inherits(sampling_f,
               what = c("function"), which = FALSE)) {
    stop("'sampling_f' must be an object of class 'function'")
  }

  data_test <- sampling_f()
  if(!inherits(data_test,
               what = c("sf", "data.frame"), which = FALSE)) {
    stop("The object return by 'sampling_f' must be of an 'sf' object")
  }

  if (!"units_m" %in% colnames(data_test)) {
    stop("The object returned by 'sampling_f' must contain a column named 'units_m'.")
  }


  data_sim <- list()
  coords_sim <- list()
  for(i in 1:n_sim) {
    data_sim[[i]] <- sampling_f()
    coords_sim[[i]] <- st_coordinates(data_sim[[i]])
    if(scale_to_km) coords_sim[[i]] <- coords_sim[[i]]/1000
    if (st_crs(data_sim[[i]]) != st_crs(pred_grid)) {
      pred_grid <- st_transform(pred_grid, st_crs(data_sim[[i]]))
    }

    # Find nearest neighbor in 'pred_grid' for each feature in 'data'
    nearest_indices <- st_nearest_feature(data_sim[[i]], pred_grid)

    # Extract variables from the nearest features in 'pred_grid'
    pred_grid_vars <- st_drop_geometry(pred_grid[nearest_indices, ])

    # Bind the extracted variables to 'data'
    data_sim[[i]] <- cbind(data_sim[[i]], pred_grid_vars)
  }


  kappa <- inter_f$gp.spec$kappa
  if(kappa < 0) stop("kappa must be positive.")

  if(family != "gaussian" & family != "binomial" &
     family != "poisson") stop("'family' must be either 'gaussian', 'binomial'
                               or 'poisson'")

  inter_lt_f <- inter_f
  inter_lt_f$pf <- update(inter_lt_f$pf, NULL ~.)
  D <- list()
  n <- list()
  for(i in 1:n_sim) {
    mf <- model.frame(inter_lt_f$pf,data=data_sim[[i]], na.action = na.fail)


    # Extract covariates matrix
    D[[i]] <- as.matrix(model.matrix(attr(mf,"terms"),data=data_sim[[i]]))
    n[[i]] <- nrow(D[[i]])
  }


  if(length(inter_f$re.spec) > 0) {
    stop("In the current impletementation of 'surf_sim' the addition of random effects
        with re() is not supported")
  }
  mf_pred <- model.frame(inter_lt_f$pf,data=pred_grid, na.action = na.fail)
  D_pred <- as.matrix(model.matrix(attr(mf_pred,"terms"),data=pred_grid))
  if(ncol(D_pred)!=ncol(D[[1]])) stop("the provided variables in 'grid_pred' do not match the number of explanatory variables used to fit the model.")

  # Number of covariates
  p <- ncol(D[[1]])

  beta <- par0$beta

  if(length(beta)!=p) stop("The values passed to 'beta' do not match the variables in 'grid_pred'")

  sigma2 <- par0$sigma2
  phi <- par0$phi

  if(is.null(par0$tau2)) {
    tau2 <- 0
  } else {
    tau2 <- par0$tau2
  }

  if(is.null(fix_var_me)) {
    sigma2_me <- 0
  } else {
    sigma2_me <- par0$sigma2_me
  }

  ID_coords <- sapply(1:n_sim, function(i) 1:n[[i]])

  grid_pred <- st_coordinates(pred_grid)
  if(scale_to_km) grid_pred <- grid_pred/1000
  n_pred <- nrow(grid_pred)

  # Simulate on the grid
  # Simulate S
  Sigma <- sigma2*matern_cor(dist(grid_pred), phi = phi, kappa = kappa,
                             return_sym_matrix = TRUE)
  if(nugget_over_grid) {
    diag(Sigma) <- diag(Sigma) + tau2
  }
  Sigma_sroot <- t(chol(Sigma))

  S_sim <- sapply(1:n_sim, function(i) Sigma_sroot%*%rnorm(n_pred))

  sim_columns <- paste0("sim_", 1:n_sim)
  S_sim_df <- as.data.frame(S_sim)
  colnames(S_sim_df) <- sim_columns

  # Combine simulations with grid_pred
  S_sim <- cbind(grid_pred, S_sim_df)
  S_sim <- st_sf(S_sim, geometry = st_geometry(pred_grid))

  coords_tot <- list()
  out <- list()

  for(i in 1:n_sim) {
    coords_tot[[i]] <- rbind(coords_sim[[i]], grid_pred)
    ind_data <- 1:n[[i]]
    ind_pred <- (n[[i]]+1):(n[[i]]+n_pred)

    D_tot <- rbind(D[[i]], D_pred)

    # Find the nearest features in grid_pred_with_sim for each location in data_sim[[i]]
    nearest_indices <- st_nearest_feature(data_sim[[i]],
                                          S_sim)

    # Extract the i-th simulation values
    sim_column <- paste0("sim_", i)  # Column name for the i-th simulation
    S_sim_data <- S_sim[nearest_indices, sim_column, drop = TRUE]

    if(tau2>0 & !nugget_over_grid) {
      S_sim_data <- S_sim_data + sqrt(tau2)*rnorm(n[[i]])
    }
    # Linear predictor
    if(include_covariates) {
      eta_sim_tot <- D_tot%*%beta +
        c(S_sim_data, S_sim[[sim_column]])
    } else {
      eta_sim_tot <- c(S_sim_data, S_sim[[sim_column]])
    }

    data_sim[[i]]$y <- NA

    if(family=="gaussian") {

      data_sim[[i]]$y <- eta_sim_tot[1:n[[i]]] + sqrt(sigma2_me)*rnorm(n)

    } else if(family=="binomial") {
      prob_i <- exp(eta_sim_tot[1:n[[i]]])/(1+exp(eta_sim_tot[1:n[[i]]]))
      data_sim[[i]]$y <- rbinom(n[[i]], size = data_sim[[i]]$units_m,
                                prob = prob_i)
    } else if(family=="poisson") {
      mean_i <- data_sim[[i]]$units_m*exp(eta_sim_tot[1:n[[i]]])
      data_sim[[i]]$y <- rpois(n[[i]], lambda = mean_i)
    }

    data_sim[[i]]$lp_data <- S_sim_data

    pred_grid[[paste0("lp_",sim_column)]] <- eta_sim_tot[-(1:n[[i]])]

  }

  out$data_sim <- data_sim
  out$lp_grid_sim <- pred_grid
  out$include_covariates <- include_covariates
  out$nugget_over_grid <- nugget_over_grid
  out$include_cov_offset <- include_cov_offset
  out$par0 <- par0
  out$family <- family
  class(out) <- "RiskMap.sim"
  return(out)
}

##' Plot simulated surface data for a given simulation
##'
##' This function plots the simulated surface data for a specific simulation from the result of `surf_sim`. It visualizes the linear predictor values on a raster grid along with the actual data points.
##'
##' @param surf_obj The output object from `surf_sim`, containing both simulated data (`data_sim`) and predicted grid simulations (`lp_grid_sim`).
##' @param sim The simulation index to plot.
##' @param ... Additional graphical parameters to be passed to the plotting function of the `terra` package.
##'
##' @return A plot of the simulation results.
##'
##' @importFrom stars st_rasterize
##'
##' @export
plot_sim_surf <-  function(surf_obj, sim, ...) {

  sf_object <- surf_obj$lp_grid_sim
  value_column <- paste0("lp_sim_",sim)
  r <- rast(st_rasterize(sf_object[,c("x",value_column)]))
  r[r == 0] <- NA

  plot(r, main = paste("Simulation no.", sim), ...)
  points(st_coordinates(surf_obj$data_sim[[sim]]), pch = 20)

}

##' @title Assess Simulations
##'
##' @description This function evaluates the performance of models based on simulation results from the `surf_sim` function.
##'
##' @param obj_sim An object of class `RiskMap.sim`, obtained as an output from the `surf_sim` function.
##' @param models A named list of models to be evaluated.
##' @param control_mcmc A control object for MCMC sampling, created with `set_control_sim()`. Default is `set_control_sim()`.
##' @param spatial_scale The scale at which predictions are assessed, either `"grid"` or `"area"`.
##' @param messages Logical, if `TRUE` messages will be displayed during processing. Default is `TRUE`.
##' @param f_grid_target A function for processing grid-level predictions.
##' @param f_area_target A function for processing area-level predictions.
##' @param shp A shapefile of class `sf` or `data.frame` for area-level analysis, required if `spatial_scale = "area"`.
##' @param col_names Column name in `shp` containing unique region names. If `NULL`, defaults to `"region"`.
##' @param pred_objective A character vector specifying objectives, either `"mse"`, `"classify"`, or both.
##' @param categories A numeric vector of thresholds defining categories for classification. Required if `pred_objective = "classify"`.
##'
##' @return A list of class `RiskMap.sim.res` containing model evaluation results.
##'
##' @export
assess_sim <- function(obj_sim,
                       models,
                       control_mcmc = set_control_sim(),
                       spatial_scale,
                       messages = TRUE,
                       f_grid_target = NULL,
                       f_area_target = NULL,
                       shp = NULL, col_names = NULL,
                       pred_objective = c("mse","classify"),
                       categories= NULL) {

  if (!inherits(obj_sim, "RiskMap.sim")) {
    stop("'obj_sim' must be an object of class 'RiskMap.sim' obtained as an output from the 'surf_sim' function")
  }
  if (length(setdiff(pred_objective, c("mse","classify")))>0) {
    stop(paste("Invalid value for pred_objective. Allowed values are:", paste(c("mse","classify"), collapse = ", ")))
  }
  if(spatial_scale != "grid" & spatial_scale != "area") {
    stop("'spatial_scale' must be set to 'grid' or 'area'")
  }
  if(spatial_scale=="area" & is.null(shp)) {
    stop("if spatial_scale='area' then a shape file of the area(s) must be passed to
         'shp'")
  }
  units_m <- NULL
  if(any(pred_objective=="classify")) {
    if(is.null(categories)) stop("if pred_objective='class', a value for 'categories' must be specified")
    if (length(categories) < 3) {
      stop("Categories vector must contain at least three unique, strictly increasing values.")
    }
  }
  n_sim <- length(obj_sim$data_sim)
  n_models <- length(models)

  if(spatial_scale == "area" & is.null(f_area_target)) {
    stop("If 'spatial_scale' is set to 'area', then 'f_area_target' must be provided")
  }
  model_names <- names(models)

  include_covariates <- obj_sim$include_covariates
  include_cov_offset <- obj_sim$include_cov_offset
  include_nugget <- obj_sim$nugget_over_grid

  fits <- list()
  preds <- list()
  if(spatial_scale=="grid") {
    type <- "marginal"
  } else if(spatial_scale=="area") {
    type <- "joint"
    n_reg <- nrow(shp)
    if(is.null(shp)) stop("If spatial_scale='area', then 'shp' must be specified")
    if(!inherits(shp,
                 what = c("sf","data.frame"), which = FALSE)) {
      stop("The object passed to 'shp' must be an object of class 'sf'")
    }

    if(is.null(col_names)) {
      shp$region <- paste("reg",1:n_reg, sep="")
      col_names <- "region"
      names_reg <- shp$region
    } else {
      names_reg <- shp[[col_names]]
      if(n_reg != length(names_reg)) {
        stop("The names in the column identified by 'col_names' do not
         provide a unique set of names, but there are duplicates")
      }
    }
    shp <- st_transform(shp, st_crs(obj_sim$lp_grid_sim))
    inter <- st_intersects(shp, obj_sim$lp_grid_sim)
  }

  for(i in 1:n_models) {
    if(messages) message("Model: ", paste(model_names[i]),"\n")

    if_i <- interpret.formula(models[[i]])
    rhs_terms <- attr(terms(if_i$pf), "term.labels")
    # Check if there are any covariates
    if (length(rhs_terms) == 0) {
      predictors_i <- NULL
    } else {
      predictors_i <- obj_sim$lp_grid_sim
    }
    for(j in 1:n_sim) {
      if(messages) message("Processing simulation no.", j)
      f_i <- update(models[[i]], y ~ .)
      if(messages) message("Estimation")
      fits[[paste(model_names[i])]][[j]] <- glgpm(formula = f_i,
                                                  den = units_m,
                                                  family = obj_sim$family,
                                                  data = obj_sim$data_sim[[j]],
                                                  control_mcmc = control_mcmc,
                                                  messages = FALSE)

      if(messages) message("Prediction over the grid")
      preds[[paste(model_names[i])]][[j]] <-
        pred_over_grid(fits[[paste(model_names[i])]][[j]],
                       grid_pred = st_as_sfc(obj_sim$lp_grid_sim),
                       predictors = predictors_i,
                       type = type, messages = FALSE)
    }
  }

  n_samples <- (control_mcmc$n_sim-control_mcmc$burnin)/control_mcmc$thin
  n_pred <- nrow(obj_sim$lp_grid_sim)


  out <- list(pred_objective = list())

  if(any(pred_objective=="mse")) {
    out$pred_objective$mse <- array(NA, c(n_models, n_sim))
    rownames(out$pred_objective$mse) <- model_names
    colnames(out$pred_objective$mse) <- paste0("sim_",1:n_sim)
  }

  if(any(pred_objective=="classify")) {

    # Ensure categories are unique and strictly increasing
    categories <- unique(sort(categories))


    # Assign classification to the output object
    out$pred_objective$classify <- setNames(vector("list", length(model_names)), model_names)

    # Correctly generate breaks and labels
    breaks <- categories  # Use categories directly as breaks
    categories_class <- factor(paste0("(", head(categories, -1), ",",
                                      categories[-1], "]"))  # Labels to match intervals



    for(i in 1:n_models) {
      out$pred_objective$classify[[model_names[i]]] <- list(by_cat = list(),
                                                            across_cat = list())
      out$pred_objective$classify[[model_names[i]]]$by_cat <- vector("list", n_sim)
      for(j in 1:n_sim) {
        out$pred_objective$classify[[paste(model_names[i])]]$by_cat[[j]] <-
          data.frame(
            Class = categories_class,
            Sensitivity = NA,
            Specificity = NA,
            PPV = NA,
            NPV = NA,
            CC = NA
          )
      }
      out$pred_objective$classify[[model_names[i]]]$CC <- rep(NA,n_sim)
    }
  }
  lp_true_sim <- st_drop_geometry(obj_sim$lp_grid_sim[, grepl("lp_sim",
                                                              names(obj_sim$lp_grid_sim))])

  if(spatial_scale == "grid") {
    true_target_sim <- f_grid_target(lp_true_sim)
  } else if(spatial_scale == "area") {
    true_target_sim <- matrix(NA, nrow = n_reg, ncol = n_sim)
    true_target_grid_sim <- f_grid_target(lp_true_sim)
    for(i in 1:n_reg) {
      for(j in 1:n_sim) {
        if(length(inter[[i]])==0) {
          warning(paste("No points on the grid fall within", shp[[col_names]][h],
                        "and no predictions are carried out for this area"))
          no_comp <- c(no_comp, h)
        } else {
          true_target_sim[i,j] <- f_area_target(true_target_grid_sim[inter[[i]],j])
        }
      }
    }
  }



  for(i in 1:n_models) {
    for(j in 1:n_sim) {
      obj_pred_ij <- preds[[paste(model_names[i])]][[j]]
      if(length(obj_pred_ij$mu_pred)==1 && obj_pred_ij$mu_pred==0 &&
         include_covariates) {
        stop("Covariates have not been provided; re-run pred_over_grid
         and provide the covariates through the argument 'predictors'")
      }


      if(!include_covariates) {
        mu_target <- 0
      } else {

        if(is.null(obj_pred_ij$mu_pred)) stop("the output obtained from 'pred_S' does not
                                     contain any covariates; if including covariates
                                     in the predictive target these shuold be included
                                     when running 'pred_S'")
        mu_target <- obj_pred_ij$mu_pred
      }

      if(!include_cov_offset) {
        cov_offset <- 0
      } else {
        if(length(obj_pred_ij$cov_offset)==1) {
          stop("No covariate offset was included in the model;
           set include_cov_offset = FALSE, or refit the model and include
           the covariate offset")
        }
        cov_offset <- obj_pred_ij$cov_offset
      }

      if(include_nugget) {
        if(is.null(obj_pred_ij$par_hat$tau2)) stop("'include_nugget' cannot be
                                                   set to TRUE if this has not been included
                                                   in the fit of the model")
        Z_sim <- matrix(rnorm(n_samples*n_pred,
                              sd = sqrt(obj_pred_ij$par_hat$tau2)),
                        ncol = n_samples)
        obj_pred_ij$S_samples <- obj_pred_ij$S_samples+Z_sim
      }

      if(is.matrix(mu_target)) {
        lp_samples_ij <- sapply(1:n_samples,
                                function(h)
                                  mu_target[,h] + cov_offset +
                                  obj_pred_ij$S_samples[,h])
      } else {
        lp_samples_ij <- sapply(1:n_samples,
                                function(h)
                                  mu_target + cov_offset +
                                  obj_pred_ij$S_samples[,h])
      }

      target_samples_ij <- f_grid_target(lp_samples_ij)

      if(spatial_scale == "grid") {
        mean_target_ij <- apply(target_samples_ij, 1, mean)
      } else if(spatial_scale == "area") {
        target_area_samples_ij <- matrix(NA, nrow = n_reg, ncol = n_samples)
        mean_target_ij <- rep(NA,n_reg)
        for(h in 1:n_reg) {
          if(length(inter[[h]])==0) {
            warning(paste("No points on the grid fall within", shp[[col_names]][h],
                          "and no predictions are carried out for this area"))
            no_comp <- c(no_comp, h)
          } else {
            ind_grid_h <- inter[[h]]
            target_area_samples_ij[h,] <-  apply(target_samples_ij[ind_grid_h,], 2,
                                                 f_area_target)
            mean_target_ij[h] <- mean(target_area_samples_ij[h,])
          }
        }
      }

      if(any(pred_objective=="mse")) {
        out$pred_objective$mse[i,j] <- mean((mean_target_ij-true_target_sim[,j])^2)
      }

      if(any(pred_objective=="classify")) {
        true_class_ij <- cut(true_target_sim[,j], breaks = categories)
        n_categories <- length(categories)-1
        if(spatial_scale == "grid") {
          prob_cat_ij <- matrix(0, nrow=n_pred, ncol = n_categories)
        } else if(spatial_scale == "area") {
          prob_cat_ij <- matrix(0, nrow=n_reg, ncol = n_categories)
        }
        for(h in 1:(n_categories)) {
          if(spatial_scale == "grid") {
            prob_cat_ij[,h] <- apply(categories[h] < target_samples_ij &
                                       categories[h+1] > target_samples_ij, 1, mean)
          } else if(spatial_scale == "area") {
            prob_cat_ij[,h] <- apply(categories[h] < target_area_samples_ij &
                                       categories[h+1] > target_area_samples_ij, 1, mean)
          }
        }

        pred_class_ij <- apply(prob_cat_ij, 1, function(x) categories_class[which.max(x)])

        # Define the confusion matrix
        conf_matrix <- table(true_class_ij, pred_class_ij)

        # Calculate metrics for each class
        for (h in 1:nrow(conf_matrix)) {
          TP <- conf_matrix[h, h]  # True Positives: diagonal entry
          FP <- sum(conf_matrix[, h]) - conf_matrix[h, h]  # False Positives: column sum minus diagonal
          FN <- sum(conf_matrix[h, ]) - conf_matrix[h, h]  # False Negatives: row sum minus diagonal
          TN <- sum(conf_matrix) - sum(conf_matrix[h, ]) - sum(conf_matrix[, h]) + conf_matrix[h, h]


          # Handle cases where division by zero could occur
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$Sensitivity <-
            ifelse((TP + FN) == 0, NA, TP / (TP + FN))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$Specificity <-
            ifelse((TN + FP) == 0, NA, TN / (TN + FP))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$PPV <-
            ifelse((TP + FP) == 0, NA, TP / (TP + FP))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$NPV <-
            ifelse((TN + FN) == 0, NA, TN / (TN + FN))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$CC <-
            ifelse((TP + FN) == 0, NA, (TP ) / (TP + FN))
        }
        out$pred_objective$classify[[paste(model_names[[i]])]]$CC[j] <-
          mean(true_class_ij==pred_class_ij)
      }
    }
  }
  if(any(pred_objective=="classify")) out$pred_objective$classify$Class <- categories_class
  class(out) <- "RiskMap.sim.res"
  return(out)
}

##' @title Summarize Simulation Results
##'
##' @description Summarizes the results of model evaluations from a `RiskMap.sim.res` object. Provides average metrics for classification by category and overall correct classification (CC) summary.
##'
##' @param object An object of class `RiskMap.sim.res`, as returned by `assess_sim`.
##' @param ... Additional arguments (not used).
##'
##' @return A list containing summary data for each model:
##' - `by_cat_summary`: A data frame with average sensitivity, specificity, PPV, NPV, and CC by category.
##' - `CC_summary`: A numeric vector with mean, 2.5th percentile, and 97.5th percentile for CC across simulations.
##'
##' @method summary RiskMap.sim.res
##' @export
summary.RiskMap.sim.res <- function(object, ...) {
  stopifnot(inherits(object, "RiskMap.sim.res"))

  # Initialize results
  results <- list()

  # Check for "mse" in pred_objective
  if ("mse" %in% names(object$pred_objective)) {
    mse_data <- object$pred_objective$mse

    # Check if mse_data is a matrix
    if (is.matrix(mse_data)) {
      # Compute mean and SD for each model (row)
      mse_summary <- data.frame(
        Model = rownames(mse_data),
        MSE_mean = rowMeans(mse_data, na.rm = TRUE),
        MSE_sd = apply(mse_data, 1, sd, na.rm = TRUE)
      )

      results$mse <- mse_summary
    } else {
      stop("mse_data must be a matrix.")
    }
  }

  # Check for "classify" in pred_objective
  if ("classify" %in% names(object$pred_objective)) {
    classify_data <- object$pred_objective$classify

    # Loop over each model (e.g., M1, M2)
    n_models <- length(classify_data)-1
    name_models <- names(classify_data)[1:n_models]
    results$classify <- list()
    for(i in 1:n_models) {

      model_data <- classify_data[[i]]

      n_sim <- length(model_data$by_cat)
      res_class <- model_data$by_cat[[1]][,-1]
      den <- 0
      for(j in 2:n_sim) {
        if(!any(is.na(model_data$by_cat[[j]][,-1]))) {
          den <- den + 1
          res_class <- res_class+model_data$by_cat[[j]][,-1]
        }
      }
      res_class <- data.frame(res_class/den)
      res_class$Class <- model_data$by_cat[[1]][,1]

      cc_summary <- list(mean = mean(model_data$CC, na.rm = TRUE),
                         lower = quantile(model_data$CC, 0.025, na.rm = TRUE),
                         upper = quantile(model_data$CC, 0.975, na.rm = TRUE))

      results$classify[[paste(name_models[i])]] <- list(classify_res = res_class,
                                            cc_summary = list(mean = mean(model_data$CC, na.rm = TRUE),
                                                                     lower = quantile(model_data$CC, 0.025, na.rm = TRUE),
                                                                     upper = quantile(model_data$CC, 0.975, na.rm = TRUE)))
    }
  }

  # Assign class for S3 print method
  class(results) <- "summary.RiskMap.sim.res"
  return(results)
}



##' @title Print Simulation Results
##'
##' @description Prints a concise summary of simulation results from a `RiskMap.sim.res` object, including average metrics by category and a summary of overall correct classification (CC).
##'
##' @param x An object of class `summary.RiskMap.sim.res`, as returned by `summary.RiskMap.sim.res`.
##' @param ... Additional arguments (not used).
##'
##' @return Invisibly returns `x`.
##'
##'
##' Print Simulation Results
##'
##' Prints a concise summary of simulation results from a `summary.RiskMap.sim.res` object,
##' including average metrics by category and a summary of overall correct classification (CC).
##'
##' @param x An object of class `summary.RiskMap.sim.res`, as returned by `summary.RiskMap.sim.res`.
##' @param ... Additional arguments (not used).
##'
##' @return Invisibly returns `x`.
##'
##' @method print summary.RiskMap.sim.res
##' @export
print.summary.RiskMap.sim.res <- function(x, ...) {
  cat("Summary of Simulation Results\n\n")

  if (!is.null(x$mse)) {
    cat("Mean Squared Error (MSE):\n")
    print(x$mse)
    cat("\n")
  }

  if (!is.null(x$classify)) {
    cat("Classification Results:\n")

    # Iterate over each model in classify results
    for (model_name in names(x$classify)) {
      model_data <- x$classify[[model_name]]

      cat(sprintf("\nModel: %s\n", model_name))

      cat("\nAverages across simulations by Category:\n")
      print(model_data$classify_res)

      cat("\nProportion of Correct Classification (CC) across categories:\n")
      cc_summary <- model_data$cc_summary
      cat(sprintf("Mean: %.3f, 95%% CI: [%.3f, %.3f]\n",
                  cc_summary$mean, cc_summary$lower, cc_summary$upper))
    }
    cat("\n")
  }

  invisible(x)
}
