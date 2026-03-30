compute_mda_effect <- function(survey_times_data, mda_times, intervention,
                               alpha, gamma, kappa) {
  n <- length(survey_times_data)
  effect <- rep(NA, n)

  mda_effect_f <- function(v, alpha, gamma, kappa) {
    alpha*exp(-(v/gamma)^kappa)
  }
  mda_effect_f <- Vectorize(mda_effect_f, "v")

  f <- function(t, mda_times, int, alpha, gamma, kappa) {
    ind_t <- which(t > mda_times)
    u_j <- mda_times[ind_t]
    if(length(u_j) > 0) {
      out <- prod((1-mda_effect_f(t-u_j, alpha, gamma, kappa))^int[ind_t])
    } else {
      out <- 1
    }
    return(out)
  }

  f <- Vectorize(f, "t")

  for(i in 1:n) {
    effect[i] <- f(survey_times_data[i], mda_times, intervention[i,],
                   alpha, gamma, kappa)
  }
  return(effect)
}


compute_mda_effect_derivatives <- function(survey_times_data, mda_times, intervention,
                                           alpha, gamma, kappa) {
  n <- length(survey_times_data)
  effect <- numeric(n)
  d_alpha <- numeric(n)
  d_gamma <- numeric(n)
  d2_alpha <- numeric(n)
  d2_gamma <- numeric(n)
  d2_alpha_gamma <- numeric(n)

  for (i in 1:n) {
    t <- survey_times_data[i]
    int <- intervention[i,]

    ind_t <- which(t > mda_times)
    u_j <- mda_times[ind_t]
    int_j <- int[ind_t]

    if (length(u_j) > 0) {
      v_j <- t - u_j
      z_j <- (v_j/gamma)^kappa
      e_j <- exp(-z_j)
      denom <- 1 - alpha * e_j
      effect[i] <- exp(sum(int_j * log(denom)))

      # First derivatives (CORRECT - matches numeric)
      term_alpha <- int_j * e_j / denom
      dlogP_dalpha <- -sum(term_alpha)

      gamma_term <- kappa * z_j / gamma
      term_gamma <- alpha * term_alpha * gamma_term
      dlogP_dgamma <- -sum(term_gamma)

      # Second derivatives
      # ∂²(log P)/∂α² (CORRECT - matches numeric)
      d2logP_dalpha2 <- -sum(int_j * e_j^2 / denom^2)

      #######################################################
      # CORRECTED ∂²(log P)/∂γ² (FIXED SIGN)
      # f = αe^{-z}, z = (v/γ)^κ
      # ∂f/∂γ = f * (κz/γ)
      # ∂²f/∂γ² = f * [κz(κ+1 - κz)/γ²]
      term1 <- (alpha * e_j * kappa * z_j * (kappa + 1 - kappa*z_j)) /
        (gamma^2 * denom)
      term2 <- (alpha * e_j * kappa * z_j / (gamma * denom))^2
      d2logP_dgamma2 <- sum(int_j * (term1 - term2))
      #######################################################

      # ∂²(log P)/∂α∂γ (CORRECT - matches numeric)
      term_cross <- e_j * gamma_term / denom
      d2logP_dalphadgamma <- sum(int_j * (-term_cross - alpha * e_j * term_cross / denom))

      # Convert to derivatives of P
      d_alpha[i] <- effect[i] * dlogP_dalpha
      d_gamma[i] <- effect[i] * dlogP_dgamma
      d2_alpha[i] <- effect[i] * (dlogP_dalpha^2 + d2logP_dalpha2)
      d2_gamma[i] <- effect[i] * (dlogP_dgamma^2 + d2logP_dgamma2)
      d2_alpha_gamma[i] <- effect[i] * (dlogP_dalpha * dlogP_dgamma + d2logP_dalphadgamma)

    } else {
      effect[i] <- 1
      d_alpha[i] <- 0
      d_gamma[i] <- 0
      d2_alpha[i] <- 0
      d2_gamma[i] <- 0
      d2_alpha_gamma[i] <- 0
    }
  }

  return(list(effect = effect,
              d_alpha = d_alpha,
              d_gamma = d_gamma,
              d2_alpha = d2_alpha,
              d2_gamma = d2_gamma,
              d2_alpha_gamma = d2_alpha_gamma))
}

dast_initial_value <- function(y, D, units_m, int_mat, survey_times_data,
                               penalty,
                               mda_times, fix_alpha, power_val) {

  p <- ncol(D)
  n <- nrow(D)

  llik <- function(par) {
    beta <- par[1:p]
    if(is.null(fix_alpha)) {
      alpha <- exp(par[p+1])/(1+exp(par[p+1]))
      gamma <- exp(par[p+2])
    } else {
      alpha <- fix_alpha
      gamma <- exp(par[p+1])
    }

    fact <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                               alpha, gamma, kappa = power_val)
    eta <- as.numeric(D%*%beta)
    prob_star <- 1/(1+exp(-eta))
    prob <- fact*prob_star

    # Apply penalty for both alpha and gamma (penalty[[4]] is the gamma penalty)
    out <- -(sum(y*log(prob/(1-prob)) + units_m*log(1-prob)) -
               penalty[[1]](alpha) - penalty[[4]](gamma))
    return(out)
  }

  start <- c(runif(p+1,-1,1), mean(dist(survey_times_data)))
  est <- nlminb(start, llik)

  est$beta <- est$par[1:p]
  if(is.null(fix_alpha)) {
    est$alpha <- exp(est$par[p+1])/(1+exp(est$par[p+1]))
    est$gamma <- exp(est$par[p+2])
  } else {
    est$gamma <- exp(est$par[p+1])
  }

  return(est)
}

##' @title Fitting of decay-adjusted spatio-temporal (DAST) model
##'
##' @description
##' The function fits a decay-adjusted spatio-temporal (DAST) model using Monte Carlo maximum likelihood.
##' The DAST model allows for the incorporation of temporal decay in disease prevalence due
##' to the impact of mass drug administration (MDA). The function requires the full MDA history as detailed in the arguments below.
##'
##' Spatial and spatio-temporal dependence is specified through the \code{gp()} term in the model formula:
##' \itemize{
##'   \item \code{gp(x, y)} fits a purely spatial Gaussian process.
##'   \item \code{gp(x, y, t_gp)} fits a spatio-temporal Gaussian process, where \code{t_gp} is used as the GP temporal index.
##' }
##' In all cases, the \code{time} argument must be specified separately and provides the observation-level survey times used
##' in modelling MDA impact. These survey times may differ from the GP temporal index.
##'
##' @param formula A model formula specifying the response variable, predictors, and the GP structure through \code{gp()}.
##' @param data A \code{data.frame} or \code{sf} object containing the dataset.
##' @param den The denominator for binomial models.
##' @param time A variable in \code{data} giving the survey times of observations (required).
##' @param mda_times A vector specifying the mass drug administration (MDA) times.
##' @param int_mat Intervention matrix specifying the timing and coverage of MDA; the dimension of the matrix
##' must be \code{n * n_mda}, where \code{n} is the number of rows of \code{data} and \code{n_mda} is the length of \code{mda_times}.
##' @param penalty Optional list of length 6 specifying penalty functions for regularization.
##'   Elements \code{[[1]]}, \code{[[2]]}, \code{[[3]]} are the penalty, its first and second
##'   derivative with respect to \code{alpha}. Elements \code{[[4]]}, \code{[[5]]}, \code{[[6]]}
##'   are the penalty, its first and second derivative with respect to \code{gamma}.
##'   Any of these six elements may be set to \code{function(x) 0} to suppress penalisation
##'   for that parameter.
##' @param drop Optional value used for fixing the "drop" parameter of the MDA impact function.
##' @param power_val Value expressing the power of the MDA impact function.
##' @param crs Optional coordinate reference system (CRS) for spatial data.
##' @param convert_to_crs CRS to which spatial data should be converted.
##' @param scale_to_km Logical; whether to scale distances to kilometers (default: \code{TRUE}).
##' @param control_mcmc A list of MCMC control parameters, typically from \code{set_control_sim()}.
##' @param par0 Optional list of initial parameter values.
##' @param S_samples Number of posterior samples to retain.
##' @param return_samples Logical; whether to return posterior samples (default: \code{TRUE}).
##' @param messages Logical; whether to print messages (default: \code{TRUE}).
##' @param start_pars List of starting values for parameters.
##'
##' @return A list containing model estimates, posterior samples, and metadata, including:
##' \itemize{
##'   \item \code{y}: Response variable values.
##'   \item \code{D}: Covariate matrix.
##'   \item \code{coords}: Unique spatial coordinates.
##'   \item \code{mda_times}: MDA time points.
##'   \item \code{survey_times_data}: Survey time data from the \code{time} argument.
##'   \item \code{time}: GP temporal index if specified in \code{gp(x,y,t_gp)}.
##'   \item \code{int_mat}: Intervention matrix.
##'   \item \code{ID_coords}: Indices of spatial locations (and time if spatio-temporal GP).
##'   \item \code{re}: Random effects levels (if applicable).
##'   \item \code{ID_re}: Indices of random effects (if applicable).
##'   \item \code{power_val}: Power of the MDA impact function.
##'   \item \code{fix_tau2}: Fixed tau-squared value (if applicable).
##'   \item \code{fix_alpha}: Fixed alpha value (if applicable).
##'   \item \code{formula}: Model formula.
##'   \item \code{crs}: Coordinate reference system.
##'   \item \code{scale_to_km}: Indicator of distance scaling.
##'   \item \code{data_sf}: Processed spatial dataset.
##'   \item \code{family}: Model family (e.g., "binomial").
##'   \item \code{sst}: Logical indicator of whether a spatio-temporal GP was used.
##'   \item \code{kappa}: Smoothness parameter.
##'   \item \code{units_m}: Denominator for binomial models.
##'   \item \code{cov_offset}: Offset for covariates.
##'   \item \code{call}: Function call.
##'   \item \code{penalty}: Penalty function details (if applicable).
##'   \item \code{posterior_samples}: Posterior samples if \code{return_samples = TRUE}.
##' }
##'
##' @seealso \code{\link{set_control_sim}}, \code{\link{summary.RiskMap}}, \code{\link{to_table}}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @export
dast <- function(formula,
                 data,
                 den = NULL,
                 time,
                 mda_times, int_mat,
                 penalty = NULL,
                 drop = NULL, power_val,
                 crs = NULL, convert_to_crs = NULL,
                 scale_to_km = TRUE,
                 control_mcmc = set_control_sim(),
                 par0=NULL,
                 S_samples = NULL,
                 return_samples = TRUE,
                 messages = TRUE,
                 start_pars = list(beta = NULL,
                                   sigma2 = NULL,
                                   tau2 = NULL,
                                   phi = NULL,
                                   psi = NULL,
                                   sigma2_re = NULL,
                                   gamma = NULL,
                                   alpha = NULL)) {

  nong <- TRUE

  if(!inherits(formula,
               what = "formula", which = FALSE)) {
    stop("'formula' must be a 'formula'
         object indicating the variables of the
         model to be fitted")
  }

  inter_f <- interpret.formula(formula)

  # --- NEW: gp/time interpretation ---
  gp_terms <- inter_f$gp.spec$term
  gp_dim   <- inter_f$gp.spec$dim

  # time= (required) -> survey_times_data
  time_name <- deparse(substitute(time))
  if (time_name == "NULL") stop("You must supply a time column via the `time` argument.")
  if (!time_name %in% names(data)) stop("`time` column not found in 'data'.")
  survey_times_data <- data[[time_name]]

  # GP temporal index (optional, only if gp has 3 terms)
  sst <- FALSE
  gp_time_obs <- NULL
  if (length(gp_terms) == 1 && gp_terms[1] == "sf") {
    sst <- FALSE
  } else if (gp_dim == 3) {
    sst <- TRUE
    gp_time_obs <- data[[gp_terms[3]]]   # may differ from survey_times_data
  } else if (gp_dim == 2) {
    sst <- FALSE
  } else {
    stop("Specify gp(x, y), gp(x, y, t), or gp(sf).")
  }

  if(length(crs)>0) {
    if(!is.numeric(crs) |
       (is.numeric(crs) &
        (crs%%1!=0 | crs <0))) stop("'crs' must be a positive integer number")
  }
  if(class(data)[1]=="data.frame") {
    if(is.null(crs)) {
      warning("'crs' is set to 4326 (long/lat)")
      crs <- 4326
    }
    if(length(gp_terms)>1 && gp_terms[1] != "sf") {
      new_x <- paste(gp_terms[1],"_sf",sep="")
      new_y <- paste(gp_terms[2],"_sf",sep="")
      data[[new_x]] <-  data[[gp_terms[1]]]
      data[[new_y]] <-  data[[gp_terms[2]]]
      data <- sf::st_as_sf(data,
                           coords = c(new_x, new_y),
                           crs = crs)
    }
  }

  if(length(gp_terms) == 1 & gp_terms[1]=="sf" &
     class(data)[1]!="sf") stop("'data' must be an object of class 'sf'")

  if(class(data)[1]=="sf") {
    if(is.na(sf::st_crs(data)) & is.null(crs)) {
      stop("the CRS of the sf object passed to 'data' is missing and and is not specified through 'crs'")
    } else if(is.na(sf::st_crs(data))) {
      data <- sf::st_as_sf(data, crs = crs)
    }
  }

  # Build default no-penalty list.
  # [[1]]-[[3]]: alpha penalty and its first/second derivatives.
  # [[4]]-[[6]]: gamma penalty and its first/second derivatives.
  if(is.null(penalty)) {
    no_penalty <- TRUE
    penalty <- list(pn        = function(x) return(0),
                    pn_d1     = function(x) return(0),
                    pn_d2     = function(x) return(0),
                    pn_gamma    = function(x) return(0),
                    pn_gamma_d1 = function(x) return(0),
                    pn_gamma_d2 = function(x) return(0))
  } else {
    no_penalty <- FALSE
    if(length(penalty) != 6)
      stop("'penalty' must be a list of length 6: penalty/d1/d2 for alpha (elements 1-3) ",
           "and penalty/d1/d2 for gamma (elements 4-6).")
  }

  kappa <- inter_f$gp.spec$kappa
  if(kappa < 0) stop("kappa must be positive.")


  mf <- model.frame(inter_f$pf,data=data, na.action = na.fail)

  # Extract outcome data
  y <- as.numeric(model.response(mf))
  n <- length(y)

  # Extract covariates matrix
  D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))

  if(is.null(inter_f$offset)) {
    cov_offset <- rep(0, nrow(data))
  } else {
    cov_offset <- data[[inter_f$offset]]
  }

  # Define denominators for Binomial and Poisson distributions

  # units_m <- data$m
  if(nong) {
    do_name <- deparse(substitute(den))
    if(do_name=="NULL") {
      units_m <- rep(1, nrow(data))
      warning("'den' is assumed to be 1 for all observations \n")
    } else {
      units_m <- data[[do_name]]
    }
    if(is.integer(units_m)) units_m <- as.numeric(units_m)
    if(!is.numeric(units_m)) stop("the variable passed to `den` must be numeric")
    if(any(y > units_m)) stop("The counts identified by the outcome variable cannot be larger
                              than `den` in the case of a Binomial distribution")
    if(!inherits(control_mcmc,
                 what = "mcmc.RiskMap", which = FALSE)) {
      stop ("the argument passed to 'control_mcmc' must be an output
                                                  from the function set_control_sim; see ?set_control_sim
                                                  for more details")

    }

  }

  # survey_times_data already set from 'time' argument above

  if(length(inter_f$re.spec) > 0) {
    hr_re <- inter_f$re.spec$term
    re_names <- inter_f$re.spec$term
  } else {
    hr_re <- NULL
  }
  if(!is.null(drop)) {
    fix_alpha <- drop
  } else {
    fix_alpha <- NULL
  }

  if(!is.null(hr_re)) {
    # Define indices of random effects
    re_mf <- sf::st_drop_geometry(data[hr_re])
    re_mf_n <- re_mf

    if(any(is.na(re_mf))) stop("Missing values in the variable(s) of the random effects specified through re() ")
    names_re <- colnames(re_mf)
    n_re <- ncol(re_mf)

    ID_re <- matrix(NA, nrow = n, ncol = n_re)
    re_unique <- list()
    re_unique_f <- list()
    for(i in 1:n_re) {
      if(is.factor(re_mf[,i])) {
        re_mf_n[,i] <- as.numeric(re_mf[,i])
        re_unique[[names_re[i]]] <- 1:length(levels(re_mf[,i]))
        ID_re[, i] <- sapply(1:n,
                             function(j) which(re_mf_n[j,i]==re_unique[[names_re[i]]]))
        re_unique_f[[names_re[i]]] <-levels(re_mf[,i])
      } else if(is.numeric(re_mf[,i])) {
        re_unique[[names_re[i]]] <- unique(re_mf[,i])
        ID_re[, i] <- sapply(1:n,
                             function(j) which(re_mf_n[j,i]==re_unique[[names_re[i]]]))
        re_unique_f[[names_re[i]]] <- re_unique[[names_re[i]]]
      }
    }
    ID_re <- data.frame(ID_re)
    colnames(ID_re) <- re_names
  } else {
    n_re <- 0
    re_unique <- NULL
    ID_re <- NULL
  }


  # Extract coordinates
  if(!is.null(convert_to_crs)) {
    if(!is.numeric(convert_to_crs)) stop("'convert_to_utm' must be a numeric object")
    data <- sf::st_transform(data, crs = convert_to_crs)
    crs <- convert_to_crs
  }

  if(messages) message("The CRS used is ", as.list(sf::st_crs(data))$input, "\n")

  coords_o <- sf::st_coordinates(data)
  if(sst) {
    coords_time <- unique(cbind(coords_o, gp_time_obs))
    coords <- coords_time[,1:2, drop = FALSE]
    time_gp <- coords_time[,3]
  } else {
    coords <- unique(coords_o)
    time_gp <- NULL
  }

  m <- nrow(coords_o)
  if(sst) {
    ID_coords <- sapply(1:m, function(i)
      which(coords_o[i,1]==coords[,1] &
              coords_o[i,2]==coords[,2] &
              gp_time_obs[i]==time_gp))
  } else {
    ID_coords <- sapply(1:m, function(i)
      which(coords_o[i,1]==coords[,1] &
              coords_o[i,2]==coords[,2]))
  }
  s_unique <- unique(ID_coords)

  fix_tau2 <- inter_f$gp.spec$nugget


  if(scale_to_km) {
    coords_o <- coords_o/1000
    coords <- coords/1000
    if(messages) message("Distances between locations are computed in kilometers ")
  } else {
    if(messages) message("Distances between locations are computed in meters ")
  }

  if(is.null(start_pars) & !is.null(par0)) {
    start_pars <- par0
    if(length(par0$beta)!=ncol(D)) stop("the values passed to `beta` in par0 do not match the
                                        variables specified in the formula")
  }

  if(is.null(start_pars$beta)) {
    aux_data <- data.frame(y=y, units_m = units_m, D[,-1])
    if(length(cov_offset)==1) cov_offset_aux <- rep(cov_offset, n)
    glm_fitted <- glm(cbind(y, units_m - y) ~ ., offset = cov_offset,
                      data = aux_data, family = binomial)
    start_pars$beta <- stats::coef(glm_fitted)
  }

  if(is.null(start_pars$gamma) | (is.null(drop) & is.null(start_pars$alpha)) ) {
    dast_i <- dast_initial_value(y, D, units_m, int_mat = int_mat, survey_times_data,
                                 fix_alpha = fix_alpha,penalty=penalty,
                                 mda_times, power_val = power_val)
    start_pars$beta <- dast_i$beta
    if(is.null(drop)) {
      start_pars$alpha <- dast_i$alpha
    }
    start_pars$gamma <- dast_i$gamma
  } else {
    if(length(start_pars$beta)!=ncol(D)) stop("number of starting values provided
                                              for 'beta' do not match the number of
                                              covariates specified in the model,
                                              including the intercept")
  }

  if(is.null(start_pars$sigma2)) {
    start_pars$sigma2 <- 1
  } else {
    if(start_pars$sigma2<0) stop("the starting value for sigma2 must be positive")
  }

  if(is.null(start_pars$phi)) {
    start_pars$phi <- stats::quantile(stats::dist(coords),0.1)
  } else {
    if(start_pars$phi<0) stop("the starting value for phi must be positive")
  }

  if(is.null(fix_tau2)) {
    if(is.null(start_pars$tau2)) {
      start_pars$tau2 <- 1
    } else {
      if(start_pars$tau2<0) stop("the starting value for tau2 must be positive")
    }
  }
  if(sst) {
    if(is.null(start_pars$psi)) {
      # base distances from GP time index (not survey time)
      psi_base <- as.numeric(unique(gp_time_obs))
      if (length(psi_base) > 1) {
        start_pars$psi <- stats::quantile(stats::dist(psi_base), 0.1)
      } else {
        start_pars$psi <- 1
      }
    } else {
      if(start_pars$psi<0) stop("the starting value for psi must be positive")
    }
  }

  if(n_re > 0) {
    if(is.null(start_pars$sigma2_re)) {
      start_pars$sigma2_re <- rep(1,n_re)
    } else {
      if(length(start_pars$sigma2_re)!=n_re) stop("starting values for 'sigma2_re' do not
                                       match the number of specified unstructured
                                       random effects")
      if(any(start_pars$sigma2_re<0)) stop("all the starting values for sigma2_re must be positive")
    }
  }



  if(is.null(par0)) {
    par0 <- start_pars
  }
  res <- dast_fit(y = y, D = D, coords = coords, time = time_gp,
                  units_m = units_m,
                  mda_times = mda_times, survey_times_data = survey_times_data,
                  sst = sst,
                  int_mat = int_mat,
                  kappa = inter_f$gp.spec$kappa,
                  penalty = penalty,
                  ID_coords = ID_coords, ID_re = ID_re, s_unique = s_unique, re_unique = re_unique,
                  fix_tau2 = fix_tau2, fix_alpha = fix_alpha,
                  return_samples = return_samples,
                  par0 = par0, cov_offset = cov_offset,
                  power_val = power_val,
                  start_beta = start_pars$beta,
                  start_alpha = start_pars$alpha,
                  start_gamma = start_pars$gamma,
                  start_cov_pars = c(start_pars$sigma2,
                                     start_pars$phi,
                                     start_pars$tau2,
                                     start_pars$sigma2_re),
                  start_psi = start_pars$psi,
                  control_mcmc = control_mcmc,
                  messages = messages)



  res$y <- y
  res$D <- D
  res$coords <- coords
  res$sst <- sst
  if(sst) res$time <- time_gp
  res$mda_times <- mda_times
  res$survey_times_data <- survey_times_data
  res$int_mat <- int_mat
  res$ID_coords <- ID_coords
  if(n_re>0) {
    res$re <- re_unique_f
    res$ID_re <- as.data.frame(ID_re)
    colnames(res$ID_re) <- names_re
  }
  res$power_val <- power_val
  res$fix_tau2 <- fix_tau2
  res$fix_alpha <- fix_alpha
  res$formula <- formula
  if(!is.null(convert_to_crs)) {
    crs <- convert_to_crs
  } else {
    crs <- sf::st_crs(data)$input
  }
  if(no_penalty) {
    res$penalty <- NULL
  } else {
    res$penalty <- penalty
  }
  res$crs <- crs
  res$scale_to_km <- scale_to_km
  res$data_sf <- data
  res$family <- "binomial"
  res$kappa <- kappa
  if(nong) res$units_m <- units_m
  res$cov_offset <- cov_offset
  res$call <- match.call()


  invlink <- NULL
  if (is.null(invlink)) {
    inv <- function(eta) stats::plogis(eta)
    d1  <- function(eta) { p <- inv(eta); p * (1 - p) }
    d2  <- function(eta) { p <- inv(eta); d <- p * (1 - 2 * p) }
    invlink <- list(inv = inv, d1 = d1, d2 = d2, name = "canonical")
  }
  res$linkf <- invlink

  return(res)
}


##' Simulation from Decay-adjusted Spatio-temporal (DAST) models
##'
##' Simulates data from a fitted DAST model object (output from \code{dast})
##' or from user-specified DAST parameters.
##'
##' @param n_sim Number of simulations.
##' @param model_fit Optional fitted DAST model object of class \code{RiskMap}.
##' If supplied, it overrides model specification arguments.
##' @param formula Model formula including a \code{gp()} term.
##' @param data Data frame or \code{sf} object used for simulation.
##' @param den Binomial denominator variable. If missing, it is assumed to be 1.
##' @param time Survey-time information. For simulations from scratch this can be a
##' column in \code{data} (unquoted name or character string) or a numeric vector
##' of length \code{nrow(data)}.
##' @param mda_times Vector of MDA times.
##' @param int_mat Intervention matrix (n x length(mda_times)) with coverage values.
##' @param power_val Power value for the MDA impact function.
##' @param cov_offset Optional offset for the linear predictor.
##' @param crs Coordinate reference system (CRS) code.
##' @param convert_to_crs Optional CRS to transform coordinates to before simulation.
##' @param scale_to_km Logical; if TRUE distances are computed in kilometers.
##' @param sim_pars List of simulation parameters. Used only when \code{model_fit}
##' is \code{NULL}. Includes \code{beta}, \code{sigma2}, \code{tau2}, \code{phi},
##' \code{psi}, \code{sigma2_re}, \code{alpha}, and \code{gamma}.
##' @param messages Logical; if TRUE print progress messages.
##'
##' @return A list with simulated outcomes in \code{data_sim} as an
##' \code{n x n_sim} matrix (rows are observations, columns are simulations),
##' simulated latent components and parameter values used for simulation.
##' @export
dast_sim <- function(n_sim,
                     model_fit = NULL,
                     formula = NULL,
                     data = NULL,
                     den = NULL,
                     time = NULL,
                     mda_times = NULL,
                     int_mat = NULL,
                     power_val = NULL,
                     cov_offset = NULL,
                     crs = NULL,
                     convert_to_crs = NULL,
                     scale_to_km = TRUE,
                     sim_pars = list(beta = NULL,
                                     sigma2 = NULL,
                                     tau2 = NULL,
                                     phi = NULL,
                                     psi = NULL,
                                     sigma2_re = NULL,
                                     alpha = NULL,
                                     gamma = NULL),
                     messages = TRUE) {

  if (!is.null(model_fit)) {
    if (!inherits(model_fit, what = "RiskMap", which = FALSE)) {
      stop("'model_fit' must be of class 'RiskMap'")
    }
    if (is.null(model_fit$power_val)) {
      stop("'model_fit' does not look like an output from 'dast'")
    }
    formula <- as.formula(model_fit$formula)
    data <- model_fit$data_sf
    crs <- model_fit$crs
    convert_to_crs <- model_fit$convert_to_crs
    scale_to_km <- model_fit$scale_to_km
    power_val <- model_fit$power_val
    mda_times <- model_fit$mda_times
    int_mat <- model_fit$int_mat
    time <- model_fit$survey_times_data
  }

  if (!inherits(formula, what = "formula", which = FALSE)) {
    stop("'formula' must be a 'formula'\n         object indicating the variables of the\n         model to be fitted")
  }

  inter_f <- interpret.formula(formula)
  gp_terms <- inter_f$gp.spec$term
  gp_dim <- inter_f$gp.spec$dim

  # Survey time and optional GP time index
  sst <- FALSE
  gp_time_obs <- NULL

  if (is.null(model_fit)) {
    time_name <- deparse(substitute(time))
    if (time_name != "NULL" && time_name %in% names(data)) {
      survey_times_data <- data[[time_name]]
    } else if (is.character(time) && length(time) == 1 && time %in% names(data)) {
      survey_times_data <- data[[time]]
    } else if (length(time) == nrow(data)) {
      survey_times_data <- as.numeric(time)
    } else {
      stop("You must supply survey times through `time` as a column in 'data' or as a numeric vector of length nrow(data).")
    }
  } else {
    survey_times_data <- as.numeric(time)
  }

  if (length(gp_terms) == 1 && gp_terms[1] == "sf") {
    sst <- FALSE
  } else if (gp_dim == 3) {
    sst <- TRUE
    gp_time_obs <- data[[gp_terms[3]]]
  } else if (gp_dim == 2) {
    sst <- FALSE
  } else {
    stop("Specify gp(x, y), gp(x, y, t), or gp(sf).")
  }

  if(length(crs)>0) {
    if (is.character(crs) && length(crs) == 1) {
      crs_chr <- trimws(crs)
      if (grepl("^EPSG:[0-9]+$", crs_chr, ignore.case = TRUE)) {
        crs_chr <- sub("^EPSG:", "", crs_chr, ignore.case = TRUE)
      }
      if (grepl("^[0-9]+$", crs_chr)) {
        crs <- as.numeric(crs_chr)
      } else {
        stop("'crs' must be a positive integer number or a string like 'EPSG:<code>'")
      }
    }
    if(!is.numeric(crs) |
       (is.numeric(crs) &
          (crs%%1!=0 | crs <0))) stop("'crs' must be a positive integer number")
  }
  if(class(data)[1]=="data.frame") {
    if(is.null(crs)) {
      warning("'crs' is set to 4326 (long/lat)")
      crs <- 4326
    }
    if(length(gp_terms)>1 && gp_terms[1] != "sf") {
      new_x <- paste(gp_terms[1],"_sf",sep="")
      new_y <- paste(gp_terms[2],"_sf",sep="")
      data[[new_x]] <-  data[[gp_terms[1]]]
      data[[new_y]] <-  data[[gp_terms[2]]]
      data <- sf::st_as_sf(data,
                           coords = c(new_x, new_y),
                           crs = crs)
    }
  }

  if(length(gp_terms) == 1 & gp_terms[1]=="sf" &
     class(data)[1]!="sf") stop("'data' must be an object of class 'sf'")

  if(class(data)[1]=="sf") {
    if(is.na(sf::st_crs(data)) & is.null(crs)) {
      stop("the CRS of the sf object passed to 'data' is missing and and is not specified through 'crs'")
    } else if(is.na(sf::st_crs(data))) {
      data <- sf::st_as_sf(data, crs = crs)
    }
  }

  kappa <- inter_f$gp.spec$kappa
  if(kappa < 0) stop("kappa must be positive.")

  mf <- model.frame(inter_f$pf,data=data, na.action = na.fail)
  D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))
  n <- nrow(D)

  if (is.null(cov_offset)) {
    if (is.null(inter_f$offset)) {
      cov_offset <- rep(0, n)
    } else {
      cov_offset <- data[[inter_f$offset]]
    }
  }

  # denominator
  if(!is.null(den)) {
    do_name <- deparse(substitute(den))
    if (do_name != "NULL") {
      units_m <- data[[do_name]]
    } else {
      units_m <- den
    }
  } else if (!is.null(model_fit)) {
    units_m <- model_fit$units_m
  } else {
    units_m <- rep(1, n)
    warning("'den' is assumed to be 1 for all observations")
  }
  if(is.integer(units_m)) units_m <- as.numeric(units_m)
  if(!is.numeric(units_m)) stop("the variable passed to `den` must be numeric")

  if (length(survey_times_data) != n) {
    stop("'time' must have length equal to the number of observations")
  }
  if (is.null(mda_times) || length(mda_times) < 1) {
    stop("'mda_times' must be provided")
  }
  if (is.null(int_mat)) {
    stop("'int_mat' must be provided")
  }
  int_mat <- as.matrix(int_mat)
  if (nrow(int_mat) != n || ncol(int_mat) != length(mda_times)) {
    stop("'int_mat' must have dimensions nrow(data) x length(mda_times)")
  }
  if (is.null(power_val)) stop("'power_val' must be provided")

  if(length(inter_f$re.spec) > 0) {
    hr_re <- inter_f$re.spec$term
  } else {
    hr_re <- NULL
  }

  if(!is.null(hr_re)) {
    re_mf <- sf::st_drop_geometry(data[hr_re])
    re_mf_n <- re_mf

    if(any(is.na(re_mf))) stop("Missing values in the variable(s) of the random effects specified through re()")
    names_re <- colnames(re_mf)
    n_re <- ncol(re_mf)

    ID_re <- matrix(NA, nrow = n, ncol = n_re)
    re_unique <- list()
    re_unique_f <- list()
    for(i in 1:n_re) {
      if(is.factor(re_mf[,i])) {
        re_mf_n[,i] <- as.numeric(re_mf[,i])
        re_unique[[names_re[i]]] <- 1:length(levels(re_mf[,i]))
        ID_re[, i] <- sapply(1:n,
                             function(j) which(re_mf_n[j,i]==re_unique[[names_re[i]]]))
        re_unique_f[[names_re[i]]] <-levels(re_mf[,i])
      } else if(is.numeric(re_mf[,i])) {
        re_unique[[names_re[i]]] <- unique(re_mf[,i])
        ID_re[, i] <- sapply(1:n,
                             function(j) which(re_mf_n[j,i]==re_unique[[names_re[i]]]))
        re_unique_f[[names_re[i]]] <- re_unique[[names_re[i]]]
      }
    }
    ID_re <- data.frame(ID_re)
    colnames(ID_re) <- names_re
  } else {
    n_re <- 0
    re_unique <- NULL
    ID_re <- NULL
  }

  if(!is.null(convert_to_crs)) {
    if(!is.numeric(convert_to_crs)) stop("'convert_to_crs' must be a numeric object")
    data <- sf::st_transform(data, crs = convert_to_crs)
    crs <- convert_to_crs
  }

  if(messages) message("The CRS used is ", as.list(sf::st_crs(data))$input, "\n")

  coords_o <- sf::st_coordinates(data)
  if(sst) {
    coords_time <- unique(cbind(coords_o, gp_time_obs))
    coords <- coords_time[,1:2, drop = FALSE]
    time_gp <- coords_time[,3]
  } else {
    coords <- unique(coords_o)
    time_gp <- NULL
  }

  m <- nrow(coords_o)
  if(sst) {
    ID_coords <- sapply(1:m, function(i)
      which(coords_o[i,1]==coords[,1] &
              coords_o[i,2]==coords[,2] &
              gp_time_obs[i]==time_gp))
  } else {
    ID_coords <- sapply(1:m, function(i)
      which(coords_o[i,1]==coords[,1] &
              coords_o[i,2]==coords[,2]))
  }

  if(scale_to_km) {
    coords <- coords/1000
    if(messages) message("Distances between locations are computed in kilometers \n")
  } else {
    if(messages) message("Distances between locations are computed in meters \n")
  }

  if (!is.null(model_fit)) {
    par_hat <- coef(model_fit)
    beta <- par_hat$beta
    sigma2 <- par_hat$sigma2
    phi <- par_hat$phi
    if (is.null(model_fit$fix_tau2)) {
      tau2 <- par_hat$tau2
    } else {
      tau2 <- model_fit$fix_tau2
    }
    if (is.null(model_fit$fix_alpha)) {
      alpha <- par_hat$alpha
    } else {
      alpha <- model_fit$fix_alpha
    }
    gamma <- par_hat$gamma
    if (sst) psi <- par_hat$psi
    if(n_re > 0) sigma2_re <- par_hat$sigma2_re
  } else {
    beta <- sim_pars$beta
    sigma2 <- sim_pars$sigma2
    tau2 <- sim_pars$tau2
    phi <- sim_pars$phi
    alpha <- sim_pars$alpha
    gamma <- sim_pars$gamma
    psi <- sim_pars$psi
    sigma2_re <- sim_pars$sigma2_re
  }

  if(is.null(beta)) stop("'beta' is missing")
  if(length(beta)!=ncol(D)) stop("the number of values provided for 'beta' does not match the number of covariates")
  if(is.null(sigma2)) stop("'sigma2' is missing")
  if(is.null(phi)) stop("'phi' is missing")
  if(is.null(tau2)) stop("'tau2' is missing")
  if(is.null(alpha)) stop("'alpha' is missing")
  if(is.null(gamma)) stop("'gamma' is missing")
  if(sst && is.null(psi)) stop("'psi' is missing for spatio-temporal DAST simulations")
  if(n_re>0) {
    if(is.null(sigma2_re)) stop("'sigma2_re' is missing")
    if(length(sigma2_re)!=n_re) stop("the values passed to 'sigma2_re' do not match the number of random effects")
  }

  if(sst) {
    Sigma <- sigma2*matern_cor(dist(coords), phi = phi, kappa = kappa,
                                                return_sym_matrix = TRUE) *
      matern_cor(dist(time_gp), phi = psi, kappa = 0.5,
                 return_sym_matrix = TRUE)
  } else {
    Sigma <- sigma2*matern_cor(dist(coords), phi = phi, kappa = kappa,
                               return_sym_matrix = TRUE)
  }
  diag(Sigma) <- diag(Sigma) + tau2

  Sigma_sroot <- t(chol(Sigma))
  S_sim <- t(sapply(1:n_sim, function(i) Sigma_sroot %*% rnorm(nrow(coords))))

  if(n_re > 0) {
    re_names <- if(!is.null(model_fit)) names(model_fit$re) else inter_f$re.spec$term
    dim_re <- sapply(1:n_re, function(j) length(re_unique[[j]]))
    re_sim <- vector("list", n_sim)
    for(i in 1:n_sim) {
      re_sim[[i]] <- list()
      for(j in 1:n_re) {
        re_sim[[i]][[paste(re_names[j])]] <- rnorm(dim_re[j]) * sqrt(sigma2_re[j])
      }
    }
  }

  eta_sim <- t(sapply(1:n_sim, function(i) as.numeric(D %*% beta + cov_offset + S_sim[i,][ID_coords])))

  if(n_re > 0) {
    for(i in 1:n_sim) {
      for(j in 1:n_re) {
        eta_sim[i,] <- eta_sim[i,] + re_sim[[i]][[paste(re_names[j])]][ID_re[,j]]
      }
    }
  }

  mda_effect <- compute_mda_effect(survey_times_data = survey_times_data,
                                   mda_times = mda_times,
                                   intervention = int_mat,
                                   alpha = alpha,
                                   gamma = gamma,
                                   kappa = power_val)

  y_sim <- matrix(NA, nrow = n_sim, ncol = n)
  for(i in 1:n_sim) {
    prob_i <- stats::plogis(eta_sim[i,]) * mda_effect
    prob_i <- pmax(0, pmin(prob_i, 1))
    y_sim[i,] <- stats::rbinom(n = n, size = units_m, prob = prob_i)
  }

  data_sim <- t(y_sim)
  response_name <- if (!is.null(inter_f$response) && nzchar(inter_f$response)) {
    inter_f$response
  } else {
    "y"
  }
  colnames(data_sim) <- paste(response_name, "_sim", 1:n_sim, sep = "")

  out <- list(data_sim = data_sim,
              S_sim = S_sim,
              lin_pred_sim = eta_sim,
              mda_effect = mda_effect,
              beta = beta,
              sigma2 = sigma2,
              tau2 = tau2,
              phi = phi,
              alpha = alpha,
              gamma = gamma,
              power_val = power_val,
              units_m = units_m)
  if (sst) out$psi <- psi
  if (n_re > 0) {
    out$sigma2_re <- sigma2_re
    out$re_sim <- re_sim
  }
  return(out)
}

##' Internal Posterior Predictive Validation for DAST Fits
##'
##' Internal utility to validate DAST model fit quality by comparing observed
##' outcomes to simulations generated from the fitted DAST data-generating process.
##' The function computes discrepancy summaries at global, temporal,
##' and IU aggregation levels, and optionally returns diagnostic plots.
##'
##' Interpretation guide:
##' \itemize{
##'   \item \code{summaries$global$pvals}: posterior predictive p-values (PPP) for global
##'   discrepancy statistics. Values close to 0 or 1 indicate mismatch between
##'   data and model-implied simulations for the corresponding summary.
##'   \item \code{summaries$global$pit_summary}: PIT mean/variance; departures from roughly
##'   Uniform(0,1) moments suggest calibration issues.
##'   \item \code{summaries$global$z_summary}: quantiles of standardized residuals
##'   \eqn{(y_i-\mu_i)/\sqrt{v_i}} based on simulation moments. Heavy tails or
##'   skewness can indicate poor dispersion/shape fit.
##'   \item \code{summaries$temporal$pvals}: year-level trajectory mismatch summaries.
##'   \item \code{summaries$iu$pvals}: IU and IU-year support-aware mismatch summaries.
##'   \item \code{summaries$iu$iu_support}: IU-level support table (number of surveys,
##'   total examined, observed pooled/unweighted prevalence and inclusion flag).
##'   \item \code{summaries$iu$diagnostics}: counts summarising IU inclusion and zero vs
##'   non-zero observed records, including out-of-interval IU counts.
##' }
##'
##' @param y_obs Numeric vector of observed positive counts.
##' @param m Numeric vector of denominators (examined counts).
##' @param y_sim Numeric matrix of simulated counts with dimension
##' \code{length(y_obs) x nsim} (rows are observations, columns are simulation draws).
##' @param year Observation year (numeric/integer/factor/character coercible to numeric).
##' @param iu IU identifier vector of same length as \code{y_obs}.
##' @param m_bins Breaks used to create denominator groups for grouped PPC plots.
##' @param iu_min_total_m Minimum total examined count required to include an IU
##' in IU-level summaries.
##' @param iu_year_min_total_m Minimum total examined count required to include an
##' IU-year cell in IU-year summaries.
##' @param prob Central predictive interval probability used in coverage summaries
##' and envelopes.
##' @param n_overlay Number of simulation draws to overlay in PPC plots.
##' @param seed Integer random seed used for randomized PIT and draw subsampling.
##' @param make_plots Logical; if \code{TRUE}, return PPC plots when plotting
##' dependencies are available.
##' @param prevalence_estimand Which prevalence estimand to report and plot:
##' \code{"pooled"}, \code{"unweighted"} or \code{"both"}.
##'
##' @return An object of class \code{"ppc_dast_dgp"} with components:
##' \describe{
##'   \item{inputs}{Input metadata and key settings.}
##'   \item{summaries}{A nested list with \code{global}, \code{temporal}, and \code{iu}
##'   summary outputs.}
##'   \item{plots}{A nested list with \code{global}, \code{temporal}, and \code{iu}
##'   diagnostic plots (possibly empty).}
##' }
##'
##' Year-level and IU-level prevalence trajectories use the examined-weighted
##' prevalence (pooled prevalence), computed as
##' \eqn{\sum_i y_i / \sum_i m_i} within each aggregation stratum.
##' This is different from the unweighted mean of location-level prevalences.
##'
##' For IU-level summaries and plots, an IU is included when it has at least one
##' observed record and total examined count \eqn{\sum_i m_i \ge iu\_min\_total\_m}.
##'
##' @details The two-sided PPP is computed as a median-centered tail area:
##' \deqn{PPP = \frac{1}{S}\sum_{s=1}^{S}
##' \mathbf{1}\left(\left|T_s^{rep}-\mathrm{median}(T^{rep})\right|
##' \ge \left|T^{obs}-\mathrm{median}(T^{rep})\right|\right)}
##' This is robust for asymmetric predictive distributions and does not assume
##' normality of discrepancy statistics.
##'
##' When \code{make_plots = TRUE}, \pkg{bayesplot} is required for PPC plot
##' primitives. If \pkg{bayesplot} is unavailable, the function returns numeric
##' summaries and sets \code{out$plot_note} with a warning message.
##' Prevalence values are computed internally on the 0-1 scale but shown in
##' plots as percentages for readability.
##'
##' The IU mean prevalence PPC plot is a forest-style diagnostic:
##' observed IU prevalence is shown as points, simulated central intervals are
##' shown as horizontal ranges, and observed point fill color encodes the number
##' of survey records contributing to each IU mean. A companion plot
##' \code{plots$iu$iu_mean_prev_misfit} shows only IUs where observed prevalence
##' falls outside the simulated central interval.
##'
##' @examples
##' \dontrun{
##' out <- dast_sim(n_sim = 250, model_fit = dast_with_pen, messages = FALSE)
##' ppc_all <- ppc_dast_dgp(
##'   y_obs = dast_with_pen$data_sf$positive,
##'   m = dast_with_pen$data_sf$examined,
##'   y_sim = out$data_sim,
##'   year = dast_with_pen$data_sf$year,
##'   iu = dast_with_pen$data_sf$IUs_NAME,
##'   prob = 0.9
##' )
##' print(ppc_all)
##' ppc_all$plots$temporal$year_traj
##' ppc_all$plots$temporal$year_diffs
##' ppc_all$plots$iu$iu_mean_prev
##' ppc_all$plots$iu$iu_mean_prev_misfit
##' ppc_all$plots$iu$iu_year_cell_ecdf
##' ppc_all$plots$global$pit_ecdf
##' ppc_all$plots$global$stat_hist$binom_deviance_plugin_mean
##' }
##'
##' @keywords internal
ppc_dast_dgp <- function(
    y_obs, m, y_sim, year, iu,
    m_bins = c(0, 23.5, 30.5, 49.5, Inf),
    iu_min_total_m = 50,
    iu_year_min_total_m = 25,
    prob = 0.9,
    n_overlay = 50,
    seed = 1L,
    make_plots = TRUE,
    prevalence_estimand = c("both", "pooled", "unweighted")
) {
  set.seed(seed)
  prevalence_estimand <- match.arg(prevalence_estimand)

  if (!is.numeric(y_obs) || !is.numeric(m) || !is.matrix(y_sim)) {
    stop("'y_obs' and 'm' must be numeric vectors and 'y_sim' must be a matrix.")
  }
  if (length(iu) != length(y_obs)) {
    stop("'iu' must have the same length as 'y_obs'.")
  }
  n <- length(y_obs)
  if (length(m) != n) stop("'m' must have the same length as 'y_obs'.")
  if (nrow(y_sim) != n) stop("'y_sim' must be n x nsim (rows are observations).")
  if (ncol(y_sim) < 2) stop("'y_sim' must contain at least two simulations (nsim >= 2).")
  nsim <- ncol(y_sim)

  if (anyNA(y_obs) || anyNA(m) || anyNA(y_sim) || anyNA(year) || anyNA(iu)) {
    stop("Missing values are not allowed; filter data before calling 'ppc_dast_dgp'.")
  }
  if (any(m <= 0)) stop("'m' must be strictly positive.")
  if (any(y_obs < 0 | y_obs > m)) stop("'y_obs' must satisfy 0 <= y_obs <= m.")
  if (any(y_sim < 0) || any(y_sim > matrix(m, n, nsim))) {
    stop("'y_sim' must satisfy 0 <= y_sim <= m for each row.")
  }
  if (!is.numeric(prob) || length(prob) != 1 || prob <= 0 || prob >= 1) {
    stop("'prob' must be a scalar in (0, 1).")
  }

  year_num <- suppressWarnings(as.numeric(as.character(year)))
  if (anyNA(year_num)) {
    stop("'year' must be coercible to numeric without introducing missing values.")
  }
  year <- year_num
  iu <- as.factor(iu)

  yrep <- t(y_sim)
  p_obs <- y_obs / m
  p_rep <- sweep(yrep, 2, m, "/")

  m_group <- cut(m, breaks = m_bins, include.lowest = TRUE, right = FALSE)

  ord <- order(year)
  year_ord <- year[ord]
  y_obs_ord <- y_obs[ord]
  m_ord <- m[ord]
  yrep_ord <- yrep[, ord, drop = FALSE]
  iu_ord <- iu[ord]

  safe_var <- function(x) {
    if (length(x) <= 1) return(0)
    stats::var(x)
  }
  assert_all_finite <- function(x, label) {
    if (any(!is.finite(x))) {
      stop(sprintf("Non-finite values found in %s. Check inputs and simulated values.", label))
    }
    invisible(TRUE)
  }

  ppp_two_sided <- function(T_rep, T_obs) {
    T_rep <- T_rep[is.finite(T_rep)]
    if (length(T_rep) == 0 || !is.finite(T_obs)) return(NA_real_)
    med <- stats::median(T_rep)
    mean(abs(T_rep - med) >= abs(T_obs - med))
  }

  estimand_levels <- switch(prevalence_estimand,
                            both = c("pooled", "unweighted"),
                            pooled = "pooled",
                            unweighted = "unweighted")

  prevalence_by_stratum <- function(y_vec, m_vec, group, levels_group, type = c("pooled", "unweighted")) {
    type <- match.arg(type)
    sapply(levels_group, function(g) {
      idx <- which(group == g)
      if (length(idx) == 0) return(NA_real_)
      if (type == "pooled") {
        den <- sum(m_vec[idx])
        if (den <= 0) return(NA_real_)
        sum(y_vec[idx]) / den
      } else {
        mean(y_vec[idx] / m_vec[idx])
      }
    })
  }

  T_wmean_obs <- sum(y_obs) / sum(m)
  T_wmean_rep <- rowSums(yrep) / sum(m)
  T_mean_obs <- mean(p_obs)
  T_mean_rep <- rowMeans(p_rep)
  T_var_obs <- safe_var(p_obs)
  T_var_rep <- apply(p_rep, 1, safe_var)
  T_zeros_obs <- mean(y_obs == 0)
  T_zeros_rep <- rowMeans(yrep == 0)
  eps0 <- 0.01
  T_zeroish_obs <- mean(p_obs <= eps0)
  T_zeroish_rep <- rowMeans(p_rep <= eps0)
  T_maxp_obs <- max(p_obs)
  T_maxp_rep <- apply(p_rep, 1, max)

  mu_i <- colMeans(yrep)
  v_i <- apply(yrep, 2, safe_var)
  v_i[v_i < 1e-12] <- 1e-12
  z_obs <- (y_obs - mu_i) / sqrt(v_i)
  T_maxabsz_obs <- max(abs(z_obs))
  T_maxabsz_rep <- apply(yrep, 1, function(yr) max(abs((yr - mu_i) / sqrt(v_i))))

  alpha <- 1 - prob
  q_lo <- apply(yrep, 2, stats::quantile, probs = alpha / 2, type = 8)
  q_hi <- apply(yrep, 2, stats::quantile, probs = 1 - alpha / 2, type = 8)
  T_cov_obs <- mean(y_obs >= q_lo & y_obs <= q_hi)
  T_cov_rep <- apply(yrep, 1, function(yr) mean(yr >= q_lo & yr <= q_hi))

  mu_p <- pmin(pmax(mu_i / m, 1e-9), 1 - 1e-9)
  binom_dev <- function(y, m, p) {
    y0 <- y
    y1 <- m - y
    a <- ifelse(y0 == 0, 0, y0 * log(y0 / (m * p)))
    b <- ifelse(y1 == 0, 0, y1 * log(y1 / (m * (1 - p))))
    2 * (a + b)
  }
  T_dev_obs <- sum(binom_dev(y_obs, m, mu_p))
  T_dev_rep <- apply(yrep, 1, function(yr) sum(binom_dev(yr, m, mu_p)))

  rand_pit_i <- function(y_i, sims_i) {
    Fy <- mean(sims_i <= y_i)
    Fym <- if (y_i <= 0) 0 else mean(sims_i <= (y_i - 1))
    Fym + stats::runif(1) * (Fy - Fym)
  }
  pit <- vapply(seq_len(n), function(i) rand_pit_i(y_obs[i], y_sim[i, ]), numeric(1))

  global_pvals <- data.frame(
    stat = c("pooled_mean_prev", "unweighted_mean_prev", "var_prev",
             "zero_fraction_counts", "zeroish_prev<=0.01", "max_prev",
             "binom_deviance(plugin_mean)", "max_abs_standardised_resid",
             sprintf("PI_coverage_%d%%", round(100 * prob))),
    obs = c(T_wmean_obs, T_mean_obs, T_var_obs,
            T_zeros_obs, T_zeroish_obs, T_maxp_obs,
            T_dev_obs, T_maxabsz_obs, T_cov_obs),
    ppp_two_sided = c(
      ppp_two_sided(T_wmean_rep, T_wmean_obs),
      ppp_two_sided(T_mean_rep, T_mean_obs),
      ppp_two_sided(T_var_rep, T_var_obs),
      ppp_two_sided(T_zeros_rep, T_zeros_obs),
      ppp_two_sided(T_zeroish_rep, T_zeroish_obs),
      ppp_two_sided(T_maxp_rep, T_maxp_obs),
      ppp_two_sided(T_dev_rep, T_dev_obs),
      ppp_two_sided(T_maxabsz_rep, T_maxabsz_obs),
      ppp_two_sided(T_cov_rep, T_cov_obs)
    ),
    stringsAsFactors = FALSE
  )

  yrs <- sort(unique(year_ord))
  temporal_summaries <- list()
  temporal_pvals <- data.frame(stat = character(0), ppp_two_sided = numeric(0), stringsAsFactors = FALSE)

  for (est in estimand_levels) {
    obs_traj_est <- prevalence_by_stratum(y_obs_ord, m_ord, year_ord, yrs, est)
    rep_traj_est <- t(apply(yrep_ord, 1, function(yr) prevalence_by_stratum(yr, m_ord, year_ord, yrs, est)))
    if (is.null(dim(rep_traj_est))) rep_traj_est <- matrix(rep_traj_est, nrow = nsim)

    rep_lo_est <- apply(rep_traj_est, 2, stats::quantile, probs = alpha / 2, type = 8, na.rm = TRUE)
    rep_md_est <- apply(rep_traj_est, 2, stats::quantile, probs = 0.5, type = 8, na.rm = TRUE)
    rep_hi_est <- apply(rep_traj_est, 2, stats::quantile, probs = 1 - alpha / 2, type = 8, na.rm = TRUE)

    T_year_sse_obs_est <- sum((obs_traj_est - rep_md_est)^2, na.rm = TRUE)
    T_year_sse_rep_est <- apply(rep_traj_est, 1, function(tr) sum((tr - rep_md_est)^2, na.rm = TRUE))

    d_obs_est <- diff(obs_traj_est)
    d_rep_est <- t(apply(rep_traj_est, 1, diff))
    if (is.null(dim(d_rep_est))) d_rep_est <- matrix(d_rep_est, nrow = nsim)
    d_md_est <- apply(d_rep_est, 2, stats::quantile, probs = 0.5, type = 8, na.rm = TRUE)

    T_year_diff_sse_obs_est <- sum((d_obs_est - d_md_est)^2, na.rm = TRUE)
    T_year_diff_sse_rep_est <- apply(d_rep_est, 1, function(dr) sum((dr - d_md_est)^2, na.rm = TRUE))

    temporal_pvals <- rbind(
      temporal_pvals,
      data.frame(
        stat = c(paste0("year_traj_SSE_", est), paste0("year_diff_SSE_", est)),
        ppp_two_sided = c(
          ppp_two_sided(T_year_sse_rep_est, T_year_sse_obs_est),
          ppp_two_sided(T_year_diff_sse_rep_est, T_year_diff_sse_obs_est)
        ),
        stringsAsFactors = FALSE
      )
    )

    temporal_summaries[[est]] <- list(
      years = yrs,
      obs_traj = obs_traj_est,
      rep_traj = rep_traj_est,
      rep_lo = rep_lo_est,
      rep_md = rep_md_est,
      rep_hi = rep_hi_est,
      d_obs = d_obs_est,
      d_rep = d_rep_est,
      d_md = d_md_est
    )
  }

  iu_n_obs <- as.integer(table(iu_ord))
  names(iu_n_obs) <- names(table(iu_ord))
  iu_tot_m <- tapply(m_ord, iu_ord, sum)
  iu_tot_y <- tapply(y_obs_ord, iu_ord, sum)
  iu_prev_unweighted <- sapply(split(y_obs_ord / m_ord, iu_ord), mean)
  iu_zero_frac <- sapply(split(y_obs_ord == 0, iu_ord), mean)
  iu_support <- data.frame(
    iu_id = names(iu_tot_m),
    n_obs = iu_n_obs[names(iu_tot_m)],
    total_m = as.numeric(iu_tot_m),
    total_y = as.numeric(iu_tot_y),
    obs_prev_pooled = as.numeric(iu_tot_y / iu_tot_m),
    obs_prev_unweighted = as.numeric(iu_prev_unweighted[names(iu_tot_m)]),
    zero_obs_frac = as.numeric(iu_zero_frac[names(iu_tot_m)]),
    stringsAsFactors = FALSE
  )
  iu_support$included <- iu_support$n_obs >= 1 & iu_support$total_m >= iu_min_total_m
  iu_keep <- iu_support$iu_id[iu_support$included]
  iu_plot_data_list <- list()

  iu_diagnostics <- list(
    n_iu_total = nrow(iu_support),
    n_iu_included = sum(iu_support$included),
    n_iu_excluded_low_m = sum(iu_support$n_obs >= 1 & iu_support$total_m < iu_min_total_m),
    n_records_zero_y = sum(y_obs_ord == 0),
    n_records_nonzero_y = sum(y_obs_ord > 0)
  )
  iu_details <- list()
  iu_pvals <- data.frame(stat = character(0), ppp_two_sided = numeric(0), stringsAsFactors = FALSE)

  cell_id <- interaction(iu_ord, year_ord, drop = TRUE)
  cell_tot_m <- tapply(m_ord, cell_id, sum)
  cell_keep <- names(cell_tot_m)[cell_tot_m >= iu_year_min_total_m]

  for (est in estimand_levels) {
    if (length(iu_keep) > 0) {
      iu_prev <- function(y_vec) {
        y_tot <- tapply(y_vec, iu_ord, sum)
        if (est == "pooled") {
          y_tot[iu_keep] / iu_tot_m[iu_keep]
        } else {
          vals <- split(y_vec / m_ord, iu_ord)
          sapply(iu_keep, function(k) mean(vals[[k]]))
        }
      }
      obs_iu_mean_est <- iu_prev(y_obs_ord)
      rep_iu_mean_est <- t(apply(yrep_ord, 1, iu_prev))
      if (is.null(dim(rep_iu_mean_est))) rep_iu_mean_est <- matrix(rep_iu_mean_est, nrow = nsim)
      assert_all_finite(obs_iu_mean_est, paste0("IU observed mean prevalence (", est, ")"))
      assert_all_finite(rep_iu_mean_est, paste0("IU replicated mean prevalence (", est, ")"))
      iu_lo_est <- apply(rep_iu_mean_est, 2, stats::quantile, probs = alpha / 2, type = 8)
      iu_md_est <- apply(rep_iu_mean_est, 2, stats::quantile, probs = 0.5, type = 8)
      iu_hi_est <- apply(rep_iu_mean_est, 2, stats::quantile, probs = 1 - alpha / 2, type = 8)
      assert_all_finite(iu_lo_est, paste0("IU replicated lower quantile (", est, ")"))
      assert_all_finite(iu_md_est, paste0("IU replicated median (", est, ")"))
      assert_all_finite(iu_hi_est, paste0("IU replicated upper quantile (", est, ")"))
      T_iu_sse_obs_est <- sum((obs_iu_mean_est - iu_md_est)^2)
      T_iu_sse_rep_est <- apply(rep_iu_mean_est, 1, function(v) sum((v - iu_md_est)^2))
      assert_all_finite(T_iu_sse_rep_est, paste0("IU SSE replicates (", est, ")"))
      p_iu_sse_est <- ppp_two_sided(T_iu_sse_rep_est, T_iu_sse_obs_est)

      ord_iu <- order(obs_iu_mean_est)
      iu_ids_ord <- names(obs_iu_mean_est)[ord_iu]
      support_ord <- iu_support[match(iu_ids_ord, iu_support$iu_id), c("n_obs", "total_m")]
      iu_plot_data_list[[est]] <- data.frame(
        estimand = est,
        iu_id = iu_ids_ord,
        rank = seq_along(iu_ids_ord),
        obs = as.numeric(obs_iu_mean_est[ord_iu]),
        md = as.numeric(iu_md_est[ord_iu]),
        lo = as.numeric(iu_lo_est[ord_iu]),
        hi = as.numeric(iu_hi_est[ord_iu]),
        n_obs = as.integer(support_ord$n_obs),
        total_m = as.numeric(support_ord$total_m),
        stringsAsFactors = FALSE
      )
    } else {
      obs_iu_mean_est <- numeric(0)
      rep_iu_mean_est <- matrix(numeric(0), nrow = nsim, ncol = 0)
      p_iu_sse_est <- NA_real_
      iu_plot_data_list[[est]] <- data.frame(
        estimand = character(0), iu_id = character(0), rank = integer(0),
        obs = numeric(0), md = numeric(0), lo = numeric(0), hi = numeric(0),
        n_obs = integer(0), total_m = numeric(0), stringsAsFactors = FALSE
      )
    }

    if (length(cell_keep) > 0) {
      cell_prev <- function(y_vec) {
        if (est == "pooled") {
          cell_y <- tapply(y_vec, cell_id, sum)
          cell_y[cell_keep] / cell_tot_m[cell_keep]
        } else {
          vals <- split(y_vec / m_ord, cell_id)
          sapply(cell_keep, function(k) mean(vals[[k]]))
        }
      }
      obs_cell_est <- cell_prev(y_obs_ord)
      rep_cell_est <- t(apply(yrep_ord, 1, cell_prev))
      if (is.null(dim(rep_cell_est))) rep_cell_est <- matrix(rep_cell_est, nrow = nsim)
      assert_all_finite(obs_cell_est, paste0("IU-year observed prevalence (", est, ")"))
      assert_all_finite(rep_cell_est, paste0("IU-year replicated prevalence (", est, ")"))
      cell_md_est <- apply(rep_cell_est, 2, stats::quantile, probs = 0.5, type = 8)
      assert_all_finite(cell_md_est, paste0("IU-year replicated median (", est, ")"))
      T_cell_sse_obs_est <- sum((obs_cell_est - cell_md_est)^2)
      T_cell_sse_rep_est <- apply(rep_cell_est, 1, function(v) sum((v - cell_md_est)^2))
      assert_all_finite(T_cell_sse_rep_est, paste0("IU-year SSE replicates (", est, ")"))
      p_cell_sse_est <- ppp_two_sided(T_cell_sse_rep_est, T_cell_sse_obs_est)
    } else {
      obs_cell_est <- numeric(0)
      rep_cell_est <- matrix(numeric(0), nrow = nsim, ncol = 0)
      p_cell_sse_est <- NA_real_
    }

    iu_pvals <- rbind(
      iu_pvals,
      data.frame(
        stat = c(paste0("IU_mean_prev_SSE_", est), paste0("IU_year_cell_prev_SSE_", est)),
        ppp_two_sided = c(p_iu_sse_est, p_cell_sse_est),
        stringsAsFactors = FALSE
      )
    )

    iu_details[[est]] <- list(
      obs_iu_mean = obs_iu_mean_est,
      rep_iu_mean = rep_iu_mean_est,
      obs_cell_prev = obs_cell_est,
      rep_cell_prev = rep_cell_est
    )
  }

  iu_module <- list(
    iu_keep = iu_keep,
    iu_tot_m = stats::setNames(iu_support$total_m, iu_support$iu_id)[iu_keep],
    cell_keep = cell_keep,
    iu_support = iu_support,
    diagnostics = iu_diagnostics,
    details = iu_details,
    pvals = iu_pvals
  )

  # Backward-compatible aliases (pooled if available, else first estimand requested)
  est_primary <- if ("pooled" %in% estimand_levels) "pooled" else estimand_levels[1]
  iu_module$obs_iu_mean <- iu_module$details[[est_primary]]$obs_iu_mean
  iu_module$rep_iu_mean <- iu_module$details[[est_primary]]$rep_iu_mean
  iu_module$obs_cell_prev <- iu_module$details[[est_primary]]$obs_cell_prev
  iu_module$rep_cell_prev <- iu_module$details[[est_primary]]$rep_cell_prev
  iu_module$plot_data_iu_mean_prev <- do.call(rbind, iu_plot_data_list)
  if (nrow(iu_module$plot_data_iu_mean_prev) > 0) {
    iu_module$plot_data_iu_mean_prev$misfit_flag <-
      iu_module$plot_data_iu_mean_prev$obs < iu_module$plot_data_iu_mean_prev$lo |
      iu_module$plot_data_iu_mean_prev$obs > iu_module$plot_data_iu_mean_prev$hi
    misfit_by_est <- tapply(iu_module$plot_data_iu_mean_prev$misfit_flag,
                            iu_module$plot_data_iu_mean_prev$estimand, sum)
    for (est in estimand_levels) {
      nm <- paste0("n_iu_misfit_", est)
      iu_diagnostics[[nm]] <- as.integer(if (est %in% names(misfit_by_est)) misfit_by_est[[est]] else 0)
    }
    iu_diagnostics$n_iu_misfit_total <- as.integer(sum(iu_module$plot_data_iu_mean_prev$misfit_flag))
  } else {
    iu_module$plot_data_iu_mean_prev$misfit_flag <- logical(0)
    for (est in estimand_levels) {
      iu_diagnostics[[paste0("n_iu_misfit_", est)]] <- 0L
    }
    iu_diagnostics$n_iu_misfit_total <- 0L
  }
  iu_module$diagnostics <- iu_diagnostics

  plots <- list()
  plot_note <- NULL
  pct_lab <- function(x) paste0(sub("\\.0$", "", sprintf("%.1f", 100 * x)), "%")
  if (isTRUE(make_plots) && !requireNamespace("bayesplot", quietly = TRUE)) {
    plot_note <- "Plots were not produced because package 'bayesplot' is not installed. Install it to enable plot outputs."
    warning(plot_note, call. = FALSE)
  } else if (isTRUE(make_plots)) {
    plots_global <- list()
    plots_temporal <- list()
    plots_iu <- list()

    S_use <- min(nsim, n_overlay)
    set.seed(seed + 10L)
    idx_draws <- sample(seq_len(nsim), S_use)
    yrep_use <- yrep[idx_draws, , drop = FALSE]
    p_rep_use <- p_rep[idx_draws, , drop = FALSE]

    plots_global$ecdf_overlay_prev <- bayesplot::ppc_ecdf_overlay(y = p_obs, yrep = p_rep_use, discrete = TRUE) +
      scale_x_continuous(labels = pct_lab) +
      labs(x = "Prevalence (%)")
    plots_global$rootogram_counts <- bayesplot::ppc_rootogram(y = y_obs, yrep = yrep_use)
    plots_global$error_hist <- bayesplot::ppc_error_hist(y = y_obs, yrep = yrep_use)
    plots_global$pit_ecdf <- suppressMessages(bayesplot::ppc_loo_pit_ecdf(pit = pit))
    plots_global$ecdf_prev_by_m <- bayesplot::ppc_ecdf_overlay_grouped(y = p_obs, yrep = p_rep_use, group = m_group, discrete = TRUE) +
      scale_x_continuous(labels = pct_lab) +
      labs(x = "Prevalence (%)")
    plots_global$violin_prev_by_m <- bayesplot::ppc_violin_grouped(y = p_obs, yrep = p_rep_use, group = m_group) +
      scale_y_continuous(labels = pct_lab) +
      labs(y = "Prevalence (%)")

    global_stat_meta <- list(
      list(key = "pooled_mean_prev", label = "Pooled mean prevalence", pval_stat = "pooled_mean_prev", pct = TRUE,
           fun = function(z) sum(z) / sum(m)),
      list(key = "unweighted_mean_prev", label = "Unweighted mean prevalence", pval_stat = "unweighted_mean_prev", pct = TRUE,
           fun = function(z) mean(z / m)),
      list(key = "var_prev", label = "Variance of prevalence", pval_stat = "var_prev", pct = TRUE,
           fun = function(z) stats::var(z / m)),
      list(key = "zero_fraction_counts", label = "Zero fraction (counts)", pval_stat = "zero_fraction_counts", pct = TRUE,
           fun = function(z) mean(z == 0)),
      list(key = "zeroish_prev_le_0_01", label = "Zeroish prevalence <= 0.01", pval_stat = "zeroish_prev<=0.01", pct = TRUE,
           fun = function(z) mean((z / m) <= 0.01)),
      list(key = "max_prev", label = "Maximum prevalence", pval_stat = "max_prev", pct = TRUE,
           fun = function(z) max(z / m)),
      list(key = "binom_deviance_plugin_mean", label = "Binomial deviance (plugin mean)", pval_stat = "binom_deviance(plugin_mean)", pct = FALSE,
           fun = function(z) sum(binom_dev(z, m, mu_p))),
      list(key = "max_abs_standardised_resid", label = "Max abs standardised residual", pval_stat = "max_abs_standardised_resid", pct = FALSE,
           fun = function(z) max(abs((z - mu_i) / sqrt(v_i)))),
      list(key = paste0("pi_coverage_", round(100 * prob)), label = sprintf("PI coverage %.0f%%", 100 * prob),
           pval_stat = sprintf("PI_coverage_%d%%", round(100 * prob)), pct = TRUE,
           fun = function(z) mean(z >= q_lo & z <= q_hi))
    )
    plots_global$stat_hist <- list()
    for (md in global_stat_meta) {
      pdat <- bayesplot::ppc_stat_data(y = y_obs, yrep = yrep, stat = md$fun)
      obs_v <- pdat$value[pdat$variable == "y"][1]
      rep_v <- pdat$value[pdat$variable != "y"]
      ppp_v <- global_pvals$ppp_two_sided[match(md$pval_stat, global_pvals$stat)]
      p_stat <- ggplot(data.frame(value = rep_v), aes(x = value)) +
        geom_histogram(bins = 30, fill = "#BFD7EA", color = "white") +
        geom_vline(xintercept = obs_v, linewidth = 1.1, color = "#08306B") +
        theme_bw() +
        labs(
          title = paste0("Global PPC: ", md$label),
          subtitle = sprintf("Observed = %.4g | PPP = %.3f", obs_v, ppp_v),
          x = "Statistic value",
          y = "Frequency"
        )
      if (isTRUE(md$pct)) p_stat <- p_stat + scale_x_continuous(labels = pct_lab) + labs(x = "Statistic (%)")
      plots_global$stat_hist[[md$key]] <- p_stat
    }

    make_line_ribbon_plot <- function(df, xvar, xlabel, ylabel, title_str, facet = FALSE, facet_var = "estimand", x_breaks = NULL, point_size = 1.8) {
      df_lines <- rbind(
        data.frame(df[c(xvar, facet_var)], value = df$obs, series = "Observed"),
        data.frame(df[c(xvar, facet_var)], value = df$md, series = "Simulated median")
      )
      interval_lab <- sprintf("Simulated %.0f%% interval", 100 * prob)
      p <- ggplot(df, aes(x = !!rlang::sym(xvar))) +
        geom_ribbon(aes(ymin = lo, ymax = hi, fill = interval_lab), alpha = 0.2) +
        geom_line(data = df_lines, aes(y = value, color = series, linetype = series), linewidth = 0.8) +
        geom_point(data = df_lines[df_lines$series == "Observed", , drop = FALSE], aes(y = value, color = series), size = point_size) +
        scale_color_manual(values = c("Observed" = "black", "Simulated median" = "#1f78b4"), name = NULL) +
        scale_linetype_manual(values = c("Observed" = "solid", "Simulated median" = "dashed"), name = NULL) +
        scale_fill_manual(values = stats::setNames("#4F93D2", interval_lab), name = NULL) +
        guides(fill = guide_legend(order = 1), color = guide_legend(order = 2), linetype = guide_legend(order = 2)) +
        scale_y_continuous(labels = pct_lab) +
        labs(x = xlabel, y = ylabel, title = title_str) +
        theme_bw()
      if (!is.null(x_breaks)) p <- p + scale_x_continuous(breaks = x_breaks)
      if (isTRUE(facet)) {
        p <- p + facet_wrap(as.formula(paste("~", facet_var)),
                            labeller = as_labeller(c(pooled = "Pooled", unweighted = "Unweighted")))
      }
      p
    }
    make_iu_forest_plot <- function(df, title_str, facet = FALSE, facet_var = "estimand") {
      interval_lab <- sprintf("Simulated %.0f%% interval", 100 * prob)
      p <- ggplot(df, aes(y = rank)) +
        ggplot2::geom_linerange(aes(xmin = lo, xmax = hi, color = interval_lab), linewidth = 0.65) +
        geom_point(aes(x = md, color = "Simulated median"), shape = 124, size = 3.2) +
        geom_point(aes(x = obs, fill = n_obs, color = "Observed"), shape = 21, size = 2.2, stroke = 0.25, alpha = 0.72) +
        scale_color_manual(
          values = c("Observed" = "black", "Simulated median" = "black", stats::setNames("grey20", interval_lab)),
          breaks = c(interval_lab, "Observed", "Simulated median"),
          name = NULL
        ) +
        ggplot2::scale_fill_gradient(name = "Number of surveys (n)", low = "#DCEAF7", high = "#0B4F8A") +
        scale_x_continuous(labels = pct_lab) +
        labs(
          x = "Prevalence (%)",
          y = "IU rank (ordered by observed prevalence)",
          title = title_str
        ) +
        theme_bw() +
        theme(panel.grid.major.y = ggplot2::element_line(color = "grey90"),
              panel.grid.minor.y = ggplot2::element_blank())
      if (isTRUE(facet)) {
        p <- p + facet_wrap(as.formula(paste("~", facet_var)),
                            labeller = as_labeller(c(pooled = "Pooled", unweighted = "Unweighted")))
      }
      p
    }

    if (identical(prevalence_estimand, "both")) {
      year_f <- do.call(rbind, lapply(c("pooled", "unweighted"), function(est) {
        tm <- temporal_summaries[[est]]
        data.frame(estimand = est, year = tm$years, obs = tm$obs_traj, lo = tm$rep_lo, md = tm$rep_md, hi = tm$rep_hi)
      }))
      plots_temporal$year_traj <- make_line_ribbon_plot(
        year_f, "year", "Year", "Prevalence (%)",
        sprintf("Year trend PPC (%.0f%% envelope)", 100 * prob),
        facet = TRUE, x_breaks = yrs
      )

      diff_f <- do.call(rbind, lapply(c("pooled", "unweighted"), function(est) {
        tm <- temporal_summaries[[est]]
        dy <- tm$years[-1]
        lo_d <- apply(tm$d_rep, 2, stats::quantile, probs = alpha / 2, type = 8, na.rm = TRUE)
        md_d <- apply(tm$d_rep, 2, stats::quantile, probs = 0.5, type = 8, na.rm = TRUE)
        hi_d <- apply(tm$d_rep, 2, stats::quantile, probs = 1 - alpha / 2, type = 8, na.rm = TRUE)
        data.frame(estimand = est, year = dy, obs = tm$d_obs, lo = lo_d, md = md_d, hi = hi_d)
      }))
      plots_temporal$year_diffs <- make_line_ribbon_plot(
        diff_f, "year", "Year", "Delta prevalence (%)",
        sprintf("Year-to-year change PPC (%.0f%% envelope)", 100 * prob),
        facet = TRUE, x_breaks = sort(unique(diff_f$year))
      ) + geom_abline(intercept = 0, slope = 0, linetype = 2)

      iu_f <- iu_module$plot_data_iu_mean_prev
      iu_f <- iu_f[iu_f$estimand %in% c("pooled", "unweighted"), , drop = FALSE]
      if (!is.null(iu_f) && nrow(iu_f) > 0) {
        plots_iu$iu_mean_prev <- make_iu_forest_plot(
          iu_f,
          sprintf("IU mean prevalence PPC (%.0f%% envelope)", 100 * prob),
          facet = TRUE
        )
        iu_f_misfit <- iu_f[iu_f$misfit_flag, , drop = FALSE]
        if (nrow(iu_f_misfit) > 0) {
          plots_iu$iu_mean_prev_misfit <- make_iu_forest_plot(
            iu_f_misfit,
            sprintf("IU prevalence PPC: out-of-interval IUs (%.0f%% envelope)", 100 * prob),
            facet = TRUE
          )
        } else {
          plots_iu$iu_mean_prev_misfit <- NULL
        }
      }
    } else {
      est <- estimand_levels[1]
      est_label <- if (est == "pooled") "Pooled prevalence (%)" else "Unweighted prevalence (%)"
      tm <- temporal_summaries[[est]]
      df_year <- data.frame(estimand = est, year = tm$years, obs = tm$obs_traj, lo = tm$rep_lo, md = tm$rep_md, hi = tm$rep_hi)
      plots_temporal$year_traj <- make_line_ribbon_plot(
        df_year, "year", "Year", est_label,
        sprintf("Year trend PPC (%s, %.0f%% envelope)", est, 100 * prob),
        x_breaks = yrs
      )

      dy <- tm$years[-1]
      lo_d <- apply(tm$d_rep, 2, stats::quantile, probs = alpha / 2, type = 8, na.rm = TRUE)
      md_d <- apply(tm$d_rep, 2, stats::quantile, probs = 0.5, type = 8, na.rm = TRUE)
      hi_d <- apply(tm$d_rep, 2, stats::quantile, probs = 1 - alpha / 2, type = 8, na.rm = TRUE)
      df_diff <- data.frame(estimand = est, year = dy, obs = tm$d_obs, lo = lo_d, md = md_d, hi = hi_d)
      plots_temporal$year_diffs <- make_line_ribbon_plot(
        df_diff, "year", "Year", paste0("Delta ", est_label),
        sprintf("Year-to-year change PPC (%s, %.0f%% envelope)", est, 100 * prob),
        x_breaks = dy
      ) + geom_abline(intercept = 0, slope = 0, linetype = 2)

      df_iu <- iu_module$plot_data_iu_mean_prev
      df_iu <- df_iu[df_iu$estimand == est, , drop = FALSE]
      if (nrow(df_iu) > 0) {
        plots_iu$iu_mean_prev <- make_iu_forest_plot(
          df_iu,
          sprintf("IU mean prevalence PPC (%s, %.0f%% envelope)", est, 100 * prob)
        )
        df_iu_misfit <- df_iu[df_iu$misfit_flag, , drop = FALSE]
        if (nrow(df_iu_misfit) > 0) {
          plots_iu$iu_mean_prev_misfit <- make_iu_forest_plot(
            df_iu_misfit,
            sprintf("IU prevalence PPC: out-of-interval IUs (%s, %.0f%% envelope)", est, 100 * prob)
          )
        } else {
          plots_iu$iu_mean_prev_misfit <- NULL
        }
      }
    }

    est_plot_primary <- if ("pooled" %in% estimand_levels) "pooled" else estimand_levels[1]
    rep_cell_primary <- iu_module$details[[est_plot_primary]]$rep_cell_prev
    obs_cell_primary <- iu_module$details[[est_plot_primary]]$obs_cell_prev
    if (!is.null(rep_cell_primary) && ncol(rep_cell_primary) > 0) {
      obs_cell_vec <- as.numeric(obs_cell_primary)
      rep_cell_use <- rep_cell_primary[idx_draws, , drop = FALSE]
      plots_iu$iu_year_cell_ecdf <- bayesplot::ppc_ecdf_overlay(y = obs_cell_vec, yrep = rep_cell_use, discrete = TRUE)
    }

    plots <- list(global = plots_global, temporal = plots_temporal, iu = plots_iu)
  }

  out <- list(
    inputs = list(
      n = n, nsim = nsim, prob = prob,
      m_bins = m_bins, iu_min_total_m = iu_min_total_m,
      iu_year_min_total_m = iu_year_min_total_m,
      prevalence_estimand = prevalence_estimand
    ),
    summaries = list(
      global = list(
        pvals = global_pvals,
        pit_summary = c(mean = mean(pit), var = safe_var(pit)),
        z_summary = stats::quantile(z_obs, probs = c(0.01, 0.05, 0.5, 0.95, 0.99), na.rm = TRUE)
      ),
      temporal = list(
        pvals = temporal_pvals,
        summaries = temporal_summaries,
        years = temporal_summaries[[est_primary]]$years,
        obs_traj = temporal_summaries[[est_primary]]$obs_traj,
        rep_lo = temporal_summaries[[est_primary]]$rep_lo,
        rep_md = temporal_summaries[[est_primary]]$rep_md,
        rep_hi = temporal_summaries[[est_primary]]$rep_hi
      ),
      iu = iu_module
    ),
    plots = plots
  )
  if (!is.null(plot_note)) out$plot_note <- plot_note
  class(out) <- "ppc_dast_dgp"
  out
}

##' Print Method for Internal DAST PPC Objects
##'
##' Prints key posterior predictive discrepancy summaries returned by
##' \code{ppc_dast_dgp()}.
##'
##' @param x An object of class \code{"ppc_dast_dgp"}.
##' @param ... Not used.
##'
##' @return Invisibly returns \code{x}.
##'
##' @method print ppc_dast_dgp
##' @export
##' @keywords internal
print.ppc_dast_dgp <- function(x, ...) {
  cat("DAST DGP PPC summary\n")
  cat(sprintf("n = %d locations; nsim = %d replicates; prob = %.2f\n",
              x$inputs$n, x$inputs$nsim, x$inputs$prob))
  cat("\nGlobal discrepancy PPC p-values:\n")
  print(x$summaries$global$pvals, row.names = FALSE)
  cat("\nTemporal discrepancy PPC p-values:\n")
  print(x$summaries$temporal$pvals, row.names = FALSE)
  cat("\nIU-level discrepancy PPC p-values:\n")
  print(x$summaries$iu$pvals, row.names = FALSE)
  cat("\nInterpretation note:\n")
  cat("PPP values near 0 or 1 indicate potential misfit. Very high PPP values can also\n")
  cat("mean the observed summary is unusually close to the predictive center.\n")
  cat("\nPIT summary:\n")
  print(x$summaries$global$pit_summary)
  cat("\nStandardised residual quantiles:\n")
  print(x$summaries$global$z_summary)
  if (!is.null(x$plot_note)) {
    cat("\nPlot note:\n")
    cat(x$plot_note, "\n")
  }
  invisible(x)
}



dast_fit <-
  function(y, D, coords, time, units_m, kappa, penalty,
           mda_times, survey_times_data,
           sst,
           int_mat,
           par0, cov_offset, power_val,
           ID_coords, ID_re, s_unique, re_unique,
           fix_tau2, fix_alpha, return_samples,
           start_beta,
           start_alpha, start_gamma,
           start_cov_pars,
           start_psi,
           control_mcmc,
           messages = TRUE) {

    beta0 <- par0$beta
    mu0 <- D%*%beta0+cov_offset

    sigma2_0 <- par0$sigma2

    phi0 <- par0$phi

    tau2_0 <- par0$tau2

    if(is.null(tau2_0)) tau2_0 <- fix_tau2

    sigma2_re_0 <- par0$sigma2_re

    n_loc <- nrow(coords)
    n_re <- length(sigma2_re_0)
    n_samples <- (control_mcmc$n_sim-control_mcmc$burnin)/control_mcmc$thin

    u = dist(coords)
    if(sst) v = dist(time)

    if(sst) {
      psi0 <- par0$psi
      Sigma0 <- sigma2_0*matern_cor(u = u, phi = phi0, kappa = kappa,
                                    return_sym_matrix = TRUE) *
        matern_cor(u = v, phi = psi0, kappa = 0.5,
                   return_sym_matrix = TRUE)
    } else {
      Sigma0 <- sigma2_0*matern_cor(u = u, phi = phi0, kappa = kappa,
                                    return_sym_matrix = TRUE)
    }

    diag(Sigma0) <- diag(Sigma0) + tau2_0

    sigma2_re_0 <- par0$sigma2_re

    alpha0 <- par0$alpha
    if(is.null(alpha0)) alpha0<- fix_alpha

    gamma0 <- par0$gamma

    mda_effect0 <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                                      alpha0, gamma0, kappa = power_val)
    if(messages) message("\n - Obtaining covariance matrix and mean for the proposal distribution of the MCMC \n")
    out_maxim <-
      maxim.integrand.dast(y = y, units_m = units_m, Sigma = Sigma0, mu = mu0,
                           mda_effect = mda_effect0,
                           ID_coords = ID_coords, ID_re = ID_re,
                           sigma2_re = sigma2_re_0,
                           hessian = FALSE, gradient = TRUE)

    Sigma_pd <- out_maxim$Sigma.tilde
    mean_pd <- out_maxim$mode

    simulation <-
      Laplace_sampling_MCMC_dast(y = y, units_m = units_m, mu = mu0,
                                 mda_effect = mda_effect0, Sigma = Sigma0,
                                 sigma2_re = sigma2_re_0,
                                 ID_coords = ID_coords, ID_re = ID_re,
                                 control_mcmc = control_mcmc,
                                 Sigma_pd = Sigma_pd, mean_pd = mean_pd,
                                 messages = messages)

    S_tot_samples <- simulation$samples$S

    # Get the number of columns in D
    p <- ncol(D)

    # Define index for beta, sigma2, and phi
    ind_beta <- 1:p
    ind_sigma2 <- p + 1
    ind_phi <- p + 2

    # Conditional indexing based on fix_tau2 and n_re
    if (!is.null(fix_tau2)) {
      if (n_re > 0) {
        ind_sigma2_re <- (p + 2 + 1):(p + 2 + n_re)
        n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[, i])))
      } else {
        n_dim_re <- NULL
        ind_re <- NULL
      }

      if (is.null(fix_alpha)) {
        ind_alpha <- p + n_re + 3
        ind_gamma <- p + n_re + 4
        ind_psi <- p + n_re + 5  # psi when alpha is not fixed
      } else {
        ind_gamma <- p + n_re + 3
        ind_psi <- p + n_re + 4  # psi when alpha is fixed
      }
    } else {
      ind_nu2 <- p + 3
      if (n_re > 0) {
        ind_sigma2_re <- (p + 3 + 1):(p + 3 + n_re)
        n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[, i])))
      } else {
        n_dim_re <- NULL
        ind_re <- NULL
      }

      if (is.null(fix_alpha)) {
        ind_alpha <- p + n_re + 4
        ind_gamma <- p + n_re + 5
        ind_psi <- p + n_re + 6  # psi when alpha is not fixed
      } else {
        ind_gamma <- p + n_re + 4
        ind_psi <- p + n_re + 5  # psi when alpha is fixed
      }
    }


    if(n_re> 0) {
      for(i in 1:n_re) {
        S_tot_samples <- cbind(S_tot_samples, simulation$samples[[i+1]])
      }
      ind_re <- list()
      add_i <- 0
      for(i in 1:n_re) {
        ind_re[[i]] <- (add_i+n_loc+1):(add_i+n_loc+n_dim_re[i])
        if(i < n_re) add_i <- sum(n_dim_re[1:i])
      }
    }


    log.integrand <- function(S_tot, val) {
      n <- length(y)
      S <- S_tot[1:n_loc]

      q.f_re <- 0
      if(n_re > 0) {
        S_re <- NULL
        S_re_list <- list()
        for(i in 1:n_re) {
          S_re_list[[i]] <- S_tot[ind_re[[i]]]
          q.f_re <- q.f_re + n_dim_re[i]*log(val$sigma2_re[i])+
            sum(S_re_list[[i]]^2)/val$sigma2_re[i]
        }
      } else {
        q.f_re <- 0
      }

      eta <- val$mu + S[ID_coords]
      if(n_re > 0) {
        for(i in 1:n_re) {
          eta <- eta + S_re_list[[i]][ID_re[,i]]
        }
      }

      prob_star <- 1/(1+exp(-eta))
      prob <- val$mda_effect*prob_star

      llik <- sum(y*log(prob)+(units_m-y)*log(1-prob))

      q.f_S <- n_loc*log(val$sigma2)+val$ldetR+t(S)%*%val$R.inv%*%S/val$sigma2
      # Subtract both alpha and gamma penalties
      out <- -0.5*(q.f_S+q.f_re)+llik - val$pen_alpha - val$pen_gamma
      return(out)
    }

    compute.log.f <- function(par,ldetR=NA,R.inv=NA) {
      beta <- par[ind_beta]
      sigma2 <- exp(par[ind_sigma2])

      if(is.null(fix_alpha)) {
        alpha <- exp(par[ind_alpha])/(1+exp(par[ind_alpha]))
      } else {
        alpha <- fix_alpha
      }

      gamma <- exp(par[ind_gamma])

      if(length(fix_tau2)>0) {
        nu2 <- fix_tau2/sigma2
      } else {
        nu2 <- exp(par[ind_nu2])
      }
      phi <- exp(par[ind_phi])

      if(sst) psi <- exp(par[ind_psi])

      val <- list()
      val$sigma2 <- sigma2
      val$mu <- as.numeric(D%*%beta)+cov_offset
      val$mda_effect <- compute_mda_effect(survey_times_data, mda_times,
                                           intervention = int_mat,
                                           alpha, gamma, kappa = power_val)
      # Alpha penalty (penalty[[1]]) and gamma penalty (penalty[[4]])
      val$pen_alpha <- penalty[[1]](alpha)
      val$pen_gamma <- penalty[[4]](gamma)
      if(n_re > 0) {
        val$sigma2_re <- exp(par[ind_sigma2_re])
      }
      if(is.na(ldetR) & is.na(as.numeric(R.inv)[1])) {
        if(sst) {
          R <- matern_cor(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)*
            matern_cor(v, phi = psi, kappa=kappa,return_sym_matrix = TRUE)
        } else {
          R <- matern_cor(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)
        }
        diag(R) <- diag(R)+nu2
        val$ldetR <- determinant(R)$modulus
        val$R.inv <- solve(R)
      } else {
        val$ldetR <- ldetR
        val$R.inv <- R.inv
      }
      sapply(1:n_samples,
             function(i) log.integrand(S_tot_samples[i,],val))
    }


    if(!is.null(fix_tau2)) {
      par0_vec <- c(par0$beta, log(c(par0$sigma2, par0$phi)))
    } else {
      par0_vec <- c(par0$beta, log(c(par0$sigma2, par0$phi, par0$tau2/par0$sigma2)))
    }

    if(n_re > 0) {
      if(is.null(fix_alpha)) {
        par0_vec <- c(par0_vec,log(par0$sigma2_re),
                      log(par0$alpha/(1-par0$alpha)), log(par0$gamma))
      } else {
        par0_vec <- c(par0_vec, log(par0$sigma2_re),
                      log(par0$gamma))
      }
    } else {
      if(is.null(fix_alpha)) {
        par0_vec <- c(par0_vec,
                      log(par0$alpha/(1-par0$alpha)), log(par0$gamma))
      } else {
        par0_vec <- c(par0_vec,
                      log(par0$gamma))
      }
    }
    if(sst) par0_vec <- c(par0_vec, log(par0$psi))

    log.f.tilde <- compute.log.f(par0_vec)

    MC.log.lik <- function(par) {
      log(mean(exp(compute.log.f(par)-log.f.tilde)))
    }

    grad.MC.log.lik <- function(par) {
      beta <- par[ind_beta]
      mu <- as.numeric(D %*% beta) + cov_offset
      sigma2 <- exp(par[ind_sigma2])
      if (is.null(fix_alpha)) {
        alpha <- exp(par[ind_alpha]) / (1 + exp(par[ind_alpha]))
      } else {
        alpha <- fix_alpha
      }
      gamma <- exp(par[ind_gamma])
      mda_effect_all <- compute_mda_effect_derivatives(survey_times_data, mda_times,
                                                       intervention = int_mat,
                                                       alpha, gamma, kappa = power_val)
      mda_effect <- mda_effect_all$effect
      mda_der_alpha <- mda_effect_all$d_alpha
      mda_der_gamma <- mda_effect_all$d_gamma

      if (length(fix_tau2) > 0) {
        nu2 <- fix_tau2 / sigma2
      } else {
        nu2 <- exp(par[ind_nu2])
      }
      phi <- exp(par[ind_phi])

      if(sst) psi <- exp(par[ind_psi])
      if (n_re > 0) {
        sigma2_re <- exp(par[ind_sigma2_re])
      }

      if (sst) {
        R_u <- matern_cor(u, phi = phi, kappa = kappa, return_sym_matrix = TRUE)
        R_v <- matern_cor(v, phi = psi, kappa = 0.5, return_sym_matrix = TRUE)
        R <- R_u * R_v
      } else {
        R <- matern_cor(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)
      }
      diag(R) <- diag(R) + nu2
      R.inv <- solve(R)
      ldetR <- determinant(R)$modulus

      if(sst) {
        psi <- exp(par[ind_psi])
        R1.psi <- R_u * matern.grad.phi(v, psi, 0.5)
        m1.psi <- R.inv %*% R1.psi
        t1.psi <- -0.5 * sum(diag(m1.psi))
        m2.psi <- m1.psi %*% R.inv
      }

      exp.fact <- exp(compute.log.f(par, ldetR, R.inv) - log.f.tilde)
      L.m <- sum(exp.fact)
      exp.fact <- exp.fact / L.m

      if (sst) {
        R1.phi <- matern.grad.phi(u, phi, kappa) * R_v
      } else {
        R1.phi <- matern.grad.phi(u, phi, kappa)
      }
      m1.phi <- R.inv %*% R1.phi
      t1.phi <- -0.5 * sum(diag(m1.phi))
      m2.phi <- m1.phi %*% R.inv

      if (is.null(fix_tau2)) {
        t1.nu2 <- -0.5 * sum(diag(R.inv))
        m2.nu2 <- R.inv %*% R.inv
      }

      gradient.S <- function(S_tot) {
        S <- S_tot[1:n_loc]
        if (n_re > 0) {
          S_re_list <- list()
          for (i in 1:n_re) {
            S_re_list[[i]] <- S_tot[ind_re[[i]]]
          }
        }

        eta <- mu + S[ID_coords]
        if (n_re > 0) {
          for (i in 1:n_re) {
            eta <- eta + S_re_list[[i]][ID_re[, i]]
          }
        }

        prob_star <- 1 / (1 + exp(-eta))
        prob <- mda_effect * prob_star
        d_S <- (y / prob - (units_m - y) / (1 - prob)) * mda_effect * prob_star / (1 + exp(eta))
        d_S <- as.numeric(d_S)

        q.f_S <- t(S) %*% R.inv %*% S
        grad.beta <- t(D) %*% d_S
        grad.log.sigma2 <- (-n_loc / (2 * sigma2) + 0.5 * q.f_S / (sigma2^2)) * sigma2
        grad.log.phi <- (t1.phi + 0.5 * as.numeric(t(S) %*% m2.phi %*% S) / sigma2) * phi

        if (is.null(fix_alpha)) {
          der.alpha <- exp(par[ind_alpha]) / ((1 + exp(par[ind_alpha]))^2)
          grad.alpha.t <- der.alpha * (
            sum((y / mda_effect - (units_m - y) * prob_star / (1 - prob)) * mda_der_alpha) -
              penalty[[2]](alpha)
          )
        }

        # Gradient wrt log(gamma): chain rule gives factor gamma.
        # Also subtract the derivative of the gamma penalty wrt log(gamma),
        # which is gamma * penalty[[5]](gamma) by the chain rule.
        grad.log.gamma <- gamma * (
          sum((y / mda_effect - (units_m - y) * prob_star / (1 - prob)) * mda_der_gamma) -
            penalty[[5]](gamma)
        )

        out <- c(grad.beta, grad.log.sigma2, grad.log.phi)

        if (is.null(fix_tau2)) {
          grad.log.nu2 <- (t1.nu2 + 0.5 * as.numeric(t(S) %*% m2.nu2 %*% S) / sigma2) * nu2
          out <- c(out, grad.log.nu2)
        }

        if (n_re > 0) {
          grad.log.sigma2_re <- rep(NA, n_re)
          for (i in 1:n_re) {
            grad.log.sigma2_re[i] <- (-n_dim_re[i] / (2 * sigma2_re[i]) +
                                        0.5 * sum(S_re_list[[i]]^2) / (sigma2_re[i]^2)) * sigma2_re[i]
          }
          out <- c(out, grad.log.sigma2_re)
        }

        if (is.null(fix_alpha)) {
          out <- c(out, grad.alpha.t, grad.log.gamma)
        } else {
          out <- c(out, grad.log.gamma)
        }

        if (sst) {
          grad.log.psi <- (t1.psi + 0.5 * as.numeric(t(S) %*% m2.psi %*% S) / sigma2) * psi
          out <- c(out, grad.log.psi)
        }

        return(out)
      }

      out <- rep(0, length(par))
      for (i in 1:n_samples) {
        out <- out + exp.fact[i] * gradient.S(S_tot_samples[i, ])
      }
      out
    }

    hess.MC.log.lik <- function(par) {
      # Unpack parameters
      beta   <- par[ind_beta]
      mu     <- as.numeric(D %*% beta) + cov_offset
      sigma2 <- exp(par[ind_sigma2])
      if (is.null(fix_alpha)) {
        alpha <- exp(par[ind_alpha]) / (1 + exp(par[ind_alpha]))
      } else {
        alpha <- fix_alpha
      }
      gamma <- exp(par[ind_gamma])

      # MDA effect + derivatives
      mda_all           <- compute_mda_effect_derivatives(
        survey_times_data, mda_times,
        intervention = int_mat,
        alpha, gamma, kappa = power_val
      )
      mda_eff        <- mda_all$effect
      mda_der_alpha  <- mda_all$d_alpha
      mda_der_gamma  <- mda_all$d_gamma
      mda_der2_alpha <- mda_all$d2_alpha
      mda_der2_gamma <- mda_all$d2_gamma
      mda_der2_alpha_gamma <- mda_all$d2_alpha_gamma

      # nu2
      if (!is.null(fix_tau2)) {
        nu2 <- fix_tau2 / sigma2
      } else {
        nu2 <- exp(par[ind_nu2])
      }

      # Build R and its derivatives
      phi <- exp(par[ind_phi])
      if (sst) {
        psi    <- exp(par[ind_psi])
        R_u    <- matern_cor(u,   phi = phi, kappa = kappa, return_sym_matrix = TRUE)
        R_v    <- matern_cor(v,   phi = psi, kappa = 0.5,    return_sym_matrix = TRUE)
        R      <- R_u * R_v
        R1.phi <- matern.grad.phi(u, phi, kappa) * R_v
        R2.phi <- matern.hessian.phi(u, phi, kappa) * R_v
        R1.psi <- R_u * matern.grad.phi(v, psi, 0.5)
        R2.psi <- R_u * matern.hessian.phi(v, psi, 0.5)
      } else {
        R      <- matern_cor(u, phi = phi, kappa = kappa, return_sym_matrix = TRUE)
        R1.phi <- matern.grad.phi(u, phi, kappa)
        R2.phi <- matern.hessian.phi(u, phi, kappa)
      }
      diag(R) <- diag(R) + nu2
      R.inv   <- solve(R)

      # Precompute phi derivatives
      m1.phi <- R.inv %*% R1.phi
      t1.phi <- -0.5 * sum(diag(m1.phi))
      m2.phi <- m1.phi %*% R.inv
      t2.phi <- -0.5 * sum(diag(R.inv %*% R2.phi - R.inv %*% R1.phi %*% R.inv %*% R1.phi))
      n2.phi <- R.inv %*% (2 * R1.phi %*% R.inv %*% R1.phi - R2.phi) %*% R.inv

      if (sst) {
        R2.psi.phi <- matern.grad.phi(u, phi, 0.5)*matern.grad.phi(v, psi, 0.5)
        m1.psi   <- R.inv %*% R1.psi
        t1.psi   <- -0.5 * sum(diag(m1.psi))
        m2.psi   <- m1.psi %*% R.inv
        t2.psi   <- -0.5 * sum(diag(R.inv %*% R2.psi - R.inv %*% R1.psi %*% R.inv %*% R1.psi))
        n2.psi   <- R.inv %*% (2 * R1.psi %*% R.inv %*% R1.psi - R2.psi) %*% R.inv

        t2.psi.phi <- -0.5*(sum(R.inv*R2.psi.phi)-
                              sum(m1.phi*t(m1.psi)))
        n2.psi.phi <- R.inv%*%(R1.phi%*%m1.psi+
                                 R1.psi%*%m1.phi-
                                 R2.psi.phi)%*%R.inv

      }

      # Precompute nu2 derivatives
      if (is.null(fix_tau2)) {
        m2.nu2     <- R.inv %*% R.inv
        t1.nu2     <- -0.5 * sum(diag(R.inv))
        t2.nu2     <-  0.5 * sum(diag(m2.nu2))
        n2.nu2     <-  2 * R.inv %*% m2.nu2
        t2.nu2.phi <-  0.5 * sum(diag(R.inv %*% R1.phi %*% R.inv))
        n2.nu2.phi <-  R.inv %*% (R.inv %*% R1.phi + R1.phi %*% R.inv) %*% R.inv
        if (sst) {
          t2.nu2.psi <-  0.5 * sum(diag(R.inv %*% R1.psi %*% R.inv))
          n2.nu2.psi <-  R.inv %*% (R.inv %*% R1.psi + R1.psi %*% R.inv) %*% R.inv
        }
      }

      # Monte Carlo weights
      ldetR <- determinant(R)$modulus
      expf  <- exp(compute.log.f(par, ldetR, R.inv) - log.f.tilde)
      w     <- expf / sum(expf)

      # Accumulators
      H_acc <- matrix(0, nrow = length(par), ncol = length(par))
      g_acc <- rep(0, length(par))

      hessian.S <- function(S_tot, ef) {
        S <- S_tot[1:n_loc]
        if (n_re > 0) {
          S_re_list <- lapply(seq_len(n_re), function(i) S_tot[ind_re[[i]]])
        }

        # Linear predictor
        eta <- mu + S[ID_coords]
        if (n_re > 0) for (i in seq_len(n_re)) eta <- eta + S_re_list[[i]][ID_re[, i]]
        prob_star <- 1 / (1 + exp(-eta))
        prob      <- mda_eff * prob_star

        # Derivatives wrt S
        d_S  <- (y/prob - (units_m-y)/(1-prob)) * mda_eff * prob_star / (1+exp(eta))
        d2_S <- (-y/prob^2 - (units_m-y)/(1-prob)^2) * (mda_eff*prob_star/(1+exp(eta)))^2 +
          (y/prob - (units_m-y)/(1-prob)) * mda_eff * exp(eta)*(1-exp(eta))/(1+exp(eta))^3
        d_S  <- as.numeric(d_S)

        # d_S_alpha
        if (is.null(fix_alpha)) {
          d_S_alpha <- (-y/prob^2 - (units_m-y)/(1-prob)^2) *
            mda_der_alpha * prob_star * mda_eff * prob_star/(1+exp(eta)) +
            (y/prob - (units_m-y)/(1-prob)) * mda_der_alpha * prob_star/(1+exp(eta))
        }

        # d_S_gamma
        d_S_gamma <- (-y/prob^2 - (units_m-y)/(1-prob)^2) *
          mda_der_gamma * prob_star * mda_eff * prob_star/(1+exp(eta)) +
          (y/prob - (units_m-y)/(1-prob)) * mda_der_gamma * prob_star/(1+exp(eta))

        qSS <- as.numeric(t(S) %*% R.inv %*% S)

        # Gradient components
        grad_beta  <- t(D) %*% d_S
        grad_ls2   <- (-n_loc/(2*sigma2) + 0.5*qSS/sigma2^2) * sigma2
        grad_lphi  <- (t1.phi + 0.5*as.numeric(t(S)%*%m2.phi%*%S)/sigma2) * phi

        g_vec <- c(grad_beta, grad_ls2, grad_lphi)

        if (is.null(fix_tau2)) {
          grad_lnu2 <- (t1.nu2 + 0.5*as.numeric(t(S)%*%m2.nu2%*%S)/sigma2)*nu2
          g_vec     <- c(g_vec, grad_lnu2)
        }

        if (is.null(fix_alpha)) {
          der_alpha  <- exp(par[ind_alpha])/(1+exp(par[ind_alpha]))^2
          grad_alpha <- der_alpha * (
            sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_alpha) -
              penalty[[2]](alpha)
          )
          g_vec <- c(g_vec, grad_alpha)
        }

        # Gradient wrt log(gamma), including gamma penalty derivative.
        # d/d(log gamma) = gamma * d/d(gamma), so subtract gamma * penalty[[5]](gamma).
        grad_lgamma <- gamma * (
          sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_gamma) -
            penalty[[5]](gamma)
        )
        g_vec <- c(g_vec, grad_lgamma)

        if (n_re > 0) {
          grad_ls2re <- sapply(seq_len(n_re), function(i)
            (-n_dim_re[i]/(2*sigma2_re[i]) + 0.5*sum(S_re_list[[i]]^2)/sigma2_re[i]^2)*sigma2_re[i]
          )
          g_vec <- c(g_vec, grad_ls2re)
        }

        if (sst) {
          grad_lpsi <- (t1.psi + 0.5*as.numeric(t(S)%*%m2.psi%*%S)/sigma2)*psi
          g_vec     <- c(g_vec, grad_lpsi)
        }

        # Build H_loc
        H_loc <- matrix(0, nrow=length(par), ncol=length(par))

        # beta-beta
        H_loc[ind_beta, ind_beta] <- t(D) %*% (D * d2_S)

        # beta-alpha
        if (is.null(fix_alpha)) {
          H_loc[ind_beta, ind_alpha] <- t(D) %*% d_S_alpha * der_alpha
          H_loc[ind_alpha, ind_beta] <- H_loc[ind_beta, ind_alpha]
        }

        # beta-gamma
        H_loc[ind_beta, ind_gamma] <- t(D) %*% d_S_gamma * gamma
        H_loc[ind_gamma, ind_beta] <- H_loc[ind_beta, ind_gamma]

        # alpha-alpha
        if (is.null(fix_alpha)) {
          der2_alpha <- exp(par[ind_alpha])*(1-exp(par[ind_alpha]))/(1+exp(par[ind_alpha]))^3
          H_loc[ind_alpha, ind_alpha] <- der2_alpha * (
            sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_alpha) -
              penalty[[2]](alpha)
          ) + der_alpha^2 * (
            sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der2_alpha +
                  (-y/mda_eff^2 - (units_m-y)*(prob_star^2)/(1-prob)^2)*mda_der_alpha^2
            ) - penalty[[3]](alpha)
          )
        }

        # alpha-gamma
        if (is.null(fix_alpha)) {
          H_loc[ind_alpha, ind_gamma] <- (gamma*der_alpha)*sum(
            (y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der2_alpha_gamma +
              (-y/mda_eff^2 - (units_m-y)*(prob_star^2)/(1-prob)^2)*mda_der_alpha*mda_der_gamma
          )
          H_loc[ind_gamma, ind_alpha] <- H_loc[ind_alpha, ind_gamma]
        }

        # gamma-gamma.
        # The Hessian of the penalized log-likelihood wrt log(gamma) is:
        #   H_lik[gamma,gamma]
        #   - gamma^2 * penalty[[6]](gamma)   [from d/d(log g) of gamma * pen'(g)]
        #   - gamma   * penalty[[5]](gamma)   [the grad_lgamma term itself]
        H_loc[ind_gamma, ind_gamma] <- gamma*sum(
          (y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_gamma
        ) + gamma^2*sum(
          (y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der2_gamma +
            (-y/mda_eff^2 - (units_m-y)*(prob_star^2)/(1-prob)^2)*mda_der_gamma^2
        ) - gamma^2 * penalty[[6]](gamma) - gamma * penalty[[5]](gamma)

        # sigma2-sigma2
        H_loc[ind_sigma2, ind_sigma2] <-
          (n_loc/(2*sigma2^2) - qSS/sigma2^3)*sigma2^2 + grad_ls2

        # sigma2-phi
        H_loc[ind_sigma2, ind_phi] <- (grad_lphi/phi - t1.phi)*(-phi)
        H_loc[ind_phi, ind_sigma2] <- H_loc[ind_sigma2, ind_phi]

        # phi-phi
        H_loc[ind_phi, ind_phi] <- (t2.phi - 0.5*as.numeric(t(S)%*%n2.phi%*%S)/sigma2)*phi^2 + grad_lphi

        # nu2-block
        if (is.null(fix_tau2)) {
          H_loc[ind_sigma2, ind_nu2] <- (grad_lnu2/nu2 - t1.nu2)*(-nu2)
          H_loc[ind_nu2, ind_sigma2] <- H_loc[ind_sigma2, ind_nu2]
          H_loc[ind_nu2, ind_nu2]   <- (t2.nu2 - 0.5*as.numeric(t(S)%*%n2.nu2%*%S)/sigma2)*nu2^2 + grad_lnu2
          H_loc[ind_phi, ind_nu2]   <- (t2.nu2.phi - 0.5*as.numeric(t(S)%*%n2.nu2.phi%*%S)/sigma2)*phi*nu2
          H_loc[ind_nu2, ind_phi]   <- H_loc[ind_phi, ind_nu2]
        }

        # psi-block
        if (sst) {
          H_loc[ind_psi, ind_psi]     <- (t2.psi - 0.5*as.numeric(t(S)%*%n2.psi%*%S)/sigma2)*psi^2 + grad_lpsi
          H_loc[ind_sigma2, ind_psi]    <- (grad_lpsi/psi - t1.psi)*(-psi)
          H_loc[ind_psi, ind_sigma2]    <- H_loc[ind_sigma2, ind_psi]
          H_loc[ind_phi, ind_psi]       <- (t2.psi.phi-0.5*t(S)%*%n2.psi.phi%*%S/sigma2)*phi*psi
          H_loc[ind_psi, ind_phi]       <- H_loc[ind_phi, ind_psi]
          if (is.null(fix_tau2)) {
            val <- (t2.nu2.psi - 0.5*as.numeric(t(S)%*%n2.nu2.psi%*%S)/sigma2)*nu2*psi
            H_loc[ind_nu2, ind_psi] <- val
            H_loc[ind_psi, ind_nu2] <- val
          }
        }

        list(
          mat1 = ef * (g_vec %*% t(g_vec) + H_loc),
          g    = ef * g_vec
        )
      }

      # Monte Carlo aggregation
      for (i in seq_len(n_samples)) {
        tmp   <- hessian.S(S_tot_samples[i, ], w[i])
        g_acc <- g_acc + tmp$g
        H_acc <- H_acc + tmp$mat1
      }

      H_acc - g_acc %*% t(g_acc)
    }

    start_cov_pars[-(1:2)] <- start_cov_pars[-(1:2)]/start_cov_pars[1]
    if(is.null(fix_alpha)) {
      start_par <- c(start_beta, log(start_cov_pars),
                     log(start_alpha/(1-start_alpha)), log(start_gamma))
    } else {
      start_par <- c(start_beta, log(start_cov_pars), log(start_gamma))
    }
    if(sst) start_par <- c(start_par, log(start_psi))


    out <- list()
    estim <- nlminb(
      start     = start_par,
      objective = function(x)   -MC.log.lik(x),
      gradient  = function(x)   -grad.MC.log.lik(x),
      hessian   = function(x)   -hess.MC.log.lik(x),
      control   = list(trace=1*messages)
    )

    out$estimate <- estim$par
    out$grad.MLE <- grad.MC.log.lik(estim$par)
    hess.MLE <- hess.MC.log.lik(estim$par)
    out$covariance <- solve(-hess.MLE)
    out$log.lik <- -estim$objective
    if(return_samples) out$S_samples <- S_tot_samples
    class(out) <- "RiskMap"
    return(out)
  }


maxim.integrand.dast <- function(y,units_m,mu,Sigma,ID_coords, ID_re = NULL,
                                 mda_effect,
                                 sigma2_re = NULL,
                                 hessian=FALSE, gradient=FALSE) {
  # Sigma <- Sigma0
  # mda_effect <- mda_effect0
  # mu <- mu0
  Sigma.inv <- solve(Sigma)
  n_loc <- nrow(Sigma)
  n <- length(y)
  if((!is.null(ID_re) & is.null(sigma2_re)) | (is.null(ID_re) & !is.null(sigma2_re))) {
    stop("To introduce unstructured random effects both `ID_re` and `sigma2_re`
           must be provided.")
  }
  n_re <- length(sigma2_re)
  if(n_re > 0) {
    n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[,i])))
    ind_re <- list()
    add_i <- 0
    for(i in 1:n_re) {
      ind_re[[i]] <- (add_i+n_loc+1):(add_i+n_loc+n_dim_re[i])
      if(i < n_re) add_i <- sum(n_dim_re[1:i])
    }
  }
  n_tot <- n_loc
  if(n_re > 0) n_tot <- n_tot + sum(n_dim_re)

  integrand <- function(S_tot) {
    S <- S_tot[1:n_loc]

    q.f_S <- as.numeric(t(S)%*%Sigma.inv%*%(S))

    q.f_re <- 0
    if(n_re > 0) {
      S_re <- NULL
      S_re_list <- list()
      for(i in 1:n_re) {
        S_re_list[[i]] <- S_tot[ind_re[[i]]]
        q.f_re <- q.f_re + sum(S_re_list[[i]]^2)/sigma2_re[i]
      }
    }

    eta <- mu + S[ID_coords]
    if(n_re > 0) {
      for(i in 1:n_re) {
        eta <- eta + S_re_list[[i]][ID_re[,i]]
      }
    }
    prob_star <- 1/(1+exp(-eta))
    prob <- mda_effect*prob_star

    llik <- sum(y*log(prob)+(units_m-y)*log(1-prob))


    out <- -0.5*q.f_S-0.5*q.f_re+llik
    return(out)
  }


  C_S <- t(sapply(1:n_loc,function(i) ID_coords==i))

  if(n_re>0) {
    C_re <- list()
    C_S_re <- list()
    C_re_re <- list()
    for(j in 1:n_re){
      C_S_re[[j]] <- array(FALSE,dim = c(n_loc, n_dim_re[j], n))
      C_re[[j]] <- t(sapply(1:n_dim_re[j],function(i) ID_re[,j]==i))
      for(l in 1:n_dim_re[j]) {
        for(k in 1:n_loc) {
          ind_kl <- which(ID_coords==k & ID_re[,j]==l)
          if(length(ind_kl) > 0) {
            C_S_re[[j]][k,l,ind_kl] <- TRUE
          }
        }
      }

      if(j < n_re) {
        C_re_re[[j]] <- list()
        counter <- 0
        for(w in (j+1):n_re) {
          counter <- counter+1
          C_re_re[[j]][[counter]] <- array(FALSE,dim = c(n_dim_re[j], n_dim_re[w], n))
          for(l in 1:n_dim_re[j]) {
            for(k in 1:n_dim_re[w]) {
              ind_lk <- which(ID_re[,j]==l & ID_re[,w]==k)
              if(length(ind_kl) > 0) {
                C_re_re[[j]][[counter]][l,k,ind_lk] <- TRUE
              }
            }
          }
        }
      }
    }
  }

  grad.integrand <- function(S_tot) {
    S <- S_tot[1:n_loc]

    if(n_re > 0) {
      S_re_list <- list()
      for(i in 1:n_re) {
        S_re_list[[i]] <- S_tot[ind_re[[i]]]
      }
    }

    eta <- mu + S[ID_coords]
    if(n_re > 0) {
      for(i in 1:n_re) {
        eta <- eta + S_re_list[[i]][ID_re[,i]]
      }
    }

    prob_star <- 1 / (1 + exp(-eta))
    prob <- mda_effect * prob_star

    # Compute derivative of log-likelihood with respect to eta
    d_S <- (y/prob - (units_m-y)/(1-prob))*mda_effect*prob_star/(1+exp(eta))
    d_S <- as.numeric(d_S)

    out <- rep(NA,n_tot)
    out[1:n_loc] <- as.numeric(-Sigma.inv%*%S+
                                 sapply(1:n_loc,function(i) sum(d_S[C_S[i,]])))
    if(n_re>0) {
      for(j in 1:n_re) {
        out[ind_re[[j]]] <- as.numeric(-S_re_list[[j]]/sigma2_re[[j]]+
                                         sapply(1:n_dim_re[[j]],
                                                function(x) sum(d_S[C_re[[j]][x,]])))
      }
    }
    return(out)
  }


  hessian.integrand <- function(S_tot) {
    S <- S_tot[1:n_loc]


    if(n_re > 0) {
      S_re <- NULL
      S_re_list <- list()
      for(i in 1:n_re) {
        S_re_list[[i]] <- S_tot[ind_re[[i]]]
      }
    }

    eta <- mu + S[ID_coords]
    if(n_re > 0) {
      for(i in 1:n_re) {
        eta <- eta + S_re_list[[i]][ID_re[,i]]
      }
    }

    prob_star <- 1 / (1 + exp(-eta))
    prob <- mda_effect * prob_star

    # Compute derivative of log-likelihood with respect to eta
    d2_S <- (-y/prob^2 - (units_m-y)/((1-prob)^2))*(mda_effect*prob_star/(1+exp(eta)))^2+
      (y/prob - (units_m-y)/(1-prob))*mda_effect*exp(eta)*(1-exp(2*eta))/((1+exp(eta))^4)
    d2_S <- as.numeric(d2_S)
    out <- matrix(0,nrow = n_tot, ncol = n_tot)

    out[1:n_loc, 1:n_loc] <-  -Sigma.inv
    diag(out[1:n_loc, 1:n_loc]) <- diag(out[1:n_loc, 1:n_loc])+
      sapply(1:n_loc,function(i) sum(d2_S[C_S[i,]]))
    if(n_re>0) {
      for(j in 1:n_re) {
        diag(out[ind_re[[j]], ind_re[[j]]]) <- -1/sigma2_re[j]
        diag(out[ind_re[[j]], ind_re[[j]]]) <- diag(out[ind_re[[j]], ind_re[[j]]])+
          sapply(1:n_dim_re[j],function(i) sum(d2_S[C_re[[j]][i,]]))

        out[1:n_loc,ind_re[[j]]]

        for(k in 1:n_dim_re[[j]]) {
          out[1:n_loc, ind_re[[j]]][,k] <- -sapply(1:n_loc,function(i) sum(d2_S[C_S_re[[j]][i,k,]]))
          out[ind_re[[j]], 1:n_loc][k,] <- out[1:n_loc,ind_re[[j]]][,k]
        }

        if(j < n_re) {
          counter <- 0
          for(w in (j+1):n_re) {
            counter <- counter + 1
            for(k in 1:n_dim_re[[w]]) {
              out[ind_re[[j]], ind_re[[w]]][,k] <- -sapply(1:n_dim_re[j],function(i) sum(d2_S[C_re_re[[j]][[counter]][i,k,]]))
              out[ind_re[[w]], ind_re[[j]]][k,] <- out[ind_re[[j]], ind_re[[w]]][,k]
            }
          }
        }
      }
    }
    return(out)
  }

  estim <- nlminb(rep(0,n_tot),
                  function(x) -integrand(x),
                  function(x) -grad.integrand(x),
                  function(x) -hessian.integrand(x))


  out <- list()
  out$mode <- estim$par
  if(hessian) {
    out$hessian <- hessian.integrand(out$mode)
  } else {
    out$Sigma.tilde <- solve(-hessian.integrand(out$mode))
  }

  if(gradient) {
    out$gradient <- grad.integrand(out$mode)
  }

  return(out)
}
Laplace_sampling_MCMC_dast <- function(y, units_m, mu, mda_effect, Sigma, ID_coords, ID_re = NULL,
                                       sigma2_re = NULL, control_mcmc,
                                       Sigma_pd = NULL, mean_pd = NULL, messages = TRUE) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Precompute values that will not change
  Sigma.inv <- solve(Sigma)
  n_loc <- nrow(Sigma)
  n <- length(y)

  if ((!is.null(ID_re) & is.null(sigma2_re)) | (is.null(ID_re) & !is.null(sigma2_re))) {
    stop("To introduce unstructured random effects both `ID_re` and `sigma2_re` must be provided.")
  }

  n_re <- length(sigma2_re)
  n_tot <- if (n_re > 0) n_loc + sum(sapply(1:n_re, function(i) length(unique(ID_re[, i])))) else n_loc

  # Precompute the matrices
  C_S <- t(sapply(1:n_loc, function(i) ID_coords == i))

  if (n_re > 0) {
    ind_re <- vector("list", n_re)
    add_i <- 0
    n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[, i])))

    for (i in 1:n_re) {
      ind_re[[i]] <- (add_i + n_loc + 1):(add_i + n_loc + n_dim_re[i])
      if (i < n_re) add_i <- sum(n_dim_re[1:i])
    }

    C_re <- vector("list", n_re)
    C_S_re <- vector("list", n_re)
    C_re_re <- vector("list", n_re)

    for (j in 1:n_re) {
      C_S_re[[j]] <- array(FALSE, dim = c(n_loc, n_dim_re[j], n))
      C_re[[j]] <- t(sapply(1:n_dim_re[j], function(i) ID_re[, j] == i))
      for (l in 1:n_dim_re[j]) {
        for (k in 1:n_loc) {
          ind_kl <- which(ID_coords == k & ID_re[, j] == l)
          if (length(ind_kl) > 0) C_S_re[[j]][k, l, ind_kl] <- TRUE
        }
      }

      if (j < n_re) {
        C_re_re[[j]] <- vector("list", n_re - j)
        counter <- 0
        for (w in (j + 1):n_re) {
          counter <- counter + 1
          C_re_re[[j]][[counter]] <- array(FALSE, dim = c(n_dim_re[j], n_dim_re[w], n))
          for (l in 1:n_dim_re[j]) {
            for (k in 1:n_dim_re[w]) {
              ind_lk <- which(ID_re[, j] == l & ID_re[, w] == k)
              if (length(ind_lk) > 0) {
                C_re_re[[j]][[counter]][l, k, ind_lk] <- TRUE
              }
            }
          }
        }
      }
    }
  }

  if (is.null(Sigma_pd) | is.null(mean_pd)) {
    out_maxim <- maxim.integrand.dast(y = y, units_m = units_m, Sigma = Sigma, mu = mu,
                                      mda_effect = mda_effect, ID_coords = ID_coords,
                                      ID_re = ID_re, sigma2_re = sigma2_re,
                                      hessian = FALSE, gradient = TRUE)

    if (is.null(Sigma_pd)) Sigma_pd <- out_maxim$Sigma.tilde
    if (is.null(mean_pd)) mean_pd <- out_maxim$mode
  }

  n_sim  <- control_mcmc$n_sim
  burnin <- control_mcmc$burnin
  thin   <- control_mcmc$thin
  h      <- control_mcmc$h %||% (1.65 / (n_tot^(1 / 6)))  # Default if not provided
  c1.h   <- control_mcmc$c1.h
  c2.h   <- control_mcmc$c2.h

  Sigma_pd_sroot <- t(chol(Sigma_pd))
  A <- solve(Sigma_pd_sroot)

  # Construct total covariance matrix (Sigma_tot)
  Sigma_tot <- if (n_re == 0) Sigma else {
    Sigma_tot <- matrix(0, n_tot, n_tot)
    Sigma_tot[1:n_loc, 1:n_loc] <- Sigma
    for (i in 1:n_re) diag(Sigma_tot)[ind_re[[i]]] <- sigma2_re[i]
    Sigma_tot
  }

  Sigma_w_inv <- solve(A %*% Sigma_tot %*% t(A))
  mu_w <- -as.numeric(A %*% mean_pd)

  cond.dens.W <- function(W, S_tot) {
    S <- S_tot[1:n_loc]

    eta <- mu + S[ID_coords]
    if (n_re > 0) {
      for (i in 1:n_re) {
        eta <- eta + S_tot[ind_re[[i]]][ID_re[, i]]
      }
    }

    prob_star <- stats::plogis(eta)             # in (0,1)
    prob <- mda_effect * prob_star              # expected in (0,1)
    eps <- .Machine$double.eps
    prob <- pmin(pmax(prob, eps), 1 - eps)      # clamp for numerical stability

    llik <- sum(y * log(prob) + (units_m - y) * log(1 - prob))
    diff_w <- W - mu_w

    as.numeric(-0.5 * crossprod(diff_w, Sigma_w_inv %*% diff_w) + llik)
  }

  lang.grad <- function(W, S_tot) {
    S <- S_tot[1:n_loc]

    eta <- mu + S[ID_coords]
    if (n_re > 0) {
      for (i in 1:n_re) {
        eta <- eta + S_tot[ind_re[[i]]][ID_re[, i]]
      }
    }

    prob_star <- stats::plogis(eta)
    prob <- mda_effect * prob_star
    eps <- .Machine$double.eps
    prob <- pmin(pmax(prob, eps), 1 - eps)

    # derivative wrt eta
    # d/deta prob = mda_effect * prob_star * (1 - prob_star)
    dprob_deta <- mda_effect * prob_star * (1 - prob_star)
    d_llik_deta <- (y / prob - (units_m - y) / (1 - prob)) * dprob_deta

    grad_S_tot_r <- rep(NA_real_, n_tot)
    # contributions to S (size n_loc)
    grad_S_tot_r[1:n_loc] <- as.numeric(sapply(1:n_loc, function(i) sum(d_llik_deta[C_S[i, ]])))

    if (n_re > 0) {
      for (j in 1:n_re) {
        # random effect prior term and likelihood term
        grad_like_re <- sapply(1:n_dim_re[[j]], function(x) sum(d_llik_deta[C_re[[j]][x, ]]))
        grad_S_tot_r[ind_re[[j]]] <- as.numeric(-S_tot[ind_re[[j]]] / sigma2_re[[j]] + grad_like_re)
      }
    }

    as.numeric(-Sigma_w_inv %*% (W - mu_w) + t(Sigma_pd_sroot) %*% grad_S_tot_r)
  }

  # MCMC loop
  W_curr <- rep(0, n_tot)
  S_tot_curr <- as.numeric(Sigma_pd_sroot %*% W_curr + mean_pd)
  mean_curr <- as.numeric(W_curr + (h^2 / 2) * lang.grad(W_curr, S_tot_curr))
  lp_curr <- cond.dens.W(W_curr, S_tot_curr)

  acc <- 0L
  n_samples <- floor((n_sim - burnin) / thin)
  sim <- matrix(NA_real_, nrow = n_samples, ncol = n_tot)

  if (messages) message("\n - Conditional simulation (burnin=", burnin, ", thin=", thin, "):")
  pb <- NULL
  if (messages && interactive()) {
    pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  h.vec <- rep(NA_real_, n_sim)
  acc_prob <- rep(NA_real_, n_sim)

  for (i in 1:n_sim) {
    W_prop <- mean_curr + h * rnorm(n_tot)
    S_tot_prop <- as.numeric(Sigma_pd_sroot %*% W_prop + mean_pd)
    mean_prop <- as.numeric(W_prop + (h^2 / 2) * lang.grad(W_prop, S_tot_prop))
    lp_prop <- cond.dens.W(W_prop, S_tot_prop)

    dprop_curr <- -sum((W_prop - mean_curr)^2) / (2 * (h^2))
    dprop_prop <- -sum((W_curr - mean_prop)^2) / (2 * (h^2))

    log_prob <- lp_prop + dprop_prop - lp_curr - dprop_curr

    if (log(runif(1)) < log_prob) {
      acc <- acc + 1L
      W_curr <- W_prop
      S_tot_curr <- S_tot_prop
      lp_curr <- lp_prop
      mean_curr <- mean_prop
    }

    if (i > burnin && (i - burnin) %% thin == 0) {
      cnt <- (i - burnin) %/% thin
      sim[cnt, ] <- S_tot_curr
    }

    acc_prob[i] <- acc / i
    h <- max(1e-19, h + c1.h * i^(-c2.h) * (acc / i - 0.57))  # actually update h
    h.vec[i] <- h

    if (messages) {
      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, i)
      } else if (i %% max(1, floor(n_sim/20)) == 0) {
        message(sprintf("   %3d%%", round(100 * i / n_sim)))
      }
    }
  }

  if (!is.null(pb)) close(pb)
  if (messages) message(" done.\n")

  out_sim <- list(samples = list(S = sim[, 1:n_loc, drop = FALSE]))
  if (n_re > 0) {
    re_names <- if (!is.null(colnames(ID_re))) colnames(ID_re) else paste0("re", seq_len(n_re))
    for (i in 1:n_re) {
      out_sim$samples[[re_names[i]]] <- sim[, ind_re[[i]], drop = FALSE]
    }
  }

  out_sim$tuning_par <- h.vec
  out_sim$acceptance_prob <- acc_prob
  class(out_sim) <- "mcmc.RiskMap"

  return(out_sim)
}
