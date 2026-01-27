#' Estimate FAIR_E(1) Models with Exogenous Shocks
#'
#' @description
#' Functional Approximation of Impulse Responses (FAIR) with exponential basis functions using shocks that are already considered exogenous.
#'
#' @param Y A numeric vector or matrix. The dependent variable.
#' @param Z A numeric vector or matrix. The exogenous shock variable.
#' @param H Integer. Specifies the horizon of the IRF.
#' @param priorvec A numeric vector of length 2 containing priors : \code{c(a, l)}
#' @param detrend Integer (0, 1, or 2). Controls data transformation:
#'   \itemize{
#'     \item \code{0}: No trend removed.
#'     \item \code{1}: Log-difference (\code{diff(log(Y))}).
#'     \item \code{2}: Linear time trend removed (residuals of \code{lm(Y ~ time)}).
#'   }
#' @param asym Logical. If \code{TRUE}, fits separate IRFs for positive and negative shocks. If \code{FALSE}, fits a single symmetric IRF.
#' @param sdpar Numeric. Standard deviation of the Normal prior for parameters b and c.
#' @param cores Integer. Number of CPU cores for Stan. If greater than \code{chains}, defaults to \code{chains}.
#' @param chains Integer. Number of Markov chains (see \code{\link[rstan:stan]{rstan::stan}}).
#' @param iter Integer. Number of iterations per chain (see \code{\link[rstan:stan]{rstan::stan}}).
#' @param adapt_delta Numeric. Target average acceptance probability (see \code{\link[rstan:stan]{rstan::stan}}).
#' @param max_treedepth Integer. Maximum tree depth (see \code{\link[rstan:stan]{rstan::stan}}).
#' @param timeout Numeric. Maximum seconds allowed for the Stan run.
#'
#' @return A list containing:
#' \describe{
#'   \item{fair}{Posterior draws of the FAIR coefficients and corresponding IRF.}
#'   \item{stan}{Summary of the Stan model fit.}
#'   \item{data}{The processed data used for estimation.}
#' }
#' @export
#' @import rstan
#' @importFrom parallel detectCores
#' @importFrom R.utils withTimeout
#' @importMethodsFrom rstan extract
fairE <- function(Y, Z, H, priorvec = c(1, 2), detrend = 0, asym = FALSE, sdpar = 20, cores = NA, chains = 4, iter = 2000, adapt_delta = 0.95, max_treedepth = 12, timeout = 180) {
  start_time <- Sys.time()

  message("Please cite the original FAIR paper and my github !")

  # cores check
  if (is.na(cores) || cores < 1) {
    cores <- 1
  } # no specified

  if (cores > chains) {
    cores <- chains
    warning(paste("'cores' argument greater than 'chains' argument. rstan will use at most as many cores as chains; defaulting to", cores, "cores."))
  } # unnecessary cores

  if (detrend > 2) {
    stop("detrend cannot be more than 2")
  }

  # Stop if priors arent a vector of length 2K
  if (nrow(as.data.frame(priorvec)) != 2 || ncol(as.data.frame(priorvec)) != 1) {
    stop("The vector of priors must be of length 2 and width 1")
  }
  priorvec <- as.vector(priorvec)

  # Yapping
  if (cores == 1) {
    if (!asym) {
      message(paste0(
        "Estimating a symmetric FAIR_E(", 1,
        ") with stan (", iter, " iterations, ", chains, " Markov Chains, 1 core). This might take a while..."
      ))
    } else {
      message(paste0(
        "Estimating a pair of FAIR(", 1,
        ") for positive and negative shocks with stan (", iter,
        " iterations, ", chains, " Markov Chains, 1 core). This might take a while..."
      ))
    }
  } else {
    if (!asym) {
      message(paste0(
        "Estimating a symmetric FAIR(", 1,
        ") with stan (", iter, " iterations,", chains, " Markov Chains, ", cores,
        " cores). Watch your CPU temperature !"
      ))
    } else {
      message(paste0(
        "Estimating a pair of FAIR(", 1,
        ") for positive and negative shocks with stan (", iter,
        " iterations, ", chains, " Markov Chains, ", cores, " cores). Watch your CPU temperature !"
      ))
    }
  }

  message(paste("Estimation will stop if stan takes longer than", timeout, "seconds. You can adjust this setting with the argument 'timeout'."))

  if (cores > chains) {
    cores <- chains
    warning(paste("Stan will use at most as many cores as there are Markov Chains. This estimation will use only", cores, "cores."))
  }

  # 1 - Data ----
  Y <- as.matrix(Y) # Dependent variable
  Z <- as.matrix(Z) # Shock

  # Check widths
  if (ncol(Y) > 1) {
    stop("Y must be a univariate time series (width = 1).")
  }
  if (ncol(Z) > 1) {
    stop("Z must be a univariate time series (width = 1).")
  }

  # Check length consistency
  if (nrow(Y) != nrow(Z)) {
    stop("Y and Z must have the same number of observations.")
  }

  # 2 - Priors ---
  priors <- list()
    priors$a <- priorvec[1]
    priors$l <- priorvec[2]

  # Partialing out the trends
  if (detrend == 0) { # No trend
    trim.df <- data.frame(Y = Y, shock = Z)
  } else {
    if (detrend == 1) { # Log Difference
      trim.y <- diff(log(Y))
      trim.df <- data.frame(Y = trim.y, shock = Z[-1])
    } else { # Linear in time trend
      t <- seq_along(Y)
      trim.y <- residuals(lm(Y ~ t))
      trim.df <- data.frame(Y = trim.y, shock = Z)
    }
  }

  # Asymmetric shocks for FAIR
  if(asym == TRUE) {
    # Positive shocks
    trim.pdf <- trim.df
    trim.pdf$shock <- ifelse(trim.df$shock > 0, trim.df$shock, 0)
    ppriors <- priors

    # Negative shocks
    trim.ndf <- trim.df
    trim.ndf$shock <- ifelse(trim.df$shock < 0, trim.df$shock, 0)
    npriors <- priors
  }


  # 3 - Functional approximation ----
  # Getting stan data
  # Defining a stan data function for convenience
  stan_data <- function(trim.df, H, priors, sdpar) {
    out <- list(
      N = nrow(trim.df),
      Y = as.vector(trim.df$Y),
      Z = as.vector(trim.df$shock),
      H = H,
      sigma_scale = sd(trim.df$Y),
      sdpar = sdpar,
      prior_a = priors$a,
      prior_l = priors$l
    )
    return(out)
  }
  # Getting the data
  if(asym == FALSE) {
    sdata <- stan_data(trim.df, H, priors, sdpar)
  } else {
    pdata <- stan_data(trim.pdf, H, ppriors, sdpar)
    ndata <- stan_data(trim.ndf, H, npriors, sdpar)
  }

  # Running Stan
  # Forcing the initial stan guess to match the Exponential priors
  get_inits <- function(chain_id) {
    list(
      a = priors$a,
      l = priors$l,
      sigma = 1
    )
  }
  inits <- lapply(1:chains, get_inits)

  # Function
  run_fair <- function(model_obj, data_list, chains, iter, cores, timeout, inits) {
    R.utils::withTimeout(
      rstan::sampling(
        object = model_obj,   # Use the compiled object, not a file path
        data = data_list,
        chains = chains,
        iter = iter,
        init = inits,
        refresh = 0,
        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
        cores = cores
      ),
      timeout = timeout,
      onTimeout = "error"
    )
  }

  # Run FAIR
    model <- stanmodels$FAIRexponential

  if (asym == FALSE) {
    message(paste0("Estimating symmetric FAIR_E(", 1, ")"))
    fit <- run_fair(model_obj = model, sdata, chains, iter, cores, timeout, inits)
  } else {
    message(paste0("Estimating FAIR_E(", 1, ") for positive shocks"))
    pfit <- run_fair(model_obj = model, pdata, chains, iter, cores, timeout, inits)

    message(paste0("Estimating FAIR_E(", 1, ") for negative shocks"))
    nfit <- run_fair(model_obj = model, ndata, chains, iter, cores, timeout, inits)
  }


  # Outputs
  # Extract posterior draws from Stan
  extract_psi <- function(fit) {
    post <- rstan::extract(fit)
    # Parameters
      params <- cbind(
        a = post$a,
        l = post$l
      )
    # Fit
    psi <- post$psi

    list(
      params = list(
        mean = colMeans(params),
        low  = apply(params, 2, quantile, 0.025),
        high = apply(params, 2, quantile, 0.975)
      ),
      psi = list(
        mean = colMeans(psi),
        low  = apply(psi, 2, quantile, 0.025),
        high = apply(psi, 2, quantile, 0.975)
      )
    )
  }

  # Assemble the final output list
  if (!asym) {
    output <- list(
      fair = extract_psi(fit),
      stan = fit,
      data = list(
        Y = trim.df$Y,
        Z = trim.df$shock
      )
    )
  } else {
    output <- list(
      fair = list(
        positive = extract_psi(pfit),
        negative = extract_psi(nfit)
      ),
      negative.stan = nfit,
      positive.stan = pfit,
      data = list(
        Y = trim.df$Y,
        neg.Z = trim.ndf$shock,
        pos.Z = trim.pdf$shock
      )
    )
  }

  end_time <- Sys.time()
  duration <- round(difftime(end_time, start_time, units = "secs"), 2)

  message(paste0("Estimation completed in ", duration, " seconds !"))

  return(output)
}
