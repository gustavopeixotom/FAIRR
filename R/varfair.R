#' Estimate FAIR_G(K) Models
#'
#' @description
#' Functional Approximation of Impulse Responses (FAIR) using Gaussian basis functions by setting up a VAR(p) and extracting orthogonal shocks and an IRF to set priors. 
#'
#' @param Y A numeric vector or matrix. The dependent variable.
#' @param Z A numeric vector or matrix. The shock variable.
#' @param H Integer. Specifies the horizon of the IRF.
#' @param K Integer (1 or 2). Specifies the number of basis functions used.
#' @param asym Logical. If \code{TRUE}, fits separate IRFs for positive and negative shocks. If \code{FALSE}, fits a single symmetric IRF.
#' @param maxvarlag Integer. Maximum number of lags for the VAR used to identify orthogonal shocks.
#' @param varci Logical. Compute the VAR's IRF confidence intervals? Defaults to \code{FALSE}.
#' @param varboot Numeric. If \code{varci = TRUE}, specifying the runs for the bootstrap (passed to \code{\link[vars:irf]{vars::irf}}).
#' @param spanloess Numeric. Smoothing parameter for the VAR IRF prior (see \code{\link[stats:loess]{loess}}).
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
#'   \item{var}{IRF and shocks from the VAR estimated.}
#'   \item{stan}{Summary of the Stan model fit.}
#'   \item{data}{The data used for estimation.}
#' }
#' @export
#' @import rstan
#' @importFrom vars VAR VARselect irf
#' @importFrom pracma findpeaks
#' @importFrom parallel detectCores
#' @importFrom R.utils withTimeout
#' @importMethodsFrom rstan extract
varfair <- function(Y, Z, X = NULL, H, K = 1, asym = FALSE, maxvarlag = 12, varci = FALSE, varboot = 500, spanloess = 0.5, sdpar = 10, cores = NA, chains = 4, iter = 2000, adapt_delta = 0.95, max_treedepth = 12, timeout = 180) {
  start_time <- Sys.time() 
  
  message("Please cite the original FAIR paper and my github link !")
  
  if (is.na(cores) || cores < 1) {
    cores <- 1
  }
  if (cores > chains) {
    cores <- chains
    warning(paste("'cores' argument greater than 'chains' argument. rstan will use at most as many cores as chains; defaulting to", cores, "cores."))
  }
  
  if (cores == 1) {
    if (!asym) {
      message(paste0(
        "Estimating a symmetric FAIR(", K, 
        ") with stan (", iter, " iterations, ", chains, " Markov Chains, 1 core). This might take a while..."
      ))
    } else {
      message(paste0(
        "Estimating a pair of FAIR(", K, 
        ") for positive and negative shocks with stan (", iter, 
        " iterations, ", chains, " Markov Chains, 1 core). This might take a while..."
      ))
    }
  } else {
    if (!asym) {
      message(paste0(
        "Estimating a symmetric FAIR(", K, 
        ") with stan (", iter, " iterations,", chains, " Markov Chains, ", cores, 
        " cores). Watch your CPU temperature !"
      ))
    } else {
      message(paste0(
        "Estimating a pair of FAIR(", K, 
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
      
      if (!is.null(X)) {
        X <- as.matrix(X)
        dimnames(X) <- NULL
        if (is.null(colnames(X))) {
          colnames(X) <- paste0("X", seq_len(ncol(X)))
        }
        var_data <- as.data.frame(cbind(Z = as.numeric(Z), X, Y = as.numeric(Y)))
      } else {
        var_data <- as.data.frame(cbind(Z = as.numeric(Z), Y = as.numeric(Y)))
      }
      

  # 2 - VARs ----
      message("Getting the VAR representation and orthogonal shocks...")
      # Fitting
        order <- as.integer(vars::VARselect(var_data, lag.max=maxvarlag, type="trend")$selection["AIC(n)"])
        svar  <- vars::VAR(var_data, p=order, type="trend")
        
      # IRF priors from VAR's
        # Function that does it 
        irf.priors <- function(irf, K, spanloess) {
          irf_smooth <- predict(loess(irf ~ seq_along(irf), span = spanloess))
          if (K == 1) {
            b <- which.max(abs(irf)) 
            a <- irf_smooth[b] 
            half <- abs(a) / 2
            # Left half-width
            left_idx <- which(abs(irf[1:b]) <= half)
            if(length(left_idx) == 0){
              h_left <- 1
            } else {
              h_left <- max(left_idx)
            }
            # Right half-width
            right_idx <- which(abs(irf[b:length(irf)]) <= half)
            if(length(right_idx) == 0){
              h_right <- length(irf)
            } else {
              h_right <- b + min(right_idx) - 1
            }
            c <- mean(c(b - h_left, h_right - b)) / sqrt(2 * log(2))
            return(list(a = a, b = b, c = c))
          } else if (K == 2) {
            peaks_pos <- findpeaks(irf_smooth, nups = 1, ndowns = 0) 
            peaks_neg <- findpeaks(-irf_smooth, nups = 0, ndowns = 1)
            peak_vals <- c(); peak_h <- c()
            if (!is.null(peaks_pos)) { peak_vals <- c(peak_vals, peaks_pos[, 1]); peak_h <- c(peak_h, peaks_pos[, 2]) }
            if (!is.null(peaks_neg)) { peak_vals <- c(peak_vals, -peaks_neg[, 1]); peak_h <- c(peak_h, peaks_neg[, 2]) }
            if (length(peak_vals) < 2) stop("K=2 specified but fewer than 2 peaks found.")
            ord_mag <- order(abs(peak_vals), decreasing = TRUE)
            peak_vals <- peak_vals[ord_mag]; peak_h <- peak_h[ord_mag]
            if(length(peak_vals) > 2){ peak_vals <- peak_vals[1:2]; peak_h <- peak_h[1:2] }
            ord_time <- order(peak_h); peak_vals <- peak_vals[ord_time]; peak_h <- peak_h[ord_time]
            a1 <- peak_vals[1]; b1 <- peak_h[1]; a2 <- peak_vals[2]; b2 <- peak_h[2]
            get_c <- function(val_a, val_b, series) {
              half_val <- abs(val_a) / 2
              idx_l <- which(abs(series[1:val_b]) <= half_val); h_l <- if(length(idx_l) > 0) max(idx_l) else 1
              idx_r <- which(abs(series[val_b:length(series)]) <= half_val); h_r <- if(length(idx_r) > 0) (val_b + min(idx_r) - 1) else length(series)
              mean(c(val_b - h_l, h_r - val_b)) / sqrt(2 * log(2))
            }
            c1 <- get_c(a1, b1, irf_smooth); c2 <- get_c(a2, b2, irf_smooth)
            return(list(a1 = a1, b1 = b1, c1 = c1, a2 = a2, b2 = b2, c2 = c2))
          } else stop("K > 2 is not supported yet")
        }
        
        # IRF from full VAR
          if (varci == FALSE) {
            sirf_obj <- vars::irf(svar, impulse = "Z", response = c("Y"), n.ahead = H, ortho = TRUE, boot = FALSE) 
          } else {
            sirf_obj <- vars::irf(svar, impulse = "Z", response = c("Y"), n.ahead = H, ortho = TRUE, boot = TRUE, ci = 0.95, runs = varboot)
          }
          sirf <- sirf_obj$irf[["Z"]][, "Y"]
          priors <- irf.priors(sirf, K, spanloess)
      
        # Orthogonal shocks
          u <- residuals(svar)
          svar.vcov <- summary(svar)$covres
          svar.chol <- t(chol(svar.vcov))
          svar.chol.inv <- solve(svar.chol)
          svar.ortho <- u %*% t(svar.chol.inv)
          
          # Alignment: Partial out X (or Trend) from Y and Trim
          if (!is.null(X)) {
            y_clean_full <- residuals(lm(Y ~ X)) # if X is NULL
          } else {
            t <- seq_along(Y)
            y_clean_full <- residuals(lm(Y ~ t)) # if not
          }
          
          trim.y <- y_clean_full[(order + 1):length(y_clean_full)] # trim
          trim.df <- data.frame(Y = trim.y, shock = svar.ortho[, "Z"])
          
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
            if(K == 1) npriors$a <- -npriors$a        # flip 'a' for negative shocks
            if(K == 2) {                               # flip both peaks
              npriors$a1 <- -npriors$a1
              npriors$a2 <- -npriors$a2
            }
          }
  
  # 3 - Functional approximation ----
      # Getting stan data
        # Defining a stan data function for convenience
        stan_data <- function(trim.df, H, priors, K, sdpar) {
          out <- list(
            N = nrow(trim.df),
            Y = as.vector(trim.df$Y),
            Z = as.vector(trim.df$shock),
            H = H,
            sigma_scale = sd(trim.df$Y),
            sdpar = sdpar
          )
          if(K == 1) {
            out$prior_a <- priors$a
            out$prior_b <- priors$b
            out$prior_c <- priors$c
          } else {
            out$prior_a1 <- priors$a1
            out$prior_b1 <- priors$b1
            out$prior_c1 <- priors$c1
            out$prior_a2 <- priors$a2
            out$prior_b2 <- priors$b2
            out$prior_c2 <- priors$c2
          }
          return(out)
        }
      # Getting the data
        if(asym == FALSE) {
          sdata <- stan_data(trim.df, H, priors, K, sdpar)
        } else {
          pdata <- stan_data(trim.pdf, H, ppriors, K, sdpar)
          ndata <- stan_data(trim.ndf, H, npriors, K, sdpar)
        }
      
    # Running Stan
        # Forcing the initial stan guess to be Gaussians computed using priors
          get_inits <- function(chain_id) {
            if(K==1){
              list(a = priors$a, b = priors$b, c = priors$c, sigma = 1)
            } else {
              list(
                a1 = priors$a1, 
                b1 = priors$b1, 
                c1 = priors$c1, 
                
                a2 = priors$a2, 
                b2 = priors$b2, 
                c2 = priors$c2, 
                
                sigma = 1
              )
            }
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
          if (K == 1) {
            model <- stanmodels$FAIRone
          } else {
            model <- stanmodels$FAIRtwo
          }
        
        if (asym == FALSE) {
          message(paste0("Estimating symmetric FAIR(", K, ")"))
          fit <- run_fair(model_obj = model, sdata, chains, iter, cores, timeout, inits)
        } else {
          message(paste0("Estimating FAIR(", K, ") for positive shocks"))
          pfit <- run_fair(model_obj = model, pdata, chains, iter, cores, timeout, inits)
          
          message(paste0("Estimating FAIR(", K, ") for negative shocks"))
          nfit <- run_fair(model_obj = model, ndata, chains, iter, cores, timeout, inits)
        }
        
    
  # Outputs
    # Extract posterior draws from Stan
        extract_psi <- function(fit, K) {
          post <- rstan::extract(fit)
          # Parameters
          if (K == 1) {
            params <- cbind(
              a = post$a,
              b = post$b,
              c = post$c
            )
          } else if (K == 2) {
            params <- cbind(
              a1 = post$a1, b1 = post$b1, c1 = post$c1,
              a2 = post$a2, b2 = post$b2, c2 = post$c2
            )
          } else {
            stop("Unsupported K")
          }
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
      
    # Orthogonalized shocks 
      shock_all <- svar.ortho[, "Z"]
      
    # Asymmetric shocks (if requested)
      if (asym) {
        shock_pos <- ifelse(shock_all > 0, shock_all, 0)
        shock_neg <- ifelse(shock_all < 0, shock_all, 0)
      }
      
    # Assemble the final output list
      if (!asym) {
        output <- list(
          fair = extract_psi(fit, K = K),
          var = list(
            irf = list(
              mean = sirf,
              low  <- sirf_obj$Lower$Z[, "Y"],
              high <- sirf_obj$Upper$Z[, "Y"]
            ),
            shocks = list(all = shock_all),
            priors = priors
          ),
          stan = summary(fit),
          data = list(
            Y = Y,
            Z = Z,
            X = X
          )
        )
      } else {
        output <- list(
          fair = list(
            positive = extract_psi(pfit, K = K),
            negative = extract_psi(nfit, K = K)
          ),
          var = list(
            irf = list(
              mean = sirf,
              low  <- sirf_obj$Lower$Z[, "Y"],
              high <- sirf_obj$Upper$Z[, "Y"]
            ),
            shocks = list(
              all      = shock_all,
              positive = shock_pos,
              negative = shock_neg
            ),
            priors = priors
          ),
          negative.stan = summary(nfit),
          positive.stan = summary(pfit),
          data = list(
            Y = Y,
            Z = Z,
            X = X
          )
        )
      }
      
    end_time <- Sys.time()
    duration <- round(difftime(end_time, start_time, units = "secs"), 2)
    
    message(paste0("Estimation completed in ", duration, " seconds !"))
    
    return(output) 
}