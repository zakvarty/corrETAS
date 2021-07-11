#' Independent Metropolis sampler for nu-parameterised generalised Pareto distribution.
#'
#' @param x vector of observed data
#' @param threshold scalar threshold value
#' @param init initial values of parameters (nu, xi)
#' @param step_sds standard deviations of normal proposal density for each of (nu, xi)
#' @param n_samples integer number of samples to draw from the Markov chain
#' @param prior function like corrETAS::gpd_flat_prior, giving the joint proir on (nu,xi)
#' @param verbose integer/logical indicating how frequently progress and acceptance probability should be printed to the console.
#' @author Zak Varty
#' @return a data.frame of sampled parameter values
#' @export
#'
gpd_mcmc_nu <- function(x, threshold, init, step_sds, n_samples, prior = corrETAS::gpd_flat_prior, verbose = FALSE){
  # make vector of threshold exceedances
  y <- x[x > threshold] - threshold
  par_curr <- init
  n_accepted <- 0
  # Check good starting values:
    # at least one threshold exceedance
    stopifnot(length(y) >= 1)
    # Impied sigma is positive
    implied_sigma_is_positive <- init[1] / (1 + init[2]) > 0
    stopifnot(implied_sigma_is_positive)
    # all values less than upper endpoint
    if(init[2] < 0){
      UEP <- -init[1] / (init[2] * (1 + init[2]))
      stopifnot(max(y) < UEP)
    }

  # make storage for output
  samples <- data.frame(
    nu = rep(NA_real_, n_samples),
    xi = rep(NA_real_, n_samples),
    lpost = rep(NA_real_, n_samples))

  # calculate initial log-posterior
  lprior_curr <- prior(nu = par_curr[1], xi = par_curr[2], mu = 0, logged = TRUE)
  llik_curr <- corrETAS::gpd_likelihood(x = y, nu = par_curr[1], xi = par_curr[2], mu = 0, logged = TRUE)
  lpost_curr <- lprior_curr + llik_curr

  # FOR each iteration
  for (i in 1:n_samples){
    # propose new (nu,xi) pair
    par_prop <- par_curr + rnorm(n = 2, mean = c(0,0), sd = step_sds)
    # Calculate new log-posterior
    lprior_prop <- prior(nu = par_prop[1], xi = par_prop[2], mu = 0, logged = TRUE)
    llik_prop <- corrETAS::gpd_likelihood(x = y, nu = par_prop[1], xi = par_prop[2], mu = 0, logged = TRUE)
    lpost_prop <- lprior_prop + llik_prop
    # IF if (runif(0,1) < exp(new logposterior - curr log-posterior))
    if (runif(1) < exp(lpost_prop - lpost_curr)){
      lpost_curr <- lpost_prop
      par_curr <- par_prop
      n_accepted <- n_accepted + 1
    }
    # Record current parameter values
    samples[i,] <- c(par_curr, lpost_curr)
    # IF verbose
    if( (verbose > 0) & (i %% verbose == 0)){
      cat("Acceptance:", round(n_accepted/i, 4), "\t Iteration: ", i, ".\n")
    }
  }
  # Return sampled parameter values
  return(samples)
}



# ----- TESTING
# debug(gpd_mcmc_nu)
# mag_mcmc <- gpd_mcmc_nu(
#   x = cat$ms[cat$b == 0],
#   threshold = 1.5,
#   init = c(1,0),
#   step_sds = c(0.05,0.05),
#   n_samples = 10000,
#   verbose = TRUE)
#
# plot(mag_mcmc$nu, type = 'l')
# abline(h = mag_values[1], col = 2, lwd = 2)
# plot(mag_mcmc$xi, type = "l")
# abline(h = mag_values[2], col = 2, lwd = 2)
# head(mag_mcmc$lpost, 20)
