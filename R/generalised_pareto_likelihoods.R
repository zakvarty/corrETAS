#' Calculate generalised Pareto likelihoods, log-likelihoods and negative log-likelihoods.
#'
#' @param x vector of observed values
#' @param sig vector of GPD scale parameters. Specify exactly one of sig and nu.
#' @param nu vector of alternative GPD scale parameters. Specify exactly one of sig and nu.
#' @param xi vector of GPD shape parameters.
#' @param mu vector of GPD thresholds.
#' @param logged Logical, stating wheter to return logged value. Default is FALSE.
#' @param negative Logical, stating wheter to return negated value. Default is FALSE.
#' @param combined Logical, stating wheter to return combined value. Default is TRUE. Liklihood values are combined by multiplication, log-likelihood values are combined by addition.
#'
#' @return vector of Liklihood, log-likelihood or negative log-likelihood value(s).
#' @export
#'
gpd_likelihood <- function(x,sig = NULL, nu = NULL, xi, mu, logged = FALSE, negative = FALSE, combined = TRUE){

  UEP_failure <- (xi < 0) &  (max(x) > mu - nu / (xi * (1 + xi)))
  if(is.null(nu)){  sigma_failure <-  sig < 0}
  if(is.null(sig)){ sigma_failure <-  (nu / (1 + xi)) < 0}

  if(UEP_failure | sigma_failure){
    L <- -1e10
  } else {
    L <- dgpd(x = x, scale = sig, shape = xi, nu = nu, mu = mu, log = logged)
  }
  if(combined & logged ) L <- sum(L)
  if(combined & !logged) L <- prod(L)
  if(negative) L <- -L
  return(L)
}

#' Evaluate flat prior density for gpd parameters
#'
#' @param sig gpd scale parameter
#' @param nu gpd alt-scale parameter
#' @param xi gpd shape parameter
#' @param mu gpd threshold
#' @param logged logical. Return log-likelihood?
#' @param optim_friendly logical. avoid -Inf in prior?
#'
#' @return (log-)prior density
#' @export
#'
#' @examples
#' gpd_flat_prior(nu = 0.1,xi = 0)
gpd_flat_prior <- function(sig = NULL, nu = NULL, xi, mu = 0, logged = FALSE, optim_friendly = FALSE){

  # error if sigma and nu both null
  stopifnot(exactly_one_null(sig,nu))

  # if only nu calculate sigmas
  if(is.null(sig)) sig <- nu / (1+ xi)

  out_of_bounds_value <- 0 + 1e-8 * optim_friendly

  # evaluate prior
  prior_val  <- rep(NA, length(sig))
  prior_val[which(sig <= 0)] <-  out_of_bounds_value
  prior_val[which(sig > 0)]  <- 1

  # log if requested
  if(logged) prior_val <- log(prior_val)

  return(prior_val)
}





#' Check that exactly one of two objects is NULL.
#'
#' @param a first object
#' @param b second object
#'
#' @return logical. Is exactly one a and b NULL?
#' @export
#'
#' @examples exactly_one_null(3, NULL)
exactly_one_null <- function(a,b){
  (!is.null(a) & is.null(b)) | (is.null(a) & !is.null(b))
}


