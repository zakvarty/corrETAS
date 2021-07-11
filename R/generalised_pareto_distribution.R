#=====================================================================
# Functions for Generalised Pareto Distribution.
# Added option to use nu parameterisation
# checks that param values are valid
#=====================================================================
# pgpd
# qgpd
# dgpd
# rgpd
#=====================================================================
#' Generalised Pareto Distribution
#'
#' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or p<mu.
#'
#' @param q vector of quantiles.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter: nu = sigma/(1+xi)
#' @param mu  location parameter
#' @return Probability of the GPD X<=q
#' @author Zak Varty
#' @examples
#' pgpd(q = c(-1,1.5,3), shape = 1, scale = 1)
#' pgpd(q = 1.5, shape = c(0,-1), scale = c(0.1,1))
#' @export

pgpd<- function(q, shape, scale = NULL, nu = NULL, mu = 0){
  # one and only one of {nu, scale} may be sepecifed
  if( is.null(scale) & is.null(nu)){
    stop('Define one of the parameters nu or scale.')
  }
  if( !is.null(scale) & !is.null(nu)){
    stop('Define only one of the parameters nu and scale.')
  }
  # Calculate scale from nu if required
  if( !is.null(nu) & is.null(scale)){
    scale <- nu/(1+shape)
    if(any(scale <= 0)){
      stop('Implied scale parameter(s) must be positive.')
    }

  }
  # Check that scale value(s) are positive
  if(any(scale<=0)){
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure q, scale, shape and mu are of same length.
  if(length(scale) == 1 & length(q) > 1){
    scale <- rep(scale, length(q))
  }
  if(length(shape) == 1 & length(q) > 1){
    shape <- rep(shape, length(q))
  }
  if(length(mu) == 1 & length(q) > 1){
    mu <- rep(mu, length(q))
  }

  #calculate probabilities
    p <- (1 - (1 + (shape * (q - mu))/scale)^(-1/shape))
  #correct probabilities below mu or above upper end point
    p[q<mu] <- 0
    p[(shape < 0) & (q >= (mu - scale/shape))] <- 1

  #correct probabilities where xi = 0
    ex <- which(shape ==0)
    p[ex] <- stats::pexp(q = q[ex]- mu[ex], rate = 1/scale[ex])

  return(p)
}

#' Generalised Pareto Distribution
#'
#' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or p is not a valid
#' probability.
#'
#' @param p vector of quantiles.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter: nu = sigma/(1+xi)
#' @param mu  location parameter
#' @return Probability of the GPD X<=x
#' @author Zak Varty
#' @examples
#' qgpd(p = 0.5, shape = 0.5, scale = 0.5)
#' \dontrun{ qgpd(p = -0.1, shape = 0, scale = 1, mu = 0.1) }
#' @export
qgpd<- function(p, shape, scale = NULL, nu = NULL, mu = 0){
  # one and only one of {nu, scale} may be sepecifed
  if( is.null(scale) & is.null(nu)){
    stop('Define one of the parameters nu or scale.')
  }
  if( !is.null(scale) & !is.null(nu)){
    stop('Define only one of the parameters nu and scale.')
  }

  # Probabilities must all be positive
  if(!all((p>=0)&(p<=1))){
    stop('Probabilities p must be in the range [0,1].')
  }

  # Calculate scale from nu if required
  if( !is.null(nu) & is.null(scale)){
    scale <- nu/(1+shape)
    if(any(scale <= 0)){
      stop('Implied scale parameter(s) must be positive.')
    }

  }

  # Check that scale value(s) are positive
  if(any(scale<=0)){
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure p, scale, shape and mu are of same length.
  if(length(scale) == 1 & length(p) > 1){
    scale <- rep(scale, length(p))
  }
  if(length(shape) == 1 & length(p) > 1){
    shape <- rep(shape, length(p))
  }
  if(length(mu) == 1 & length(p) > 1){
    mu <- rep(mu, length(p))
  }

  #calculate quantiles
  q <- mu + (scale/shape) * ((1 - p)^(-shape) - 1)

  #correct quantiles where xi = 0
  ex <- which(shape ==0)
  q[ex] <- mu[ex] + stats::qexp(p = p[ex],rate = 1/scale[ex])

  return(q)
}



#' Generalised Pareto Distribution
#'
#' Density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or x is
#' outside of the domain of the given distribution.
#'
#' @param x vector of values as which to evaluate density.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter
#' @param mu  location parameter
#' @param log  locical. Return log
#' @return density of the GPD at x
#' @author Zak Varty
#' @examples
#' dgpd(x = c(-1,0.5,1,1.9,5),shape = -0.5, scale = 1)
#' @export
#'

dgpd<- function(x, shape, scale = NULL, nu = NULL, mu = 0, log = FALSE){
  # one and only one of {nu, scale} may be sepecifed
  if( is.null(scale) & is.null(nu)){
    stop('Define one of the parameters nu or scale.')
  }
  if( !is.null(scale) & !is.null(nu)){
    stop('Define only one of the parameters nu and scale.')
  }

  # Calculate scale from nu if required
  if( !is.null(nu) & is.null(scale)){
    scale <- nu/(1+shape)
    if(any(scale <= 0)){
      stop('Implied scale parameter(s) must be positive.')
    }
  }

  # Check that scale value(s) are positive
  if(any(scale<=0)){
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure x, scale, shape and mu are of same length.
  if(length(scale) == 1 & length(x) > 1){
    scale <- rep(scale, length(x))
  }
  if(length(shape) == 1 & length(x) > 1){
    shape <- rep(shape, length(x))
  }
  if(length(mu) == 1 & length(x) > 1){
    mu <- rep(mu, length(x))
  }

  if(log == FALSE){
    out <- (scale^(-1)) * pmax((1 + shape * (x - mu)/scale),0)^((-1/shape) - 1)
    #ammend values below threshold
    out[which(x < mu)] <- 0
    #ammend values above upper endpoint (if it exists)
    out[which((shape < 0) & (x >= (mu - scale/shape)))] <- 0
    #ammend values where xi = 0 (if they exist)
    ex <- which(shape ==0)
    out[ex] <- stats::dexp(x = x[ex]-mu[ex], rate = 1/scale[ex])
  } else {
    out <-  -log(scale) + ((-1/shape) - 1)*log(pmax((1 + shape * (x - mu)/scale),0))
    #ammend values below threshold
    out[which(x < mu)] <- -Inf
    #ammend values above upper endpoint (if it exists)
    out[which((shape < 0) & (x >= (mu - scale/shape)))] <- -Inf
    #ammend values where xi = 0 (if they exist)
    ex <- which(shape ==0)
    out[ex] <- stats::dexp(x = x[ex]-mu[ex], rate = 1/scale[ex],log = TRUE)
  }
  return(out)
}

#' Generalised Pareto Distribution
#'
#' Sample the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0.
#'
#' @param n sample size.
#' @param shape shape parameter (xi).
#' @param scale scale parameter (sigma).
#' @param nu  alternative scale parameter.
#' @param mu  location parameter.
#' @return Random sample from generalised pareto distirbution.
#' @author Zak Varty
#' @examples
#' rgpd(n = 100, shape = 0, scale = 1:100)
#' @export
rgpd<- function(n, shape, scale = NULL, nu = NULL, mu = 0){
  ## Input checks
  # one and only one of {nu, scale} may be sepecifed
  if( is.null(scale) & is.null(nu)){
    stop('Define one of the parameters nu or scale.')
  }
  if( !is.null(scale) & !is.null(nu)){
    stop('Define only one of the parameters nu and scale.')
  }
  # Calculate scale from nu if required
  if( !is.null(nu) & is.null(scale)){
    scale <- nu/(1+shape)
    if(any(scale <= 0)){
      stop('Implied scale parameter(s) must be positive.')
    }
  }
  # Check that scale value(s) are positive
  if(any(scale<=0)){
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure q, scale, shape and mu are of same length.
  if((length(scale) == 1) & (n > 1)){
    scale <- rep(scale, n)
  }
  if((length(shape) == 1) & (n > 1)){
    shape <- rep(shape, n)
  }
  if((length(mu) == 1) & (n > 1)){
    mu <- rep(mu, n)
  }

  #simulate sample
    sample <- mu + (scale/shape) * ((1 - stats::runif(n))^(-shape) - 1)
  #correct sample values where xi = 0
    ex <- which(shape ==0)
    sample[ex] <- mu[ex] +
      stats::rexp(n = length(ex),rate = 1/scale[ex])

  return(sample)
}


# -----
#' Calculate the log-likelihood of iid generalised Pareto data
#'
#' @param x Data vector
#' @param nuxi Parameters of GPD
#' @param mu Threshold of GPD
#' @param neg Logical. Return negative log-likelihood?
#'
#' @return (negative) log-likelihood at (nu,xi)
llh_gpd <- function(x, nuxi, mu, neg = FALSE){
  nu <- nuxi[1]
  xi <- nuxi[2]

  sig <- nu2sig(nu,xi)[1]
  if(sig<=0){return(1e10)}
  if(xi<0 & any(x > mu - sig/xi)){return(1e10)}

  llh <- dgpd(
    x = x,
    nu = nu,
    shape = xi,
    mu = mu,
    log = TRUE
  )
  llh <- sum(llh)

  return( (-1)^neg * llh)
}


#' Find the maximum likelihood estimate for iid generalised Pareto data
#'
#' @param nuxi Starting point for optimisation of the likelihood function
#' @param x Observations from the generalised pareto distribution
#' @param mu threshold for the generealised pareto distribution
#' @param hessian Logical. Return hessian in output? Default is false.
#' @param ... Additional arguments to optim.
#'
#' @return list of maximum likelihood estimate and log-likelihood value at the mle.
#' @export
#'
#' @examples
#' dat <- rgpd(n = 300, nu = 1.2, shape = 0, mu = 0)
#' maxLikelihoodGPD(x = dat, mu = 0, nuxi = c(0.1,0))
maxLikelihoodGPD <- function(nuxi, x, mu, hessian = FALSE, ...){

  nu <- nuxi[1]
  xi <- nuxi[2]
  stopifnot(all(x>=mu))
  sig <- nu2sig(nuxi,xi)[1]
  stopifnot(sig>0)
  if(xi<0 & any(x > mu - sig/xi)){stop("Invalid starting point. Data above upper end point of distribution.")}


  temp <- optim(fn = llh_gpd,
                par = nuxi,
                x = x,
                mu = mu,
                neg = TRUE,
                hessian = hessian,
                ...)
  out <- list(params = temp$par, loglik = -temp$value)

  if(hessian) out$hessian <- temp$hessian

  return(out)
}

