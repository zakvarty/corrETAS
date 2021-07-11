#' Testing C++ GPD cdf calculation
#'
#' @param q Numeric vector of values at which to print CDF
#' @param nu gpd alt-scale parameter
#' @param xi gpd shape parameter
#' @param mu gpd threshold
#'
#' @return NULL prints CDF value(s) at q.
#' @export
#'
cpp_pgpd <- function (q, nu, xi ,mu){
  n <- length(q)
  out <- rep(NA, n)

  for (i in 1:n){
    res <- .C('pgpdC',
              q = as.double(q[i]),
              nu = as.double(nu[i]),
              xi = as.double(xi[i]),
              mu = as.double(mu[i]),
              PACKAGE = 'corrETAS')
  }
}

#-------------------------------------------------------------------------------
#'  Conversion from GPD to Gaussian margin in C++
#'
#' @param x_gpd Numeric vector of gpd observations
#' @param nu gpd alt-scale parameter
#' @param xi gpd shape parameter
#' @param mu gpd threshold
#'
#' @return NULL prints result of Cpp function to put x_gpd on gaussian margins
#' @export
#'
cpp_gpd_to_gauss <- function (x_gpd, nu, xi ,mu){
  n <- length(x_gpd)
  out <- rep(NA_real_, n)

  for (i in 1:n){
    res <- .C('gpd_to_gaussC',
                 x_gpd = as.double(x_gpd[i]),
                 nu = as.double(nu[i]),
                 xi = as.double(xi[i]),
                 mu = as.double(mu[i]),
                 PACKAGE = 'corrETAS')
  }
  return(out)
}

#-------------------------------------------------------------------------------
#' Conversion from GPD to Gaussian margin in R
#'
#' @param x Numeric vector of gpd observations
#' @param sig gpd alt-scale parameter
#' @param xi gpd shape parameter
#' @param mu gpd threshold
#'
#' @return NULL prints result of Cpp function to put x_gpd on gaussian margins
#' @export
gpd_to_gaus <- function(x, sig, xi, mu){
  xu <- pgpd(q = x - mu, scale = sig, shape = xi)
  tail <- which(xu > 0.999)
  xg <- stats::qnorm(p = xu)
  xg[tail] <- gpd_to_gaus_tail(x = x[tail],sig = sig, xi = xi, mu = mu)
  return(xg)
}

#-------------------------------------------------------------------------------
#' Conversion from GPD to Gaussian margin for tail values in R (internal function)
#'
#' @param x Numeric vector of gpd observations
#' @param sig alt-scale parameter
#' @param xi gpd shape parameter
#' @param mu gpd threshold
#'
#' @return x on gaussian margins
#'
gpd_to_gaus_tail<- function(x, sig, xi, mu){
  xp  <- x - mu

  if(xi != 0){
    xn0 <- sqrt((2/xi)*log(1 + (xi/sig)*xp))
    epsilon <- (-log((4*pi/xi)*log(1+(xi/sig)*xp))) / ((4/xi)*log(1 + (xi/sig)*xp))
  } else {
    xn0 <- sqrt(2*xp/sig)
    epsilon <- log(2*pi*xp/sig) / (4*xp/sig)
  }

  xn <- xn0 + epsilon
  return(xn)
}
#-------------------------------------------------------------------------------
