#' Convert Paramters from powerlaw (c,p) to gpd (nu,xi)
#'
#' @param c power law length scale
#' @param p power law decay
#'
#' @return vector of converted parameters
#' @author Zak Varty
#' @export
#'
pow2nu <- function(c,p){
  nu <- (c/(p-1))*(1 + 1/(p-1))
  xi <- 1/(p-1)
  out <- c(nu,xi)
  return(out)
}
#' Convert parameters from gpd (nu, xi) to pwoer law (c,p)
#'
#' @param nu gpd proxy scale
#' @param xi gpd shape
#'
#' @return vector of converted parameters
#' @author Zak Varty
#' @export
#'
nu2pow <- function(nu, xi){
  c <- nu / (xi * (1+ xi))
  p <- 1/xi + 1
  out <- c(c,p)
  return(out)
}
#' Convert parameters from gpd (nu,xi) to gpd (sigma,xi)
#'
#' @param nu gpd proxy scale
#' @param xi gpd shape
#'
#' @return vector of converted parameters
#' @author Zak Varty
#' @export
#'
nu2sig <- function(nu, xi){
  sig <- nu / (1+ xi)
  xi <- xi
  out <- c(sig, xi)
  return(out)
}
#' Convert prameters from gpd (sigma,xi) to gpd (nu,xi)
#'
#' @param sig gpd scale
#' @param xi gpd shape
#'
#' @return vector of converted parameters
#' @author Zak Varty
#' @export
#'
sig2nu <- function(sig, xi){
  nu <- sig * (1 + xi)
  xi <- xi
  out <- c(nu, xi)
  return(out)
}
#' Convert prameters from powerlaw (c,p) to gpd (sigma,xi)
#'
#' @param c power law length scale
#' @param p power law decay
#'
#' @return vector of converted parameters
#' @author Zak Varty
#' @export
#'
pow2sig <- function(c, p){
  sig <- c/(p-1)
  xi <- 1/(p-1)
  out <- c(sig, xi)
  return(out)
}

#' Convert prameters from gpd (sigma,xi) to powerlaw (c,p)
#'
#' @param sig gpd scale
#' @param xi  gpd shape
#'
#' @return vector of converted parameters
#' @author Zak Varty
#' @export
#'
sig2pow <- function(sig, xi){
  c <- sig/xi
  p <- 1/xi + 1
  out <- c(c,p)
  return(out)
}

#-------------------------------------------------------------------------------

#' Convert ETAS parameter vector from Omori-law to cGPD parameterisation.
#'
#' @param OmoriPars Vector of Omori-law ETAS parameters (mu, K, alpha, c, p)
#' @param Mbar Mean magnitude of catalogue
#' @param M0 Threshold magnitude of catalogue
#'
#' @return Vector of cGPD parameters (mu, C, alpha, nu_t, xi_t)
#' @export
Omori2cGPD <- function(OmoriPars,Mbar,M0){
  stopifnot(length(OmoriPars) == 5)

  cGPDPars <- OmoriPars
  #K -->> C
  cGPDPars[2] <- OmoriPars[2]*exp(OmoriPars[3]*(Mbar-M0))
  #c,p -->> nu,xi
  cGPDPars[4:5] <- pow2nu(OmoriPars[4],OmoriPars[5])

  return(cGPDPars)
}


#' Convert ETAS parameter vector from cGPD to Omori-law parameterisation.
#'
#' @param cGPDPars Vector of Omori-law ETAS parameters (mu, C, alpha, nu_t, xi_t)
#' @param Mbar Mean magnitude of catalogue
#' @param M0 Threshold magnitude of catalogue
#'
#' @return Vector of cGPD parameters (mu, K, alpha, c, p)
#' @export
#'
cGPD2Omori <- function(cGPDPars,Mbar,M0){
  stopifnot(length(cGPDPars) == 5)

  OmoriPars <- cGPDPars
  #C -->> K
  OmoriPars[2] <- cGPDPars[2]*exp(OmoriPars[3]*(M0-Mbar))
  #nu,xi -->> c,p
  OmoriPars[4:5] <- nu2pow(cGPDPars[4],cGPDPars[5])

  return(cGPDPars)
}
