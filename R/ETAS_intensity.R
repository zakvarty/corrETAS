#' Calculate the conditional intenisty function of an ETAS point process
#'
#' @param t  Vector of times at which to calculate the intensity.
#' @param ts Vector of observed event times.
#' @param ms Vector of observed event magnitudes.
#' @param cGPDPars Vector of ETAS parameters (mu,C,a,nu_t,xi_t).
#' @param OmoriPars Alternative vector of ETAS parameters (mu,K,a,c,p).
#' @param M0 Minimum triggering magnitude.
#'
#' @return Vector of conditional intensities.
#' @export
#'
ETAS_intensity <- function(t, ts, ms,cGPDPars = NULL, OmoriPars = NULL, M0 = NULL){
  L <- outer(t,ts,FUN = '-')

  if(is.null(OmoriPars) & is.null(cGPDPars)){
    stop('Specify one of cGPDPars or OmoriPars.')
  } else if (!is.null(OmoriPars) & !is.null(cGPDPars)){
    stop('Specify only one of cGPDPars or OmoriPars.')
  } else if (is.null(OmoriPars) & !is.null(cGPDPars)){
   mu <- cGPDPars[1]
   C  <- cGPDPars[2]
   a  <- cGPDPars[3]
   nu_t  <- cGPDPars[4]
   xi_t  <- cGPDPars[5]

   mbar <- mean(ms)
   h <- dgpd(x = L, shape = xi_t, nu = nu_t, mu = 0)
   kappa <- C * exp(a*(ms-mbar))

  } else {
    if(is.null(M0)) stop('M0 must be specified to use OmoriPars.')
    mu <- OmoriPars[1]
    K  <- OmoriPars[2]
    a  <- OmoriPars[3]
    c_t  <- OmoriPars[4]
    p_t  <- OmoriPars[5]
    nuxi_t <- pow2nu(c_t,p_t)

    h <- dgpd(x = L, shape = nuxi_t[2], nu = nuxi_t[1], mu = 0)

    kappa <- K * exp(a*(ms-M0))
 }

   int1 <- rep(mu,length(t))
   int2 <- t(kappa) %*% t(h)
   int <-  int1 + as.vector(int2)

 return(data.frame(t = t, int = int))
}
