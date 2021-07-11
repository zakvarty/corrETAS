#' Estimate the parameters of the centred-GPD ETAS model using maximum likelihood.
#'

#'The Epidemic Type Aftershock Sequence (ETAS) model is widely used to quantify the degree of seismic activity in a geographical region, and to forecast the occurrence of future mainshocks and aftershocks (Ross 2016). The temporal ETAS model is a point process where the probability of an earthquake occurring at time \eqn{t} depends on the previous seismicity \eqn{H_t}{Ht}, and is defined by the conditional intensity function:

#' \deqn{ \lambda(t|H_t) = \mu + \sum_{t[i] < t}  \kappa(m[i]|C,\alpha) h(t[i]|\nu,\xi)}{ \lambda(t|Ht) = \mu + \sum \kappa(m[i]|C,\alpha) h(t[i]|\nu,\xi).}

#' This function uses an orthogonal parameterisation where:

#' \deqn{\kappa(m_i|C,\alpha) = Ce^{\alpha \left( m_i-\bar m \right)}}{\kappa(m[i]|C,\alpha) = C * exp(\alpha(m[i]- mean(m)))}

#' and \eqn{ h(t_i|\nu,\xi)} is the generalised Pareto density function reparameterised with \eqn{\nu = \sigma(1+\xi)}.

#' \deqn{}

#' The summation is over all previous earthquakes that occurred in the region, with the i'th such earthquake occurring at time \eqn{ t_i}{t[i]} and having magnitude \eqn{ m_i}{m[i]}. The quantity \eqn{M_0}{M0} denotes the magnitude of completeness of the catalog, so that \eqn{m_i \geq  M_0}{m[i] \ge M0} for all i. The temporal ETAS model has 5 parameters: \eqn{\mu} controls the background rate of seismicity, \eqn{ C} and \eqn{ \alpha} determine the productivity (average number of aftershocks) of an earthquake with magnitude \eqn{m}, and \eqn{\nu} and \eqn{ \xi}  are the parameters of the generalised Pareto distribution and represent the speed at which the aftershock rate decays over time. Each earthquake is assumed to have a magnitude which is an independent draw from the Gutenberg-Richter law \eqn{ p(m_i) = \beta e^{\beta(m_i-M_0)}}{ p(m) = \beta * exp(\beta(m-M0)}.

#' \deqn{}

#' This function estimates the parameters of the cGPD ETAS model using maximum likelihood

#' @param ts Vector containing the earthquake times
#' @param magnitudes Vector containing the earthquake magnitudes
#' @param M0 Magnitude of completeness.
#' @param maxTime Length of the time window [0,T] the catalog was observed over. If not specified, will be taken as the time of the last earthquake.
#' @param initval Initial value at which to start the estimation. A vector, with elements (mu, C, alpha, nu, xi)
#' @param displayOutput If TRUE then prints the out the likelihood during model fitting.
#' @param constrainOmori If TRUE then (nu,xi) are constrined to give valid implied (c,p) values.
#' @return A list consisting of
#'   \item{params}{A vector containing the estimated parameters, in the order (mu,C,alpha,nu,xi,beta)}
#'   \item{loglik}{The corresponding loglikelihood}

#' @examples
#'\dontrun{
#' beta <- 2.4; M0 <- 3; maxTime <- 500
#' catalog <- simulateETAS(0.2, 0.2, 1.5, 0.5, 2, beta = beta,M0 = M0,maxTime = maxTime)
#' maxLikelihoodETAS(catalog$ts, catalog$magnitudes, M0, 500)
#'}
#' @author  Gordon J. Ross & Zak Varty
#' @references Gordon J. Ross - Bayesian Estimation of the ETAS Model for Earthquake Occurrences (2016)
#' @export
maxLikelihoodETAScGPD <- function(ts,magnitudes,M0,maxTime=NA,initval=NA,displayOutput=TRUE, constrainOmori=FALSE){
  # check inputs
  if (is.na(maxTime)) {maxTime <- max(ts)}
  marks <- magnitudes
  if (any(marks < M0)){ stop('Error: All magnitudes must be at least M0.')}

  # Calcualte lag matrix
  L <- outer(X = ts,Y = ts,FUN = "-")
  L <- L * lower.tri(L)

  # Calculate mean magnitude
  Mbar <- mean(marks)

  # generic llh function
  markedHawkeslogLhood <- function(L, marks, cumint, mu, kappa, h, H, maxTime) {
    n <- length(ts)
    if (length(mu)==1) {mu <- rep(mu,n)}
    if (length(mu) != n) {print("mu wrong length"); return(NA)}


    llh <- log(mu + h(L) %*% kappa(marks))
    llh <- sum(llh)
    llh <- llh - cumint
    llh <- llh - sum(kappa(marks) * H(maxTime-ts))
    llh
  }

  # nllh as function of params
  nllh_fn <- function(params) {
    # separate parameters
    mu <- params[1]
    C <- params[2]
    alpha <- params[3]
    nu <- params[4]
    xi <- params[5]

    #check parameter validity
    if (mu <= 0 || nu/(1+xi) <= 0 || alpha < 0 || C < 0) {return(Inf)} #need more than these
    if (constrainOmori){
      if(xi <= 0 || nu<=0) {return(Inf)} #ensures implied p>1 and c>0.
    }
    cumint <- mu*maxTime

    #define component functions for current parameters
      # additional intensity function (takes numeric or vector)
      h <- function(z) {
        dgpd(x = z, shape = xi, nu = nu, mu = 0)
      }
      # additional intensity cdf (takes numeric or vector)
      H <- function(z) {
        pgpd(q = z, shape = xi, nu = nu, mu = 0)
      }
      # rescaling of each peak (takes vector )
      kappa  <- function(m) {
        res <- length(m)
        res[m<M0] <- 0
        inds <- which(m>=M0)
        if (length(inds)>0) {
          res[inds] <- C * exp(alpha*(m[inds]-Mbar))}
        return(res)
      }
    # calculate negative log-likelihood
    nllh <- -markedHawkeslogLhood(L, marks, cumint, mu, kappa, h, H, maxTime)
    if (displayOutput==TRUE) {print(c(mu,C,alpha,nu,xi,nllh))}
    nllh
  }

  #--- Optimisation of nllh_fn

  if (is.na(initval[1])) {
   # initval <- c(length(ts)/maxTime,0.5,0.5,1,2)     # power-law starting vals
    C0 <- 0.5 * exp(0.5*(mean(marks)-M0))
    initval <- c(length(ts)/maxTime, C0, 0.5, 2, 1)   #cGPD starting vals
  } else if (length(initval)==6) {
    initval <- initval[1:5]
  }

  opt <- optim(initval, nllh_fn, control=list(trace=FALSE))

  #now add the GR estimator....
  beta <- 1/mean(marks-M0)
  opt$par <- c(opt$par, beta)
  opt$value <- opt$value-sum(dexp(marks-M0,beta,log=TRUE))

  return(list(params=opt$par, loglik= -opt$value))
}

