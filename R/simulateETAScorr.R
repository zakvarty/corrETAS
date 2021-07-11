#' Draw aftershock magnitudes from a distribution that is dependent on the mainshock magnitude.
#'
#' @param M_main Magnitude of main shock
#' @param n_after Number of aftershock magnitudes to draw
#' @param nuxi_main Parameters of the marginal mainshock GPD.
#' @param nuxi_after Parameters of the marginal aftershock GPD.
#' @param rho Correlation coefficient of gaussian copula.
#' @param mu_main GPD threshold for mainshocks, must be smaller than M_main.
#' @param mu_after GPD threshold for aftershocks.
#' @return A vector of sampled aftershock magnitudes.
#' @export
#'
#' @examples
#' nuxi_main = c(2,0.2)
#' nuxi_after = c(0.5,0)
#' draw_M_after(M_main = 2.12, n_after = 30, nuxi_main = nuxi_main, nuxi_after = nuxi_after, rho =0.7)
draw_M_after <- function(M_main, n_after, nuxi_main, nuxi_after, rho = 0, mu_main = 0, mu_after = 0){
  mg <- M_main - mu_main
  mu <- pgpd(q = mg, shape = nuxi_main[2],nu = nuxi_main[1])
  mn <- qnorm(p = mu)

  an <- rnorm(n = n_after,mean = rho*mn, sd = sqrt(1-rho^2))
  au <- pnorm(q = an)
  ag <- qgpd(p = au,shape = nuxi_after[2],nu = nuxi_after[1])

  return(ag + mu_after)
}



#'Simulates synthetic data from the ETAS model
#'

#' This function simulates sample data from the ETAS model over a particular interval [0,T].

#'The Epidemic Type Aftershock Sequence (ETAS) model is widely used to quantify the degree of seismic activity in a geographical region, and to forecast the occurrence of future mainshocks and aftershocks (Ross 2016). The temporal ETAS model is a point process where the probability of an earthquake occurring at time \eqn{t} depends on the previous seismicity \eqn{H_t}{Ht}, and is defined by the conditional intensity function:

#' \deqn{ \lambda(t|H_t) = \mu + \sum_{t[i] < t}  \kappa(m[i]|K,\alpha) h(t[i]|c,p)}{ \lambda(t|Ht) = \mu + \sum \kappa(m[i]|K,\alpha) h(t[i]|c,p)}

#' where

#' \deqn{\kappa(m_i|K,\alpha) = Ke^{\alpha \left( m_i-M_0 \right)}}{\kappa(m[i]|K,\alpha) = K * exp(\alpha(m[i]-M0))}

#' and

#' \deqn{ h(t_i|c,p) = \frac{(p-1)c^{p-1}}{(t-t_i+c)^p}}{ h(t[i]|c,p) = (p-1) * c^(p-1) * (t-t[i]+c)^(-p)}

#' where the summation is over all previous earthquakes that occurred in the region, with the i'th such earthquake occurring at time \eqn{ t_i}{t[i]} and having magnitude \eqn{ m_i}{m[i]}. The quantity \eqn{M_0}{M0} denotes the magnitude of completeness of the catalog, so that \eqn{m_i \geq  M_0}{m[i] \ge M0} for all i. The temporal ETAS model has 5 parameters: \eqn{\mu} controls the background rate of seismicity, \eqn{ K} and \eqn{ \alpha} determine the productivity (average number of aftershocks) of an earthquake with magnitude \eqn{m}, and \eqn{c} and \eqn{ p}  are the parameters of the Modified Omori Law (which has here been normalized to integrate to 1) and represent the speed at which the aftershock rate decays over time. Each mainshock is assumed to have a magnitude which is an independent draw from a Generalised Pareto Distribution, Aftershocks can follow a separate Generalised Pareto Distribution and can be made dependent on the magnitude of the triggering event through the parameter rho.

#' \deqn{}

#' This function simulates sample data from the ETAS model over a particular interval [0,T]. Background and aftershock magnitudes follow marginal generalised pareto distributions linked by a gaussian copula.

#' @param mu Parameter of the ETAS model as described above.
#' @param K Parameter of the ETAS model as described above.
#' @param alpha Parameter of the ETAS model as described above.
#' @param gpd_t Parameter vector (nu_t,xi_t) for generalised Pareto distribution of aftershock delays
#' @param rho Correlation coefficient linking background and triggered magnitudes.
#' @param gpd_m0 Parameter vector (nu_m0, xi_m0) for generalised pareto distribution of background magnitudes.
#' @param gpd_m1 Parameter vector (nu_m1, xi_m1) for generalised pareto distribution of triggered magnitudes.
#' @param m_0 Magnitude of completeness.
#' @param t_max Length of the time window [0,t_max] to simulate the catalog over.
#' @param displayOutput If TRUE then prints the number of earthquakes simulated so far.
#' @return A list consisting of
#'   \item{ts}{The simulated earthquake times}
#'   \item{ms}{The simulated earthquake magnitudes}
#'   \item{b}{The simulated branching structure, where branching[i] is the index of the earthquake that triggered earthquake i, or 0 if earthquake i is a background event}
#' @author Zak Varty
#' @export
simulateETAScorr <- function(mu, K, alpha, gpd_t, rho = 0, gpd_m0, gpd_m1 = NULL, m_0, t_max, displayOutput = TRUE){

  ## Magnitude distribution checks
  #check a mainshock distribution is provided
  stopifnot(gpd_t[1] / (1 + gpd_t[2]) > 0)
  stopifnot(gpd_m0[1] / (1 + gpd_m0[2]) > 0)

  if(is.null(gpd_m1)){gpd_m1 <- gpd_m0}
  stopifnot(gpd_m1[1] / (1 + gpd_m1[2]) > 0)

  nu_m0 <- gpd_m0[1]
  xi_m0 <- gpd_m0[2]
  nu_m1 <- gpd_m1[1]
  xi_m1 <- gpd_m1[2]
  nu_t  <- gpd_t[1]
  xi_t  <- gpd_t[2]

  n <- rpois(n = 1, lambda = mu * t_max)
  x <- runif(n = n, min = 0,max = t_max)
  x <- sort(x)

  m <- rgpd(n = length(x), nu = nu_m0, shape = xi_m0, mu = m_0)
  b <- rep(0,n)

  count <- 1

  while(count < length(x)) { #add new events onto x...
    if (displayOutput==TRUE) {
      print(sprintf("%s events generated so far",count,length(x)))
    }

    if (count > length(x)) {break}
    pt <- x[count]
    pm <- m[count]
    pb <- b[count]

    E_after <- K* exp(alpha*(pm-m_0))
    n_after <- rpois(n = 1,lambda = E_after)

    if (n_after > 0){
      x_after <- pt + rgpd(n = n_after, shape = xi_t, nu = nu_t)
      m_after <- draw_M_after(
        M_main = pm,
        n_after = n_after,
        nuxi_main = gpd_m0 * (pb == 0) + gpd_m1 * (pb !=  0),
        nuxi_after = gpd_m1,
        rho = rho,
        mu_main = m_0,
        mu_after = m_0
      )

      x <- c(x, x_after)
      m <- c(m,m_after)
      b <- c(b, rep(x[count],n_after))
    }
    count <- count + 1
  }

  ord <- order(x)
  x <- x[ord]; m <- m[ord]; b <- b[ord]

  incl <- which(x <= t_max)
  x <- x[incl]; m <- m[incl]; b <- b[incl]

  for (i in 1:length(b)) {
    if (b[i]==0) {next}
    b[i] <- which(x==b[i])
  }
  return(list(ts=x,ms=m,b=b))
}



