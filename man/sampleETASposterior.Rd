% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampleETASposterior.R
\name{sampleETASposterior}
\alias{sampleETASposterior}
\title{Draws samples from the posterior distribution of the cGPD ETAS model}
\usage{
sampleETASposterior(
  ts,
  ms,
  m_0,
  t_max = NULL,
  init_ETAS = NULL,
  init_B = rep(0, length(ts)),
  init_mag = NULL,
  sims = 5000,
  B_samples = FALSE,
  B_fixed = FALSE,
  nuxi_t_steps = 500,
  Calpha_steps = 100,
  mag_steps = 100,
  etas_sds = rep(0.1, 4),
  mag_sds = rep(0.1, 2),
  mu_prior = c(0.1, 0.1)
)
}
\arguments{
\item{ts}{Vector containing the earthquake times}

\item{ms}{Vector containing the earthquake magnitudes}

\item{m_0}{Magnitude of completion. Only required when initval not specified. Is then used to calculate MLE starting point by the direct method.}

\item{t_max}{Length of the time window [0,T] the catalog was observed over. If not specified, will be taken as the time of the last earthquake.}

\item{init_ETAS}{Initial ETAS parameters at which to start the estimation. If specified, should be a vector, with elements (mu, K, alpha, c, p). If unspecified, the sampler will be initialized at the maximum likelihood estimate of the model parameters}

\item{init_B}{Initial branching vector at which to start the estimation. If unspecified, all events will be labelled as process 0.}

\item{init_mag}{Inital Magnitude parameters at which to start the estimation.}

\item{sims}{Number of posterior samples to draw}

\item{B_samples}{Logical. Return samples of branching vector?}

\item{B_fixed}{Logical. Sample branching vector as part of mcmc or keep it fixed at init_B?}

\item{nuxi_t_steps}{Number of MH steps to attempt at each (nu,xi) update}

\item{Calpha_steps}{Number of MH steps to attempt at each (C,alpha) update}

\item{mag_steps}{Number of MH steps to attempt at each update of the magnitude distribution}

\item{etas_sds}{Vector of standard deviations for Normal transition kernels in (logC, alpha, nu_t, xi_t)}

\item{mag_sds}{Vector of standerd deviations for Normal transition kernels in (nu_m,xi_m)}

\item{mu_prior}{Vector of parameters (alpha, beta) for the gamma prior distribution on mu}
}
\value{
A list consisting of
    \item{etas}{A  matrix containing sampled Omori ETAS parameter values. Each row is a single sample, and the columns correspond to (mu, C, alpha, nu, xi)}
    \item{b}{A matrix containing sampled Branching vectors. Each row is a single sample,and the columns correspond to earthquakes}
    \item{mag}{A matrix containing sampled nu-GPD magnitude parameter values. Each row is a single sample and the columns correspond to (nu,xi)}
}
\description{
This function implements a reparameterised version of the latent variable MCMC scheme from (Ross 2016) which draws samples from the Bayesian posterior distribution of the Epidemic Type Aftershock Sequence (ETAS)  model.
\deqn{}
The ETAS model is widely used to quantify the degree of seismic activity in a geographical region, and to forecast the occurrence of future mainshocks and aftershocks (Ross 2016). The temporal ETAS model is a point process where the probability of an earthquake occurring at time \eqn{t} depends on the previous seismicity \eqn{H_t}{Ht}, and is defined by the conditional intensity function:
\deqn{ \lambda(t|H_t) = \mu + \sum_{t[i] < t}  \kappa(m[i]|C,\alpha) h(t[i]|\nu,\xi)}{ \lambda(t|Ht) = \mu + \sum \kappa(m[i]|C,\alpha) h(t[i]|\nu,\xi)}
where
\deqn{\kappa(m_i|C,\alpha) = Ce^{\alpha \left( m_i-\bar m \right)}}{\kappa(m[i]|C,\alpha) = C * exp(\alpha(m[i]- mean(m)))}
and \eqn{ h(t_i|\nu,\xi)} is the generalised Pareto density function reparameterised with \eqn{\nu = \sigma(1+\xi)}.
\deqn{}
The summation is over all previous earthquakes that occurred in the region, with the i'th such earthquake occurring at time \eqn{ t_i}{t[i]} and having magnitude \eqn{ m_i}{m[i]}. The quantity \eqn{M_0}{M0} denotes the magnitude of completeness of the catalog, so that \eqn{m_i \geq  M_0}{m[i] \ge M0} for all i. The temporal ETAS model has 5 parameters: \eqn{\mu} controls the background rate of seismicity, \eqn{ C} and \eqn{ \alpha} determine the productivity (average number of aftershocks) of an earthquake with magnitude \eqn{m}, and \eqn{\nu} and \eqn{ \xi}  are the parameters of the generalised Pareto distribution and represent the speed at which the aftershock rate decays over time. Each earthquake is assumed to have a magnitude which is an independent draw from the Gutenberg-Richter law \eqn{ p(m_i) = \beta e^{\beta(m_i-M_0)}}{ p(m) = \beta * exp(\beta(m-M0)}.
}
\author{
Zak Varty
}
