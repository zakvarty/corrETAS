% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampleETASposterior_dualmag.R
\name{sampleETASposterior_dualmag}
\alias{sampleETASposterior_dualmag}
\title{Draws samples from the posterior distribution of the cGPD ETAS model with two independent magnitude distributions.}
\usage{
sampleETASposterior_dualmag(
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
  mag_sds = rep(0.1, 4),
  mu_prior = c(0.1, 0.1)
)
}
\arguments{
\item{ts}{Vector containing the earthquake times}

\item{ms}{Vector containing the earthquake magnitudes}

\item{m_0}{Magnitude of completion. Only required when initval not specified. Is then used to calculate MLE starting point by the direct method.}

\item{t_max}{Length of the time window [0,T] the catalog was observed over. If not specified, will be taken as the time of the last earthquake.}

\item{init_ETAS}{Initial ETAS parameters at which to start the estimation. If specified, should be a vector, with elements (mu, C, a, nu_t, xi_t). If unspecified, the sampler will be initialized at the maximum likelihood estimate of the model parameters}

\item{init_B}{Initial branching vector at which to start the estimation. If unspecified, all events will be labelled as process 0.}

\item{init_mag}{Inital Magnitude parameters at which to start the estimation (nu_0m, xi_0m, nu_1m, xi_1m).}

\item{sims}{Number of posterior samples to draw}

\item{B_samples}{Logical. Return samples of branching vector?}

\item{B_fixed}{Logical. If TRUE, branching vector treated as known. If FALSE it is estimated as part of the MCMC.}

\item{nuxi_t_steps}{Number of MH steps to attempt at each (nu_t,xi_t) update}

\item{Calpha_steps}{Number of MH steps to attempt at each (C,alpha) update}

\item{mag_steps}{Number of MH steps to attempt at each update of the magnitude distributions}

\item{etas_sds}{Vector of standard deviations for Normal transition kernels in (logC, alpha, nu_t, xi_t)}

\item{mag_sds}{Vector of standerd deviations for Normal transition kernels in (nu_m0,xi_m0,nu_m1,xi_m1)}

\item{mu_prior}{Vector of parameters (alpha, beta) for the gamma prior distribution on mu}
}
\description{
Draws samples from the posterior distribution of the cGPD ETAS model with two independent magnitude distributions.
}
\author{
Zak Varty
}
