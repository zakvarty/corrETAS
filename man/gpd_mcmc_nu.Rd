% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gpd_mcmc_nu.R
\name{gpd_mcmc_nu}
\alias{gpd_mcmc_nu}
\title{Independent Metropolis sampler for nu-parameterised generalised Pareto distribution.}
\usage{
gpd_mcmc_nu(
  x,
  threshold,
  init,
  step_sds,
  n_samples,
  prior = corrETAS::gpd_flat_prior,
  verbose = FALSE
)
}
\arguments{
\item{x}{vector of observed data}

\item{threshold}{scalar threshold value}

\item{init}{initial values of parameters (nu, xi)}

\item{step_sds}{standard deviations of normal proposal density for each of (nu, xi)}

\item{n_samples}{integer number of samples to draw from the Markov chain}

\item{prior}{function like corrETAS::gpd_flat_prior, giving the joint proir on (nu,xi)}

\item{verbose}{integer/logical indicating how frequently progress and acceptance probability should be printed to the console.}
}
\value{
a data.frame of sampled parameter values
}
\description{
Independent Metropolis sampler for nu-parameterised generalised Pareto distribution.
}
\author{
Zak Varty
}
