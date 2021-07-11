#' Draws samples from the posterior distribution of the cGPD ETAS model with two correlated magnitude distributions.
#'
#' @param ts Vector containing the earthquake times
#' @param ms Vector containing the earthquake magnitudes
#' @param m_0 Magnitude of completion. Only required when initval not specified. Is then used to calculate MLE starting point by the direct method.
#' @param t_max Length of the time window [0,T] the catalog was observed over. If not specified, will be taken as the time of the last earthquake.
#' @param init_ETAS Initial ETAS parameters at which to start the estimation. If specified, should be a vector, with elements (mu, C, a, nu_t, xi_t). If unspecified, the sampler will be initialized at the maximum likelihood estimate of the model parameters
#' @param init_B Initial branching vector at which to start the estimation. If unspecified, all events will be labelled as process 0.
#' @param init_mag Inital Magnitude parameters at which to start the estimation (nu_0m, xi_0m, nu_1m, xi_1m, rho).
#' @param sims Number of posterior samples to draw
#' @param B_samples Logical. Return samples of branching vector?
#' @param B_fixed Logical. If TRUE, branching vector treated as known. If FALSE it is estimated as part of the MCMC.
#' @param nuxi_t_steps Number of MH steps to attempt at each (nu_t,xi_t) update
#' @param Calpha_steps Number of MH steps to attempt at each (C,alpha) update
#' @param mag_steps Number of MH steps to attempt at each update of the magnitude distributions
#' @param etas_sds Vector of standard deviations for Normal transition kernels in (logC, alpha, nu_t, xi_t)
#' @param mag_sds Vector of standerd deviations for Normal transition kernels in (nu_m0,xi_m0,nu_m1,xi_m1,rho)
#' @param mu_prior Vector of parameters (alpha, beta) for the gamma prior distribution on mu
#' @param logC_prior Vector of parameters (mu,sigma) for the normal prior distribtuion on logC
#' @author Zak Varty
#' @export
sampleETASposterior_corrmag <- function(ts, ms, m_0, t_max=NULL, init_ETAS=NULL, init_B = rep(0,length(ts)), init_mag=NULL,  sims=5000, B_samples = FALSE, B_fixed = FALSE, nuxi_t_steps=500, Calpha_steps=100, mag_steps=100, etas_sds=rep(0.1,4), mag_sds=rep(0.1,4), mu_prior=c(0.1,0.1), logC_prior = c(0.3,0.01)) {
  #--- input checks
  if (is.null(t_max)) {t_max <- max(ts)}

  if(any(etas_sds < 0, !is.numeric(etas_sds))){
    stop("Error in `etas_sds`: Standard deviations must be positive and numeric.")
  }

  if(any(mag_sds < 0, !is.numeric(mag_sds))){
    stop("Error in `mag_sds`: Standard deviations must be positive and numeric.")
  }

  if(any(mu_prior <= 0, !is.numeric(mu_prior))){
    stop("Error in `mu_prior`: Gamma parameters must be positive and numeric. ")
  }

  if(any(ms < m_0)){
    stop("Error in `ms`: All values must be at least m_0.")
  }
  #--- find starting values if missing
  if(is.null(init_mag)){
    cat("No initial magnitude parameters passed. \n  Will intialise as independent.\n Computing maximum likelihood estimate for background events... \n")

    mag0_mle <- maxLikelihoodGPD(
      x = ms[init_B==0],
      nuxi = c(1,0),
      mu = m_0
    )$params
    cat("Computing maximum likelihood estimate for triggered events... \n")
    if( sum(init_B>0) > 0 ){
      mag1_mle <- maxLikelihoodGPD(
        x = ms[init_B>0],
        nuxi = c(1,0),
        mu = m_0
      )$params
    } else {
      cat("No initial triggered events. Will use background parameters. \n")
      mag1_mle <- mag0_mle
    }
    mag_mle <- c(mag0_mle, mag1_mle, 0)
    cat("Maximum likelihood estimate for magnitude parameters assuming independence is: \n")
    print(mag_mle)
    init_mag <- mag_mle
    cat("Will use this for inital value. \n\n")
  }

  if (is.null(init_ETAS)) {
    cat("No initial ETAS parameters passed. Computing maximum likelihood estimate... \n")
    n <- length(ts)
    t_max_optim <- t_max
    if (n > 500) {
      n <- 500
      t_max_optim <- ts[500]
    }
    init_ETAS <- maxLikelihoodETAScGPD(
      ts = ts[1:n],
      magnitudes = ms[1:n],
      M0 = m_0,
      maxTime = t_max_optim,
      displayOutput= FALSE,
      constrainOmori = FALSE)$params
    cat("Maximum likelihood estimate for ETAS parameters is: \n")
    print(init_ETAS[1:5])
    cat("Will use this for initial value \n\n")
  }

  #--- MCMC setup
  branching <- init_B
  mu    <- init_ETAS[1]
  C     <- init_ETAS[2]
  alpha <- init_ETAS[3]
  nu_t  <- init_ETAS[4]
  xi_t  <- init_ETAS[5]
  nu_mb <- init_mag[1] ## nu_m0 === nu_mb === nu for background magnitudes
  xi_mb <- init_mag[2]
  nu_ma <- init_mag[3] ## nu_m1 === nu_ma === nu for aftershock magnitudes
  xi_ma <- init_mag[4]
  rho   <- init_mag[5]

  B_size <-  B_samples * sims * length(ts)

  cat(sprintf("Starting Gibbs sampler. Will draw %s posterior samples... \n",sims))
  out <- NULL

  #--- MCMC Exported to C++ via C
    res <- .C("estimateETAScorr_C",
              ts=as.double(ts),
              marks=as.double(ms),
              branching=as.integer(branching),
              n=as.integer(length(ts)),
              M0 = as.double(m_0),
              maxTime = as.double(t_max),
              sims=as.integer(sims),
              Bsamples=as.integer(B_samples),
              Bsize=as.integer(B_size),
              Bfixed=as.integer(B_fixed),
              nuxiSteps=as.integer(nuxi_t_steps),
              CaSteps=as.integer(Calpha_steps),
              magSteps=as.integer(mag_steps),
              mu=as.double(mu),
              logC=as.double(log(C)),
              alpha=as.double(alpha),
              nut=as.double(nu_t),
              xit=as.double(xi_t),
              numb=as.double(nu_mb),
              ximb=as.double(xi_mb),
              numa=as.double(nu_ma),
              xima=as.double(xi_ma),
              rho=as.double(rho),
              etasstepSds=as.double(etas_sds),
              magstepSds=as.double(mag_sds),
              muprior=as.double(mu_prior),
              logCprior=as.double(logC_prior),
              constrainOmori = as.integer(FALSE),
              mus=as.double(numeric(sims)),
              logCs=as.double(numeric(sims)),
              alphas=as.double(numeric(sims)),
              nuts=as.double(numeric(sims)),
              xits=as.double(numeric(sims)),
              numbs=as.double(numeric(sims)),
              ximbs=as.double(numeric(sims)),
              numas=as.double(numeric(sims)),
              ximas=as.double(numeric(sims)),
              rhos=as.double(numeric(sims)),
              Bs=as.integer(numeric(B_size)),
              PACKAGE="corrETAS")

  #--- Formatting samples into output
    etas_samples <- data.frame(
      mu =res$mus,
      C =res$logCs,
      alpha = res$alphas,
      nu_t = res$nuts,
      xi_t = res$xits)
    branching_samples <- matrix(res$Bs, ncol = length(ts), byrow = TRUE)
    mag_samples <- data.frame(
      nu_m0 = res$numbs,
      xi_m0 = res$ximbs,
      nu_m1 = res$numas,
      xi_m1 = res$ximas,
      rho = res$rhos)

    etas_samples[,2] <- exp(etas_samples[,2])

    out <- list(
      etas = etas_samples,
      mag = mag_samples,
      b = branching_samples)

  return(out)
}

