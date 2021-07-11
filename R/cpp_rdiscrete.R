#' Sample from a discrete distribution in C++
#'
#' @param probs vector of non-negative sampling weights
#' @param n number of samples to draw from distribtuion
#'
#' @return n samples from discrete distribution
#' @export
#'
#' @examples cpp_rdiscrete(probs = c(1,2.5,5), n = 20)
cpp_rdiscrete <- function (probs, n){
  N <- length(probs)
  out <- rep(NA, n)

  for (i in 1:n){
    sample <- -2
    res <- .C('rdiscreteC',
              probs = as.double(probs),
              n = as.integer(N),
              sample = as.integer(sample),
              PACKAGE = 'corrETAS')
    out[i] <- res$sample
  }
  return(out)
}
