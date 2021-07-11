
#===============================================================================
#' Sample uniformly from the set of branching vectors of length n.
#'
#' @param n integer. Number of elements in branching vector.
#'
#' @return Vector of length n that is a valid branching vector.
#' @export
#'
#' @examples
#' rbranch(n = 100)
rbranch <- function(n){
  out <- rep(NA_integer_, n)
  for (i in seq_along(out)){
    out[i] <- sample(x = 0:(i-1), size = 1)
  }
  return(out)
}
#===============================================================================


#===============================================================================

#===============================================================================
