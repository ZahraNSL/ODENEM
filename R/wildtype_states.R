#' WildeType State
#'
#' @param model list :W, Q, prior.theta, prior.W, orig.dat, r0
#' @param times number of times points to simulate
#' @param ts time steps length
#' @param flag flag to plot states
#'
#' @return
#' @export
wildtype_states <- function(model,
                            r0, #vector of starting points
                            times,
                            ts,
                            flag) {
  out = compute_states(model$W,
                       rep(0,ncol(model$Q)),
                       model$r0,
                       times,
                       method,
                       ts,
                       flag)

  return(out)

}
