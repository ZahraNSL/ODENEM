#' compute log-prior
#'
#' @param W wights of hidden structure
#' @param prior.W list of rho, lambda, nu, tau = hyperparameters of prior over W
#'
#' @return
#' @export
log_prior = function(W,
                     prior.W) {


  # seperate upper lower triangles from diag for W, rho
  W <- W[lower.tri(W) | upper.tri(W)]

  prior.W$rho <-
    prior.W$rho[lower.tri(prior.W$rho) | upper.tri(prior.W$rho)]

  logp <-   sum(log((prior.W$rho) +
                      (
                        (1 - prior.W$rho) *
                          LaplacesDemon::dlaplace(
                            x = W,
                            location = 0, #mu
                            scale = prior.W$rate #1 in case of high peak
                          )
                      )
                    ))
  return(logp)
}
