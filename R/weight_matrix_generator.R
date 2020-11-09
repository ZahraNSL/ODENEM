#' W generator
#'
#' @param alpha confidence of edge existance
#' @param prior.W list of rho, lambda, nu, tau = hyperparameters of prior over W
#' @param small_vall diag value
#'
#' @return
#' @export
weight_matrix_generator = function(alpha,
                                   prior.W,
                                   small_val) {
  W <- alpha * matrix(
    rnorm(nrow(alpha),
          mean = prior.W$nu,
          sd = prior.W$tau),
    nrow = nrow(alpha),
    ncol = ncol(alpha)
  ) +
    ((1 - alpha) * matrix(0, nrow = nrow(alpha),
                          ncol = ncol(alpha)))

  diag(W) <- small_val
  return(W)
}
