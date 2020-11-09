#' density of no effect expected
#'
#' @param d a numeric quantile
#' @param trunc_mean mean of truncated normal
#' @param trunc_sd standard deviation of truncated normal
#'
#' @return
#' @export
density_no_effect = function(d,
                      trunc_mean,
                      trunc_sd) {
  return(msm::dtnorm(x = d,
                mean = trunc_mean,
                sd = trunc_sd))
}
