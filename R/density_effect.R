#' density of effect expected
#'
#' @param d a numeric quantile
#' @param gumb_loc location parameter of extereme distribution
#' @param gumb_scl scale parameter of extereme distribution
#' @param gumb_shp shape parameter of extereme distribution
#'
#' @return
#' @export
density_effect <- function(d,
                           gumb_loc,
                           gumb_scl,
                           gumb_shp) {
  return(fExtremes::dgev(
    x = d,
    mu = gumb_loc,
    beta = gumb_scl,
    xi = gumb_shp
  ))
}
