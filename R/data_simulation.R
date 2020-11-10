#' simulate data for a model
#'
#' @param model list :W, Q, prior.theta, prior.W, orig.dat, r0
#' @param nE number of observations to simulate
#' @param times number of times points to simulate
#' @param theta S_E connection prior
#' @param ts time steps length
#' @param method time series or steady state
#' @param depletion.factor.wt diag of W in WT
#' @param stimulation.value constant value for stimulation function
#' @param stimulation.value.wt constant value for stimulation function in WildType case
#'
#' @return data array (nE x #number of experiments x times)
#' @export
data_simulation = function(model,
                           nE,
                           times,
                           theta,
                           ts,
                           method,
                           method.wt,
                           depletion.factor.wt,
                           stimulation.value,
                           stimulation.value.wt
                           ) {
  npert = NROW(model$Q) #number of P nodes
  R = list() #output of ODE in timeseries matrix with dimension S_nodes* time (list of matrices for report?)
  gamma_val = list() # activity levels based on R matrix, dimension S_nodes* time

  # strcture of observation data = array with nE*nPert*times dimensions
  O = array(
    0,
    dim = c(nE, npert, length(seq(1, times, ts))),
    dimnames = list(
      paste0('E_', as.character(1:nE)),
      paste0('pert_', as.character(1:npert)),
      paste0('T_', 1:length(seq(1, times, ts)))
    )
  )

  if(!is.null(method.wt)){
    W=model$W
    diag(W) = depletion.factor.wt
    R0 = compute_states(W = W,
                        q = rep(0,ncol(model$Q)),
                        r0 = model$r0,
                        times = times,
                        method = method.wt,
                        ts = ts,
                        flag =TRUE,
                        stimulation.value = stimulation.value.wt)

  }else{
    R0=model$r0
  }

  # draw O[i,q,t] from mixture distribution
  for (q in 1:npert) {
    # for a given true model(with W, starting points,...), using corresponding ODE function to solve ODE

  R = compute_states(model$W,
                       model$Q[q,],
                       R0,
                       times,
                       method,
                       ts,
                       flag =TRUE,
                       stimulation.value)
       #Activety level of hidden node = gamma function = tanh|x(t)-x(WT)|

    gamma_val = gammas(R = R,
                       R_0 = R0,
                       flag = TRUE)
     # draw O[i,q,t] from mixture distribution (Eq. 3)
    for (i in 1:dim(O)[1]) {
      # iterate over all E-nodes
      for (t in 1:length(seq(1, times, ts))) {
        # iterate over times points

        # indicator is a variable with bionomial dist. for each observation regarding
        # activity level of parent in hidden nodes (q= gamma_val(t,which(theta[i,]==1)))
        indicator = sample(
          x = c(FALSE, TRUE),
          size = 1,
          prob = c(1 - gamma_val[t, which(theta[i,] ==
                                            1)],
                   gamma_val[t, which(theta[i,] == 1)])
        )
        # simulating observation from mixture of two densities
        #observation are abs(LFC)
        O[i, q, t] = ifelse(indicator,
                              abs(random_generator_effect(model$data.para$gumb.loc,
                                                        model$data.para$gumb.scl,
                                                        model$data.para$gumb.shp
                              # model$data.para[3],
                              #                           model$data.para[4],
                              #                           model$data.para[5]
                                                        )),
                            abs(random_generator_no_effect(model$data.para$trunc.mean,
                                                           model$data.para$trunc.sd
                              # model$data.para[1],
                              #                              model$data.para[2]
                              ))) #O[i,q,t]= indicator in binary case

      }

    }

  }

  return(O)
}
