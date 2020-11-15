#' main ODE solver
#'
#' @param t time
#' @param x ODE answer
#' @param parms list of input parameters: experiment vector,W
#' @param input stimulation function
#'
#' @return dynamic of hidden nodes under special experiment
#' @export
hidden_layer_ODEsolver <- function(t,
                                   x,
                                   parms,
                                   input) {
  with(as.list(c(parms, x)), {

    #####################################
    #tanh(t(W) %*% x_t + Q_t) - alpha*x_t
    #####################################


    W <- parms[[2]]
    S_val <- diag(W)[1]

    r <-as.integer(parms[[1]] != 0) # direct target of experiment
    W[, r != 0] <- 0 #zero indegree for inhibition receptors
    #diag(W)= S_val # set back diag with alpha value

    Q = stimulation( parms[[3]],
                     t,
                     length(parms[[1]]))

    #ndown with Q(t)
    if(any(r!=0) && any(Q<0)){
       Q[r==0]=0#rep(0,ncol(W))
       dx = phi(t(W) %*% x + Q) - (x * diag(W)) #propagation
       list(dx)

    }
    #in stimulation
    if(any(r!=0) && any(Q>=0)){
      Q[r!=0]=0#rep(0,ncol(W))
      x[r != 0] <- 0 #stable state for inhibited targets
      dx = phi(t(W) %*% x + Q) - (x * diag(W)) #propagation
      dx[r != 0] <- 0 #stable state for inhibited targets

    }
    if(all(r==0)){

      dx = phi(t(W) %*% x + Q) - (x * diag(W)) #propagation
    }


    list(dx)

  })
}
