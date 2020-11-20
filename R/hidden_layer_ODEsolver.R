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
    indegree_zero <- apply(W,2,function(x) (sum(x)==S_val)) #nodes with no incoming edge(for setting Q)
    diag(W)= S_val # set back diag with alpha value otw diag of inhibited node is zero

    Q = stimulation( parms[[3]],
                     t,
                     length(parms[[1]]))

    #dx/dt of nockdown with big negative Q(t)
    if(any(r!=0) && any(Q<0)){

       Q[r==0]=0#just stimulation of direct target is non-zero
       dx = phi(t(W) %*% x ) + Q - (x * diag(W)) #propagation
       list(dx)

    }
    #dx/dt of experiment of ideal inhibition, stimulation function can e added to system
    if(any(r!=0) && any(Q>=0)){

      Q[!indegree_zero]=0#just stimulation of source node is nonzero
      Q[r!=0]=0 #if a node is inhibited, stimulation of that is zero
      dx =  phi(t(W) %*%   x ) +Q- (x * diag(W)) #propagation #Q and W inside or outside?

    }
    #dx/dt of wild type
    if(all(r==0)){

      Q[!indegree_zero]=0 #just stimulation of source node is nonzero
      dx =   phi( t(W) %*%  x ) +Q- (x * diag(W)) #propagation
    }

    list(dx)

  })
}
