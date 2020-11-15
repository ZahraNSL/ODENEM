#' call to ODE solver
#'
#' @param W wighted matrix
#' @param q perturbation vector (1 = perturbed, 0 = unperturbed)
#' @param r0 initial state of unperturbed system
#' @param times number of time points
#' @param method which case (SteadyState or TimeSeries)
#' @param ts time step length
#' @param flag flag to plot states
#' @param stimulation.value constant value for stimulation function
#'
#' @return matrix of states for given experiment
#' @export
compute_states = function(W,
                          q,
                          r0,
                          times,
                          method,
                          ts,
                          flag,
                          stimulation.value) {
  nstates = ncol(W)

  #state matrix T * S
  states = matrix(
    0,
    nrow = times,
    ncol = nstates,
    dimnames = list(paste0('T_', as.character(1:times)),
                    paste0(colnames(W)[1:nstates]))
  )

  if (method == "timeSeries") {
    #ODE approach

    out <- tryCatch({

      if(any(q!=0)){

        if(class(r0)=="numeric"){

          r0[q == 1] = 0

          states <- deSolve::ode(
            y = r0 ,
            times = seq(1, times, ts),
            func = hidden_layer_ODEsolver,#ODE
            parms = list(q, W, stimulation.value),
            input = stimulation # Q(t)
          )

          }else{

            r0[1,q == 1] = 0

            states <- deSolve::ode(
              y = r0[1,] ,
              times = seq(1, times, ts),
              func = hidden_layer_ODEsolver,#ODE
              parms = list(q, W, stimulation.value),
              input = stimulation # Q(t)
            )

            }
      }else{

        states <- deSolve::ode(
          y = r0 ,
          times = seq(1, times, ts),
          func = hidden_layer_ODEsolver,#ODE
          parms = list(q, W, stimulation.value),
          input = stimulation # Q(t)
        )

        }

    #plot state of the dynamics
      if(flag){


        df_states=reshape2:::melt.matrix(states[,-1])

        colnames(df_states) = c("Time","Hidden_node","State_Value")
        Hidden_nodes=as.factor(df_states$Hidden_node)
        Times=as.factor(df_states$Time)
        qual_col_pals = RColorBrewer::brewer.pal.info[
          c("Set1","Dark2"),]
        col_vector = unlist(mapply(RColorBrewer::brewer.pal,
                                   qual_col_pals$maxcolors,
                                   rownames(qual_col_pals)))


         p=ggplot2::ggplot(df_states,
                          ggplot2::aes(x = Times,
                                       y = State_Value,
                                       fill=Hidden_nodes,
                                       shape=Hidden_nodes,
                                       color=Hidden_nodes)) +
          ggplot2::geom_point(alpha = 1,
                              size = 2) +
           #theme_minimal() +
           ggplot2::theme(legend.position = "top",
                 #aspect.ratio=1,
                 axis.text=ggplot2::element_text(size=6),
                 axis.title=ggplot2::element_text(size=6),
                 plot.title = ggplot2::element_text(size=10),
                 legend.title = ggplot2::element_text(size = 6),
                 legend.text = ggplot2::element_text(size = 6))+
           ggplot2::ggtitle(paste0("Experiment_"
                   ,ifelse(any(q!=0),yes = which(q==1),
                           no = "Wild Type States"))) +
           #ggplot2::geom_jitter()+
           ggplot2::scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,
                                                  9,10,11,12,13,14,15))+
           ggplot2::scale_color_manual(values=col_vector)+
           ggplot2::scale_fill_manual(values=col_vector)

        print(p)
      }

    states = states[,-1]
    },
    warning = function(cond) {
      #message(paste("non converged situation"))
      #message(cond)
      # Choose a return value in case of warning
      return(NULL)

    })

    return(states)
  }
  else if (method == "difference") {
    # linear difference equation

    times = 1:length(times) # just indices of time steps
    ###?
  }
  else if (method == "steadyState") {
    # steady state case

    out <- tryCatch({
      states = rootSolve::runsteady(y =r0,
                         fun = hidden_layer_ODEsolver,
                         parms = list(q,W,stimulation.value),
                         time = c(1, Inf),# run system long enough #Inf being tested
                         input=stimulation)$y

      #print state of the wt if it is steady state method
      if(flag){
        if(all(q==0)){print("Wild Type States")}else{
          print(paste0("Experiment_"),which(q==1))
        }

        print(states)

      }


    },
    warning = function(cond) {
      #message(paste("non converged situation"))
      #message(cond)
      # Choose a return value in case of warning
      return(NULL)
    })

    return(states)
  }
  return(states)

}
