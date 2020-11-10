#' constructing Model Object
#'
#' @param stype name of motif : "Diamond","Bifan","FeedForward","FeedBackward","ProteinCascade"
#' @param ti_num number of time points
#' @param prior.theta prior over observations and hidden nodes connections
#' @param prior.W prior over W
#' @param data.para list of parameter for data simulation: mu, sd, location, scale shape
#' @param W.factor W=1 * W.factor
#' @param depletion.factor diag of W
#' @param depletion.factor.wt diag of W in WT
#' @param init.ODE initial value of ODE
#' @param E_num number of observations
#' @param method.wt method for WildType state calculation: steadyState or timeSeries
#' @param stimulation.value constant value for stimulation function
#' @param stimulation.value.wt constant value for stimulation function in WildType case
#' @param tstep time step length
#'
#' @return
#' @export
model_generator <- function(stype,
                            ti_num,
                            prior.theta,
                            prior.W,
                            data.para,
                            W.factor,
                            depletion.factor,
                            depletion.factor.wt,
                            init.ODE,
                            E_num,
                            method.wt,
                            stimulation.value,
                            stimulation.value.wt,
                            tstep) {
  load(paste0(stype, ".RData")) #load the true structure of the Motif
  igraph::plot.igraph(g)#vertex.size=30

  #make matrix from true network graph
  #g = igraph object of the Motif------------------------
  mat_g = as.matrix(
    igraph::as_adjacency_matrix(
      g,
      type = c("both", "upper", "lower"),
      attr = NULL,
      edges = FALSE,
      names = TRUE,
      sparse = igraph::igraph_opt("sparsematrices")
    )
  )


  #nodes
  S_num = nrow(mat_g)#hidden nodes
  P_num = S_num #perturbations node


  #NULL Object
  model = list(
    W = NULL,
    #true_W
    Q = NULL,
    #Experiments set
    prior.theta = prior.theta,
    #probability of attachements between Obs-nodes and Hidd-nodes
    prior.W = prior.W,
    orig.dat = NULL,
    r0 = NULL,
    data.para = data.para,
    method.wt = method.wt,
    depletion.factor.wt = depletion.factor.wt,
    stimulation.value = stimulation.value,
    stimulation.value.wt = stimulation.value.wt

  )



  #initialize elements------------------------------

  #Experiments set------------------------------
  #matrix #hidden nodes * #Perturbation nodes
  #1 in each row is the direct targe of perturbaton
  model$Q <- matrix(
    0,
    nrow = S_num,
    ncol = P_num,
    dimnames = list(rownames(mat_g),
                    rownames(mat_g))
  )
  diag(model$Q) = 1 #single pert in each experiment



  #prior probabilities over S_E connections------------------------------
  if(is.null(prior.theta)){

    #for Simulation (1:E_num to S_1) ,(E_num+1 : 2*E_num to S_2), ...
    X = diag(x = 1, nrow = S_num)
    sim.theta = do.call(rbind, lapply(1:S_num ,
                                      function(s)
                                        rbind(
                                          matrix(
                                            X[s, ],
                                            nrow = E_num/S_num,
                                            ncol = S_num ,
                                            byrow = TRUE
                                          )
                                        )))

    #for true_model (later being passed for learn process)
    #we assume there is no prior over SE connection---> 1/S_num
    Y = matrix(data = (1 / S_num),
               nrow = E_num ,
               ncol = S_num)
    model$prior.theta = t(apply(Y, 1, function(x)
      x / sum(x)))

  }else{

    sim.theta = prior.theta
  }



  # parameters for W ------------------------------
  if (is.null(prior.W)) {
    model$prior.W$rate = 1    #scale parameter for dlaplace, Prior = log(sum(rho+ (1-rho)* dlaplace(0,rate)))
    model$prior.W$nu = 1       #Mu parameter for rnorm function(W!=0)
    model$prior.W$tau = 0.25    #sd parameter for rnorm function(W!=0)
    model$prior.W$rho = matrix(0.5, S_num, S_num)#it should be changed in case of prior knowledge of edges
  }

  #set W------------------------------
  model$W <- as.matrix(igraph::get.adjacency(g))[1:S_num, 1:S_num] *
    #matrix(data = runif(S_num*S_num,0.5,1),S_num,S_num) *
    W.factor
  diag(model$W) <- depletion.factor


  #start states of hidden nodes------------------------------
  if (length(init.ODE) == 1) {
    model$r0 = rep(init.ODE, S_num) #start vector of WT and Experiments,can change.but should not effect
  } else if (length(init.ODE) > 1) {
    model$r0 = init.ODE
  }

  #ODE method SteadyState or TimeSeries?----------------------------------
  if (ti_num == 1){
    method = "steadyState"
  } else{
    method = "timeSeries"
  }

  #Log Fold changes parameters----------------------------------
  if (is.null(data.para)){
    model$data.para$trunc.mean<-0
    model$data.para$trunc.sd<- 0.15
    model$data.para$gumb.loc<-1.6
    model$data.para$gumb.scl <- 0.4
    model$data.para$ gumb.shp <- 0.4

  }
  #simulate data------------------------------
  model$orig.dat = data_simulation(
    model = model,
    nE = E_num ,
    times = ti_num ,
    theta = sim.theta,
    ts = tstep,
    method = method,
    method.wt = method.wt,
    depletion.factor.wt = depletion.factor.wt,
    stimulation.value = stimulation.value,
    stimulation.value.wt = stimulation.value.wt
  )



  return(model)

}
