#' learn W from data
#'
#' @param dat feeded data of correspondings experiments,dimension nObs.nExp.times
#' @param Q matrix of perturbation vectors Experiments*hidden nodes
#' @param r_0 start states of hidden nodes
#' @param rho prior over edge existance
#' @param hyper.W list of rho, lambda, nu, tau = hyperparameters of prior over W
#' @param prior.theta prior probabilities over attachements E-nodes to S-nodes
#' @param depletion.factor diag of W
#' @param depletion.factor.wt diag of W in WT
#' @param method.wt method for WildType state calculation: steadyState or timeSeries
#' @param stimulation.value constant value for stimulation function
#' @param stimulation.value.wt constant value for stimulation function in WildType case
#' @param seed random seed
#' @param ts Timesteps
#' @param itr number of iterations for MCMC
#'
#' @return returns sampled Ws and LogLiklihoods
#' @export
MCMC_sampling <- function(dat,
                          Q,
                          r_0,
                          rho,
                          hyper.W,
                          prior.theta,
                          data.para,
                          depletion.factor,
                          depletion.factor.wt,
                          method.wt,
                          stimulation.value,
                          stimulation.value.wt,
                          seed = seed,
                          ts,
                          itr) {
  # set parameters
  npert = NROW(Q)
  nstates = NCOL(Q)
  times = dim(dat)[3]
  # tryCatch(as.numeric(dimnames(dat)[[3]]),
  #          error = function(e)
  #            stop("dimnames[[3]] has to indicate measurement time points"))

  lambda = 0.05 #modification scale in kernel of MCMC W=W+lambda

  #method of ODE
  if (times == 1) {
    method = "steadyState"
  } else{
    method = "timeSeries"
  }

  if (is.null(rho))
    rho = matrix(
      0.1,
      ncol = nstates,
      nrow = nstates,
      dimnames = list(colnames(Q), colnames(Q))
    )

  if (is.null(prior.theta))
    prior.theta = matrix(1 / (nstates + 1), nrow = NROW(dat), ncol = nstates)

  if (is.null(hyper.W)) {
    hyper.W$lambda = 10 # scale parameter of exponential distribution
    hyper.W$nu = hyper.W$tau = 10 # shape and rate parameters of gamma distribution
  }
  if (is.null(data.para)) {
    print("define para for densities")
  }


  #generate model_0 with parameters and W_0
  alpha = alpha_val(rho)
  W = weight_matrix_generator(alpha = alpha ,
            prior.W = hyper.W,
            small_val = depletion.factor)

  model = list(
    W = W,
    Q = Q,
    prior.theta = prior.theta,
    prior.W = list(
      rho = rho,
      nu = hyper.W$nu,
      tau = hyper.W$tau,
      rate = hyper.W$rate
    ),
    orig.dat = dat,
    r0 = r_0,
    data.para= data.para,
    method.wt = method.wt,
    depletion.factor.wt = depletion.factor.wt,
    stimulation.value = stimulation.value,
    stimulation.value.wt = stimulation.value.wt
  )

  #lists of parameters result by MCMC chain
  acc_res = list() #accepted Ws
  acc_ratio = c() #accepted ds
  all_res = list() #all sampled W
  all_ratio = c() #all ds


  # calculate and save log-d, log-prior of initial model_0
  ## control oscilation in steady state case
  # in case of oscillation d function returns epsilon
  pre_ll <-  rcppLikelihood(model,
                             compute_states ,
                             gammas ,
                             density_effect,
                             density_no_effect,
                             wildtype_states,
                             stimulation,
                             method,
                             ts
                            )#likelihood(model,method,ts)



  # while(test_ll==-1e-20){
  #   W=true.w(alpha,hyper.W,diag_W)
  #   model$W <- W
  #   pre_ll <- Likelihood(model,dat,compute_states ,
  #                        gammas ,ddens_1,ddens_0,
  #                        method ,ts)  # control oscilation
  # }
  all_ratio[1] = pre_ll  #save initial LL
  all_res[[1]] = W #save initial W
  pre_p = log_prior(model$W, model$prior.W)



  #perform MCMC
  #profvis({ #check speed of learning process (c++ vs R for Likelihood function)
  for (it in 1:itr)
  {
    valid <- 0
    # sample movement
    #insertion, deletion, reverse, sign change,weight modification
    # probability of movements can change
    move_type = sample(
      x = c(1, 3, 2, 4, 5),
      size = 1,
      prob = c(1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5)#c(1/14,1/14,1/14,1/14,5/7)
    )


    #vectorized vector of W(-diag)
    #later be used to check if selected move_type is possible
    #always one of triangles need to be transposed (for checking the reverse case)
    up_down_vectorized = as.vector(c(as.numeric(model$W[lower.tri(model$W,
                                                                  diag = FALSE)] !=
                                                  0),
                                     as.numeric(t(model$W)[lower.tri(t(model$W),

                                                                     diag = FALSE)] !=
                                                  0)))
    #check possibility of move_type on the W
    if (move_type == 1) {
      #insertion

      if (0 %in% up_down_vectorized) {
        #any zero weight in vectorized W?

        valid = 1

        #indexes of W_ij with zero weights
        indxs <-
          which(model$W == 0, arr.ind = TRUE)
        if (class(indxs) == "matrix" &&
            (nrow(indxs) != 0)) {
          indx = sample(1:nrow(indxs), 1)  #sample one non-edge from severals
          i_indx = indxs[indx,][1]
          j_indx = indxs[indx,][2]
        } else if (class(indxs) == "integer") {
          i_indx = indxs[1]
          j_indx = indxs[2]
        }

      } else{
        valid <- 0 #else reject movement
      }


    } else if ((move_type == 2) || #deletion
               (move_type == 4) || #sign *-1
               (move_type == 5)#weight modification
    ) {

      if (any(as.logical(up_down_vectorized))) {
        #any nonzero weight in vectorized W?

        valid <- 1

        #indexes of W_ij with nonzero weights
        indxs <-
          which(model$W != 0, arr.ind = TRUE)
        indxs <-
          indxs[which(indxs[, 1] != indxs[, 2], arr.ind = TRUE),]
        if (class(indxs) == "matrix" &&
            (nrow(indxs) != 0)) {
          indx = sample(1:nrow(indxs), 1)  #sample one edge from several possibles
          i_indx = indxs[indx,][1]
          j_indx = indxs[indx,][2]
        } else if (class(indxs) == "integer") {
          i_indx = indxs[1]
          j_indx = indxs[2]
        }
      } else{
        valid <- 0
      }


    } else {
    #reverse direction

    #if the reverse edge was not exist before
    if (any(xor(up_down_vectorized[1:(length(up_down_vectorized) / 2)] ,
                up_down_vectorized[((length(up_down_vectorized) / 2) +
                                    1):length(up_down_vectorized)]))) {
      valid <- 1

      #indexes of W_ij=0 and W_ji!=0
      indxs0 <-
        which(model$W == 0 , arr.ind = TRUE)
      indxs1 <-
        which(model$W != 0 , arr.ind = TRUE)
      indxs1 <-
        indxs1[which(indxs1[, 1] != indxs1[, 2], arr.ind = TRUE),]
      indxs <-
        rbind(indxs1, cbind(indxs0[, 2], indxs0[, 1]))
      indxs <- indxs[duplicated(indxs),]
      if (class(indxs) == "matrix" &&
          (nrow(indxs) != 0)) {
        indx = sample(1:nrow(indxs), 1)
        i_indx = indxs[indx,][1]
        j_indx = indxs[indx,][2]
      } else if (class(indxs) == "integer") {
        i_indx = indxs[1]
        j_indx = indxs[2]
      }
    } else{
      valid <- 0
    }
  }



  #make new model for alteration
  #apply the valid change
  new_model = model

  if (valid) {
    if (move_type == 1) {
      #insertion of W regarding prior.W

      new_model$W[i_indx, j_indx] = (rnorm(
        1,
        mean  = model$prior.W$nu,
        sd  = model$prior.W$tau
      ))
    } else if (move_type == 2) {
      #deletion

      new_model$W[i_indx, j_indx] = 0

    } else if (move_type == 3) {
      #reverse an edge ji is index of nonedge in old_W

      new_model$W[i_indx, j_indx] = 0
      new_model$W[j_indx, i_indx] = model$W[i_indx, j_indx]

    }  else if (move_type == 4) {
      #change sign of the edge weight

      new_model$W[i_indx, j_indx] = -1 * model$W[i_indx, j_indx]

    }  else {
      #modify edge weight!=0 with lambda

      new_model$W[i_indx, j_indx] = model$W[i_indx, j_indx] +
        rnorm(1, mean = 0, sd = lambda)

    }


    # logL logP of new_model in case of valid movements
    prop_ll =  rcppLikelihood(new_model,
                              compute_states ,
                               gammas,
                               density_effect,
                               density_no_effect,
                               wildtype_states,
                               stimulation,
                               method ,
                               ts)#likelihood(new_model,method,ts)

    prop_p = log_prior(new_model$W,
                   new_model$prior.W)
    # reject new_W in case of Oscilation
    #stability_chk <-c()
    # if(method == "S"){
    #   # for(q in 1:nrow(Q)){
    #   ###################check diagonalization##############
    #   # rout<- runsteady(y = new_model$r0, fun = rfunc,
    #parms = list(new_model$Q[q,],new_model$W),
    #   #                     time = c(1,100),input=input)$y
    #
    #
    #   Q_star <- new_model$Q#[q,]
    #   Q_star [Q_star==Inf] <- 0
    #   stability_chk <- jacMat(x_star = steady_S,
    #                           nstates = nstates,
    #                           W = new_model$W,
    #                           Q = Q_star,
    #                           depl_const = 0.1666667)
    #
    #   #}
    # }
    #



    #Vectorize new_W
    up_down_vectorized_new  = as.vector(c(as.numeric(new_model$W[lower.tri(new_model$W,
                                                                           diag = FALSE)] !=
                                                       0),
                                          as.numeric(t(new_model$W)[lower.tri(t(new_model$W),
                                                                              diag = FALSE)] !=
                                                       0)))

    # count number of possible structure you can reach from W_old and W_new
    #number of possible reverse edge direction
    #from old_W
    rev_N_old <-
      sum(xor(up_down_vectorized[1:(length(up_down_vectorized) / 2)] ,
              up_down_vectorized[((length(up_down_vectorized) /
                                     2) + 1):length(up_down_vectorized)]))
    N_old <- sum(model$W == 0) +
      (sum(model$W != 0) - nstates) + #diag is not changing
      rev_N_old
    #from new_W
    rev_N_new <-
      sum(xor(up_down_vectorized_new[1:(length(up_down_vectorized_new) / 2)] ,
              up_down_vectorized_new[((length(up_down_vectorized_new) /
                                         2) + 1):length(up_down_vectorized_new)]))

    N_new <- sum(new_model$W == 0) +
      (sum(new_model$W != 0) - nstates) + #diag is not changing
      rev_N_new


    #ratio of Propose Probabilities or Hasting ratio
    #Q(M_new|M_old)/Q(Mold|M_new) =
    #( 1/(N(M_new)) / 1/(N(M_old)) )= N(M_old) /N(M_new)
    #M_new number of graphs that are neighbours of the
    #proposed structure

    # Metropolis-Hasting Ratio min{new_score/old_score,1}
    #(transition of old to new and new to old is symmetric?no)
    #A= (Likelihood(Data|M_new)*Prior(M_new)*(Q(M_old|M_New))/...
    #...Likelihood(Data|M_old)*Prior(M_old)*(Q(M_new|M_old)))
    #logscale
    # logratio= new_score - pre_score + Q
    new_score = prop_p + prop_ll + log(N_new)
    pre_score = pre_p + pre_ll + log(N_old)



    randratio = (runif(1, min = 0, max = 1))



    if ((new_score >= pre_score)) {
      #accept M_new

      acc_res = c(acc_res, list(new_model$W))
      acc_ratio = c(acc_ratio, prop_ll)

      all_res = c(all_res, list(new_model$W))
      all_ratio = c(all_ratio, prop_ll)

      model = new_model
      pre_ll = prop_ll
      pre_p = prop_p

    } else if (exp(new_score - pre_score) >= randratio) {
      #accept

      acc_res = c(acc_res, list(new_model$W))
      acc_ratio = c(acc_ratio, prop_ll)

      all_res = c(all_res, list(new_model$W))
      all_ratio = c(all_ratio, prop_ll)

      model = new_model
      pre_ll = prop_ll
      pre_p = prop_p
    } else{
      all_res = c(all_res, list(model$W))
      all_ratio = c(all_ratio, pre_ll)
    }

  } else{
    #invalid move_type on W reject W

    all_res = c(all_res, list(model$W))
    all_ratio = c(all_ratio, pre_ll)
  }
}



# returns model structures and LL
return(list(
  all_W = all_res,
  all_ll = all_ratio,
  acc_W = acc_res,
  acc_ll = acc_ratio
))
}
