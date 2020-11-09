#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


//Remove the first maximum element from a vector

// [[Rcpp::export]]
NumericVector maxexceptfor(NumericVector vec, int ele) {
  int vs = vec.size();
  NumericVector res(vs-1);
  int j=0;
  for(int i=0; i<vs; i++) {
    if( i == ele ) {continue;}
    if(i!= ele) {res[j] = vec[i];
      j=j+1;}
  }
  return res;
}
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double rcppLikelihood(List model,
                      Function compute_states,
                      Function gammas ,
                      Function ddens1,
                      Function ddens0,
                      Function WT_R,
                      Function input,
                      String method,
                      double ts
                      ) {

  NumericMatrix m = model["Q"];
  NumericVector r0 = model["r0"];
  NumericMatrix theta =model["prior.theta"];
  List Para= model["data.para"];
  arma::mat W = model["W"];
  arma::cube dat = model["orig.dat"];
  
  std::string cx = Rcpp::as<std::string>(model["method.wt"]);
  CharacterVector methodWt= model["method.wt"];
  CharacterVector S = "steadyState";
  CharacterVector T = "timeSeries";

  int nstates = m.ncol();
  int times = dat.n_slices;


  //array to tube, tube to rowvec (send to ddens functions)
  arma::cube mem_ex;
  arma::rowvec cub_subset;
  //output of ddens
  NumericVector dat_vec1;
  NumericVector dat_vec0;

  //index holders
  NumericVector id(nstates);
  int id_max; //index of first max val


  double a= 0; //max holder
  double  ll=0 ; // LogLiklihood val


  //temporal vars
  NumericVector temp(nstates) ; // values from S loop
  double X; //float? check for 0 case (!!corrections)


  //R an gamma holders
  NumericVector gamma_vec;
  NumericVector R_vec;
  NumericMatrix gamma_mat;
  NumericMatrix R_mat;
  NumericVector R_vec_0;
  NumericMatrix R_mat_0;
  NumericVector Q_WT(nstates);

  Q_WT.fill(0);

  for(int i=0; i<nstates;i++){
    W.diag(0)[i]=model["depletion.factor.wt"];
  }
  if(methodWt[0][0]==S[0][0]){
    R_vec_0 = compute_states(W,
                             Q_WT,
                             r0,
                             (times),
                             model["method.wt"],
                             ts,
                             false,
                             model["stimulation.value.wt"]);

  }else if(methodWt[0][0]==T[0][0]){

     R_mat_0 = compute_states(W,
                              Q_WT,
                              r0,
                              (times),
                              model["method.wt"],
                              ts,
                              false,
                              model["stimulation.value.wt"]);
   }
  

  //******
  //iterate over all perturbation experiments
  for(int q = 0 ; q < (m.nrow()); q = q + 1){


    //calculate hidden states by solving ODE
    //calculate gammas
    if(times==1){

    } else{
      

      if(methodWt[0][0]==S[0][0]){
       R_mat = compute_states(model["W"],
                               m(q ,_),
                               R_vec_0,
                               (times),
                               method,
                               ts,
                               false,
                               model["stimulation.value"]);
       gamma_mat=gammas(R_mat,
                        R_vec_0,
                        false);

      }else if(methodWt[0][0]==T[0][0]){
        R_mat = compute_states(model["W"],
                               m(q ,_),
                               R_mat_0,
                               (times),
                               method,
                               ts,
                               false,
                               model["stimulation.value"]);
        gamma_mat=gammas(R_mat,
                         R_mat_0,
                         false);
      }
    }



    //compute Eq. (3), ensuring you never take log(0)
    //plug into formula in supplements:
    //iterate over all E-nodes
    //add the second formula in the Supplements to ll
    for(int i=0; i < (dat.n_rows); i++){

      //vectorize data of Obs_i under pert_q
      mem_ex=dat.tube(i,q);
      cub_subset= arma::rowvec(mem_ex.memptr(),mem_ex.n_elem,1,false );


      //f0 ,f1
      if(times==1){
     //   dat_vec1 = ddens1(as<double>(wrap(cub_subset)));
      //  dat_vec0 = ddens0(as<double>(wrap(cub_subset)));
      }else{
        dat_vec1= ddens1(as<NumericVector>(wrap(cub_subset)),
                         as<double>(Para["gumb.loc"]),
                         as<double>(Para["gumb.scl"]),
                         as<double>(Para["gumb.shp"]));
        dat_vec0= ddens0(as<NumericVector>(wrap(cub_subset)),
                         as<double>(Para["trunc.mean"]),
                         as<double>(Para["trunc.sd"]));
      }



      //iterate over all S-nodes
      //log(sum_S(exp(sum_T(log p(o_iq(t)| W;q; theta_is = 1) + logP(theta_is = 1)))))
      //a+ log (summ_S(exp(Temp-a)))--------->correction
      //a + log1p(sum(exp(Temp[-a]-a)))
      temp.fill(0);
      for(int s=0; s<nstates;s++)
      {

        //temp for time loop
        NumericVector temp_t(times);
        temp_t.fill(0);

        //temp_t=(log p(o_iq(t)| W;q; theta_is = 1) + logP(theta_is = 1))
        if(times==1){

          // temp_t[0]= log(
          //   ((gamma_vec(s)*dat_vec1[0]) +
          //   ((1-gamma_vec(s))*dat_vec0[0]))*
          //   theta(i,s)
          // );

        }else{
          for (int t=0; t<times;t++)
          {

            //value before log
            X=((gamma_mat(t,s)*dat_vec1[t]) +
              ((1-gamma_mat(t,s))*dat_vec0[t]))*
              theta(i,s);

            //check for 0 instead of adding epsilon!
            if(X!=0){
              temp_t[t]= log(X);
            }else{
              temp_t[t]= NA_REAL;
            }
          }}

        //xs in supplemet
        //sum over T(log p(o_iq(t)| W;q; theta_is = 1) + logP(theta_is = 1))
        // if all temp_t elements are real sum them for temp[s]
        // else temp[s] will be NA
        //it means for S all temp_t should be real val that temp become real!
        NumericVector v0 =  temp_t[!is_na(temp_t)]; //Rcout<<temp_t<<"\n";//Rcout<<v0<<"Len"<<v0.length()<<"\n";
        if(v0.length()==times){

          temp[s] = sum(temp_t);
        }else{

          temp[s] = NA_REAL;
        }

      }
      //Rcout<<temp<<"\n";

      //v1 =check which S has temp val != NA
      //if (length(v1!=NA)==0) return(-Inf)
      //if (length(v1!=NA)==1) return(that nonNA element)
      //else max(v1) + log1p(sum(exp(Temp[-maxv1]-maxv1)))
      NumericVector v1 = temp[!is_na(temp)];
      //Rcout<<v1<<"\n";
      if(v1.length()==0){
        return( -1 * arma::datum::inf);

      }else if(v1.length()==1){

        ll = ll+v1[0];
      }else{
        NumericVector temp_maxless(v1.length()-1);
        a=max(v1);
        id_max=which_max(v1);
        temp_maxless = maxexceptfor(v1,id_max);

        ll=ll+a+log1p(sum(exp(temp_maxless-a)));//(temp(-id_max))

      }




    }
  }



  //ret["ll"] = ll;
  //ret["Rvec"] = R_vec;
  return (ll);
}
