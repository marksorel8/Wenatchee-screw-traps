#include <TMB.hpp>

//  This model estimates production of juvenile emigrants exprexssing different life history strategies as a density -depenent function
//  of spawwner abundance with process error as a function of environmental covariates and synchronoys and idiosyncratic random effects. 
// 
// Copyright (C) 2021  Mark Sorel
// 
// This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU Affero General Public License as
//   published by the Free Software Foundation, either version 3 of the
//   License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU Affero General Public License for more details.
//   
//   You should have received a copy of the GNU Affero General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//     
//     I can be contacted at marks6@uw.edu or at:
//       Mark Sorel
//       1122 NE Boat Street,
//       Seattle, WA 98105
// 

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_sl );                // number of streams x life histories
  DATA_INTEGER( n_i );                 // number of years x streams x life histories
  DATA_INTEGER( n_t );                 // number of years 
  DATA_INTEGER( n_l );                 // number of life histories 
  DATA_INTEGER( n_s );                 // number of streams 
  DATA_INTEGER( n_f );                 // number of latent variable factors
  
  DATA_VECTOR(J_obs);                  // Observed log juveniles (posterior mean from juvenile model)
  DATA_VECTOR(J_obs_sd);               // Observed log juveniles standard error (posterior mean from juvenile modle)
  DATA_VECTOR(log_S_obs);              // log mean of observed spawners (log(redd count)-.5 CV^2)
  DATA_SCALAR(S_obs_CV);               // Observed spawners (from redd counts)
  
  DATA_SCALAR(log_J_max_prior);        // lognormal prior mean for Jmax hyper-means
  
  DATA_FACTOR(st_i);                   // stream by year index for observation i
  DATA_FACTOR(sl_i);                   // stream and life history index for observation i
  DATA_FACTOR(t_i);                    // year index for observation i
  DATA_SPARSE_MATRIX(X);               // Design matrix of covariates for process error model
  //---------------------------------------------------------------------------------------------------
  
  //Parameters
  PARAMETER_VECTOR(beta_alpha);        // log alpha hyper means (LHP-specific)
  REPORT(beta_alpha); 
  PARAMETER_VECTOR(beta_Jmax);         // log Jmax hyper means (LHP-specific)
  REPORT(beta_Jmax); 
  PARAMETER_VECTOR(beta_gamma);        // log gamma hyper means (LHP-specific)
  REPORT(beta_gamma); 
  
  PARAMETER_VECTOR(log_sigma_alpha);   // log standard deviations of stream-specific log alpha random effects (LHP-specific)
  REPORT(log_sigma_alpha);
  PARAMETER_VECTOR(log_sigma_Jmax);    // log standard deviations of stream-specific log Jmax random effects (LHP-specific)
  REPORT(log_sigma_Jmax);
  PARAMETER_VECTOR(log_sigma_gamma);   // log standard deviations of stream-specific log gamma random effects (LHP-specific)
  REPORT(log_sigma_gamma);
  
  PARAMETER_VECTOR(beta_e);            // environmental variable coefficients           
  PARAMETER_VECTOR(Loadings_vec);      // latent variable factor loadings
  PARAMETER_VECTOR(log_sigma_eta);     // log idiosyncratic process errors standard deviations
  vector<Type> sigma_eta = exp(log_sigma_eta);
  REPORT(sigma_eta);
  
  PARAMETER_VECTOR(log_sigma_Bgamma);  // standard deviation of regularizing prior on hyper-means of log gamma
  PARAMETER_VECTOR(log_sigma_BJmax);   // standard deviation of regularizing prior on hyper-means of log Jmax


  // Random effects
  PARAMETER_VECTOR(eps_alpha);         // random effect for stream- and LHP-specific alpa
  REPORT(eps_alpha);
  PARAMETER_VECTOR(eps_Jmax);          // random effect for stream- and LHP-specific Jmax
  REPORT(eps_Jmax);
  PARAMETER_VECTOR(eps_gamma);         // random effect for stream- and LHP-specific gamma
  REPORT(eps_gamma);
  PARAMETER_VECTOR(log_S_hat);         // log latent female spawner abundance
  vector<Type> S_hat =exp(log_S_hat);
  REPORT(S_hat);
  PARAMETER_MATRIX(Omega_xf);          // latent factor variables
  REPORT(Omega_xf);
  PARAMETER_VECTOR(eta);               // latent stream by year by life history (idiosyncratic) process errors in juvenile recruitment
  //---------------------------------------------------------------------------------------------------  
  
  //Variables
  
  //// Initialize likelihood
  vector<Type> jnll_comp(3);   
  jnll_comp.setZero();   
  
  //// Unpack factor loadings matrix  (taken from J. Thorson spatial class example)
  matrix<Type> Loadings_pf((n_sl), n_f);
  int Count = 0;
  for(int fac=0; fac<n_f; fac++){
    for(int p=0; p<(n_sl); p++){
      if(p==fac){
        Loadings_pf(p,fac) = exp(Loadings_vec(Count));
        Count++;
      }else{
      if(p>=fac){
        Loadings_pf(p,fac) = Loadings_vec(Count);
      
        Count++;
      }else{
        Loadings_pf(p,fac) = 0.0;
      }
    }
  }
  }
  REPORT(Loadings_pf);
    
  

  Type alpha_bias_correction = 0;// bias correction for lognormal mean  
  Type gamma_bias_correction = 0;
  Type Jmax_bias_correction = 0;

  //// calculate alphas, gammas, and Jmaxes
  //////linear predictors on log scale
  for(int j=0; j<n_l; j++){          // loop over life histories
    alpha_bias_correction = pow(exp(log_sigma_alpha(j)),2)/2.0;
    Jmax_bias_correction = pow(exp(log_sigma_Jmax(j)),2)/2.0;
    gamma_bias_correction = pow(exp(log_sigma_gamma(j)),2)/2.0;
    for( int i=0; i<n_s; i++){       // loop over streams
      
      jnll_comp(0) -= dnorm( Type(eps_alpha(i*n_l+j)),
                Type(beta_alpha(j) -alpha_bias_correction),
                Type(exp(log_sigma_alpha(j))),true );
      
      jnll_comp(0) -= dnorm( Type(eps_Jmax(i*n_l+j)),
                Type(beta_Jmax(j) -Jmax_bias_correction),
                Type(exp(log_sigma_Jmax(j))),true );
      
      jnll_comp(0) -= dnorm( Type(eps_gamma(i*n_l+j)),
                Type(beta_gamma(j) -gamma_bias_correction),
                Type(exp(log_sigma_gamma(j))),true );
      
    }
  }
  
  
  
  ////// transform to positive real space
  vector<Type> alpha = exp(eps_alpha);
  vector<Type> Jmax = exp(eps_Jmax);
  vector<Type> gamma =  exp(eps_gamma);

  
  // calculate latent juvenile emigrants
  vector<Type> J_hat(n_i);      // vector to hold latent juveniles expectations (without process error)
  vector<Type> J_pred(n_i);     // vector to hold latent juveniles (with process error)
  
  //// covariate effects on process errors
  vector<Type> cov_e = X * beta_e; // design matrix * coefficients
  
  // latent variable factor effects on process errors
  matrix<Type> LV_effects = Loadings_pf* Omega_xf.transpose();
  REPORT(LV_effects);
  
  for( int i=0; i<n_i; i++){      // loop over all observations of streams, life-histories, and years
   
   ///// calculate expected juveniles based on spawners and functional form (without process error)
      //depensatory Beverton holt II (Myers et al 1995)
      J_hat(i)=  (alpha(sl_i(i))* pow(S_hat(st_i(i)),gamma(sl_i(i)))) /
        (Type(1.0)+   alpha(sl_i(i))* pow(S_hat(st_i(i)),gamma(sl_i(i))) / Jmax(sl_i(i)));

    
   
        J_pred(i) = J_hat(i) * exp(cov_e(i)+ LV_effects(sl_i(i), t_i(i))+
    eta(i)*sigma_eta(sl_i(i)));

  }
  REPORT(J_pred);
//---------------------------------------------------------------------------------------------------  

  // Probabilities of random effects
  // 
  PARAMETER_VECTOR(rate);
    

  //regulatizing priors on hyper-standard deviations of random effects on log alphas, log gammas, and log Jmaxes
  jnll_comp(0) -= (dexp(exp(log_sigma_alpha),Type(exp(rate(0))),true).sum() + 
    log_sigma_alpha.sum());
  
  jnll_comp(0) -= (dexp(exp(log_sigma_gamma),Type(exp(rate(0))),true).sum() + 
    log_sigma_gamma.sum());
  
  jnll_comp(0) -= (dexp(exp(log_sigma_Jmax),Type(exp(rate(0))),true).sum() +
    log_sigma_Jmax.sum());

  

    //reluarizing prior on hyper-means for log gammas
    jnll_comp(0) -= dnorm(beta_gamma,Type(0),exp(log_sigma_Bgamma),true).sum();
    jnll_comp(0) -= (dexp(exp(log_sigma_Bgamma),Type(exp(rate(1))),true).sum() +
      log_sigma_Bgamma.sum());

    //regularizing prior on hyper-mean of log-Jmaxes
    jnll_comp(0) -= dnorm(beta_Jmax,log_J_max_prior,exp(log_sigma_BJmax),true).sum();
    jnll_comp(0) -= (dexp(exp(log_sigma_BJmax),Type(exp(rate(2))),true).sum() +
      log_sigma_BJmax.sum());
    
    //regularizing prior on idiosyncratic process error standard deviations
    jnll_comp(0) -= (dexp(exp(log_sigma_eta),Type(exp(rate(3))),true).sum() +
    log_sigma_eta.sum());

   
  
  //// latent factor variables
  for(int t = 0; t<(n_t); t++){ // loop over years
    jnll_comp(1)-= dnorm(vector<Type>(Omega_xf.row(t)),Type(0),Type(1),true).sum();
  }

  // idiosyncratic process error
  jnll_comp(1)-=dnorm(eta ,
            Type(0.0) ,
            //            sigma_eta,
            Type(1.0),
            true ).sum();

  
  //Likelihood
  
  ////  Latent juveniles
  vector<Type> juv_like= dnorm(J_obs ,vector<Type>(log(J_pred)) ,J_obs_sd , true );
  REPORT(juv_like);
  jnll_comp(2) -= juv_like.sum();  
  
  ////  Latent spawners
  jnll_comp(2) -= dnorm(log_S_obs ,log_S_hat  ,S_obs_CV , true ).sum();
  
  
  

//Reporting
  REPORT(gamma);
  REPORT(Jmax);
  REPORT(alpha);
  REPORT(beta_e); 
  REPORT(jnll_comp);
  ADREPORT(log(J_pred));   //get standard errors for log unobserved juveniles
  ADREPORT(log_S_hat);     //get standard errors for log unobserved spaweners
  Type obj_fun = jnll_comp.sum();
  return(obj_fun); 
}
