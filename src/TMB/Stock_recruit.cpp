#include <TMB.hpp>

// dlognorm (copied from J. Thorson DFA model)
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=false){
  Type Return;
  if(give_log==false) Return = dnorm( log(x), meanlog, sdlog, false) / x;
  if(give_log==true) Return = dnorm( log(x), meanlog, sdlog, true) - log(x);
  return Return;
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_sl );                // number of streams x life histories
  DATA_INTEGER( n_i );                 // number of years x streams x life histories
  DATA_INTEGER( n_t );                 // number of years 
  DATA_INTEGER( n_st );                 // number of streams X years
  DATA_INTEGER( n_l );                 // number of life histories 
  DATA_INTEGER( n_s);                 // number of streams 
  DATA_INTEGER( n_beta_e);            // number of environemntal covariates 
  
  DATA_IVECTOR( mod );                 // process model (e.g. BH) 
  
  DATA_VECTOR(R_obs);                  // Observed log recruits (posterior mean from juvenile modle)
  DATA_VECTOR(R_obs_sd);               // Observed log recruits standard error (posterior mean from juvenile modle)
  DATA_VECTOR(log_S_obs);              // Observed spawners (from redd counts)
  DATA_SCALAR(S_obs_CV);               // Observed spawners (from redd counts)
  DATA_FACTOR(l_i);                    // life history index for observation i
  DATA_FACTOR(st_i);                   // stream by year index for observation i
  DATA_FACTOR(lt_i);                   // life history by year index for observation i
  DATA_FACTOR(sl_i);                   // stream and life history index for observation i
  DATA_FACTOR(t_i);                    // year index for observation i
  DATA_FACTOR(s_i);                    // year index for observation i
  DATA_FACTOR(s_theta_i);              // stream index for theta i
  DATA_FACTOR(beta_e_i);               // stream by covariate index for rand covariate coefficients
  DATA_SPARSE_MATRIX(X);               // Design matrix of covariates for process error model
  DATA_SCALAR(num);
  
  //Parameters
  PARAMETER_VECTOR(beta_a);           // intrinsic productivity intercepts
  REPORT(beta_a); 
  PARAMETER_VECTOR(beta_rm);          // intrinsic productivity intercepts and coefficients
  REPORT(beta_rm); 
  PARAMETER_VECTOR(beta_e);    // covariate coefficients for process error model 
  REPORT(beta_e);
  //PARAMETER_VECTOR(hyper_beta_e); //hyper-mean of covariate effects across streams
  // REPORT(hyper_beta_e);
  // ADREPORT(hyper_beta_e);
  // // PARAMETER_VECTOR(log_sigma_e);  //sd of covariate effects accross streams
  // vector<Type> sigma_e = exp(log_sigma_e);
  // ADREPORT(log_sigma_e);
  
  //vector<Type> sigma_e = exp(log_sigma_e);
  // PARAMETER_VECTOR(zeta_a);          // 
  // REPORT(zeta_a);
  // 
  // PARAMETER_VECTOR(eta_a);          // 
  // REPORT(eta_a);
  //
  // PARAMETER_VECTOR(zeta_rm);          // 
  // REPORT(zeta_rm);
  // 
  // PARAMETER_VECTOR(eta_rm);          // 
  // REPORT(eta_rm);
  // 
  // PARAMETER_VECTOR(zeta_f);          // 
  // REPORT(zeta_f);
  // 
  // PARAMETER_VECTOR(eta_f);          // 
  // REPORT(eta_f);
  // 
  
   // PARAMETER_VECTOR(log_sigma_a);      // productivity random effect standard deviations
   // vector<Type> sigma_a = exp(log_sigma_a);
   // REPORT(sigma_a);
   // PARAMETER_VECTOR(log_sigma_rm);     // maximum recruitment random effect standard deviations
   // vector<Type> sigma_rm = exp(log_sigma_rm);
   // REPORT(sigma_rm);
   // PARAMETER_VECTOR(log_sigma_f);      // productivity random effect standard deviations
   // vector<Type> sigma_f = exp(log_sigma_f);
   
   // 
  // PARAMETER(log_sigma_zeta_a);     // 
  // Type sigma_zeta_a = exp(log_sigma_zeta_a);
  // REPORT(sigma_zeta_a);
  // 
  // 
  // PARAMETER(log_sigma_eta_a);     // 
  // Type sigma_eta_a = exp(log_sigma_eta_a);
  // REPORT(sigma_eta_a);
  
  
  // PARAMETER(log_sigma_zeta_rm);     // 
  // Type sigma_zeta_rm = exp(log_sigma_zeta_rm);
  // REPORT(sigma_zeta_rm);
  // 
  // 
  // PARAMETER(log_sigma_eta_rm);     // m
  // Type sigma_eta_rm = exp(log_sigma_eta_rm);
  // REPORT(sigma_eta_rm);
  
  // 
  // 
  // PARAMETER(log_sigma_zeta_f);     // 
  // Type sigma_zeta_f = exp(log_sigma_zeta_f);
  // REPORT(sigma_zeta_f);
  // 
  // 
  // PARAMETER(log_sigma_eta_f);     // m
  // Type sigma_eta_f = exp(log_sigma_eta_f);
  // REPORT(sigma_eta_f);
  // 
  

  PARAMETER_VECTOR(beta_f);          // intrinsic productivity intercepts and coefficients
  ADREPORT(beta_f);
  // PARAMETER(log_sigma_f);      // productivity random effect standard deviations
  // Type sigma_f = exp(log_sigma_f);
  // 
  //PARAMETER_VECTOR(par_theta);//
  
  // PARAMETER(log_proc_sigma);    // annual process errors standard deviation
  // Type proc_sigma = exp(log_proc_sigma);
  // REPORT(proc_sigma);
  // 

  PARAMETER_VECTOR(log_proc_sigma);    // year by stream by life history process errors standard deviations
  vector<Type> proc_sigma = exp(log_proc_sigma);
  REPORT(proc_sigma);

  // PARAMETER(log_sigma_delta);               // year sds
  // Type sigma_delta = exp(log_sigma_delta);
  // REPORT(sigma_delta);
  // 
  PARAMETER_VECTOR(log_sigma_theta);               // streamn by year sds
  vector<Type> sigma_theta = exp(log_sigma_theta);
  REPORT(sigma_theta);
  
  // PARAMETER_VECTOR(log_sigma_e);               // stream by year sds
  // vector<Type> sigma_e = (log_sigma_e);
  // REPORT(sigma_e);
  // 
  
 

  // PARAMETER_VECTOR(log_sigma_log_fifty);
  // vector<Type> sigma_log_fifty = exp(log_sigma_log_fifty);
  // 
  // PARAMETER_VECTOR(log_sigma_log_rm);
  // vector<Type> sigma_log_rm = exp(log_sigma_log_rm);
  // 
  // // PARAMETER_VECTOR(log_sigma_kappa);               // life history by year sds
  //vector<Type> sigma_kappa = exp(log_sigma_kappa);
  //REPORT(sigma_kappa);
  // PARAMETER_VECTOR(proc_theta);    // annual process errors standard deviation
  // //vector<Type> proc_theta = invlogit(logit_proc_theta);
  // REPORT(proc_theta);
  
  
  
  // Random effects
  PARAMETER_VECTOR(log_S_hat);
  vector<Type> S_hat =exp(log_S_hat);
  REPORT(S_hat);
  PARAMETER_VECTOR(eps_a);           // random effect for stream- and life-stage-specific productivity
  REPORT(eps_a);
  PARAMETER_VECTOR(eps_rm);          // random effect for stream- and life-stage-specific maximum recruitment
  REPORT(eps_rm);
  PARAMETER_VECTOR(eps_f);           // random effect for stream- and life-stage-specific maximum recruitment
  REPORT(eps_f);

  
//PARAMETER_MATRIX(gamma);           
PARAMETER_VECTOR(gamma);           // stream by year by life history process errors in juvenile recruitment

 REPORT(gamma);
  matrix<Type>gamma_sim(gamma.size(),int(100));
 SIMULATE{
   for ( int i=0 ;i<gamma.size();i++){
     for(int j=0 ;j<int(100);j++){
     gamma_sim(i,j)=rnorm(Type(0),Type(1));
   }
   }
 }
  
  //PARAMETER_VECTOR(delta);           // year process errors in juvenile recruitment
  
 // REPORT(delta);
  
 PARAMETER_VECTOR(theta);           // stream by year process errors in juvenile recruitment
  
  REPORT(theta);
    matrix<Type>theta_sim(theta.size(),int(100));  
  SIMULATE{

    for ( int i=0 ;i<theta.size();i++){
      for(int j=0 ;j<int(100);j++){
        theta_sim(i,j)=rnorm(Type(0),Type(1));
      }
    }
    
  }
  
  //PARAMETER_VECTOR(kappa);            // life history by year  process errors in juvenile recruitment
  
  //REPORT(kappa);
  
  
  //objective function 
  using namespace density;
  vector<Type> jnll_comp(3);   //initialize likelihood components
  jnll_comp.setZero();   
  
  
  //derived quantitites
  vector<Type> log_alpha(n_sl); // vector to hold intrinsic productivites
  vector<Type> log_R_max(n_sl); // vector to hold asymptotic maximum recruitments
  vector<Type> log_fifty(n_sl);
  
  

  DATA_IVECTOR(pen_flag);
  DATA_VECTOR(prior_mean);
  DATA_VECTOR(exp_rates);
  //

  //calculate intrinsic productivity and asymptotic maximum recruitment
    for(int j=0; j<n_l; j++){    // loop over life histories
      for( int i=0; i<n_s; i++){     // loop over streams


      log_alpha(i*n_l+j) =          // log intrinsic productivity
        beta_a(j) +               // life-history intercept
        // zeta_a(i)+            //
        // eta_a(j)+             //
        eps_a(i*n_l+j);             // stream*LH random effect


      log_R_max(i*n_l+j) =          // log maximum recruitment
        beta_rm(j) +              // life-history intercept

        //zeta_rm(i)+            //
        //eta_rm(j)+             //
        //X_rm(i) +    // life-history effect of habitat area
        eps_rm(i*n_l+j);            // stream random effect ~ N(0,sigma_rm[life history])


      log_fifty(i*n_l+j) =         // depensation (change to proportion of R_max?)
        beta_f(j) +
       // zeta_f(i)+            //
       // eta_f(j)+             //
        eps_f(i*n_l+j);


      //probabilities of random effects in life-history and stream- specific
      // productivity and maximum recruitment
    //
     //UNSTRUCTURED_CORR_t<Type> mvnorm_proc(par_theta);
     // matrix<Type> test2 = mvnorm_proc.cov();
     // REPORT(test2);
     //  vector<Type> par_sigma(3);
     //  par_sigma(0)= sigma_a;
     //  par_sigma(1)= sigma_rm;
     //  par_sigma(2)= sigma_f;
    //
      // vector<Type> par(3);
      // par(0)= eps_a(i*4+j);
      // par(1)= eps_rm(i*4+j);
      // par(2)= eps_f(i*4+j);
    //
    //
     //jnll_comp(0)+=VECSCALE(mvnorm_proc,par_sigma)(par);
    //jnll_comp(0)-=dnorm(zeta_a(i), Type(0.0), sigma_zeta_a, true);
    //jnll_comp(0)-=dnorm(zeta_rm(i), Type(0.0), sigma_zeta_rm,true);
    //jnll_comp(0)-=dnorm(eta_a(j), Type(0.0), sigma_eta_a, true);
    //jnll_comp(0)-=dnorm(eta_rm(j), Type(0.0), sigma_eta_rm,true);
    //jnll_comp(0)-=dnorm(zeta_f(i), Type(0.0), sigma_zeta_f,true);
    //jnll_comp(0)-=dnorm(eta_f(j), Type(0.0), sigma_eta_f, true);
      //
      // jnll_comp(0)-=dnorm(eps_a(i*n_l+j), Type(0.0), sigma_a(j), true);
      // jnll_comp(0)-=dnorm(eps_rm(i*n_l+j), Type(0.0), sigma_rm(j),true);
      
       

       // jnll_comp(0)-= (dexp(sigma_log_fifty(i*n_l+j),Type(1),true)
       //                       +log_sigma_log_fifty(i*n_l+j));
// DATA_SCALAR(neg_cor);
//        matrix<Type> Sigma(2,2);
//        Sigma.fill(neg_cor);
//        Sigma(0,0)=sigma_log_fifty(i*n_l+j);
//          Sigma(1,1)=sigma_log_rm(i*n_l+j);
// 
// MVNORM_t<Type> N_0_Sigma(Sigma);
// REPORT(N_0_Sigma.cov());
// 
// DATA_SCALAR(rmax_prior);

// vector<Type>tmp(2);
// tmp(0)=eps_f(i*n_l+j)-log(f_prior);
// tmp(1)=eps_rm(i*n_l+j)-log(Type(rmax_prior));

// jnll_comp(0)+=N_0_Sigma(tmp);
// // 

// 
// if(pen_flag(0)){
// 
// jnll_comp(0)-=dnorm(eps_f(i*n_l+j),Type(log(prior_mean(0))),sigma_log_fifty(i*n_l+j), true);  //normal (centered on zero) priors on coefficients
// 
// 
// DATA_VECTOR(exp_rates);
//  jnll_comp(0)-= (dexp(sigma_log_fifty(i*n_l+j),Type(exp_rates(0)),true)
//                    +log_sigma_log_fifty(i*n_l+j));
// }
// 
// if(pen_flag(1)){
//   
//   jnll_comp(0)-=dnorm(eps_rm(i*n_l+j),Type(log(prior_mean(1))),sigma_log_rm(i*n_l+j), true);  //normal (centered on zero) priors on coefficients
//   
//   
// 
//   jnll_comp(0)-= (dexp(sigma_log_rm(i*n_l+j),Type(exp_rates(1)),true)
//                     +log_sigma_log_rm(i*n_l+j));
// }


//
//
//  jnll_comp(0)-=dnorm(eps_rm(i*n_l+j),log(Type(4e6)),sigma_log_rm(i*n_l+j), true);  //normal (centered on zero) priors on coefficients
//
//  jnll_comp(0)-= (dexp(sigma_log_rm(i*n_l+j),Type(exp_rates(1)),true)
// +log_sigma_log_rm(i*n_l+j));

// jnll_comp(0)-= (dexp(sigma_log_fifty(i*n_l+j),Type(1.0),true)
//                  +log_sigma_log_fifty(i*n_l+j));



    }
      //jnll_comp(0)-= (dexp(sigma_log_fifty(j),Type(5),true)
      //                  +log_sigma_log_fifty(j));
  }


  // for(int j=0; j<n_l; j++){    // loop over life histories
  //

  //
  // }
  //
  //  REPORT(sigma_log_fifty);
//pen alized complexity
    // jnll_comp(0)+=sum(dexp(sigma_a,Type(1),true))
    // -sum(sigma_a);
    // jnll_comp(0)+=sum(dexp(sigma_rm,Type(1),true))
    // -sum(sigma_rm);




  //transform productivities and maximum recruitment to positive real space
  vector<Type> alpha = exp(log_alpha);
  REPORT(alpha);
  vector<Type> R_max = exp(log_R_max);
  REPORT(R_max);
  vector<Type> fifty =  exp(log_fifty);
   REPORT(fifty);

   
   





  //derived quantitites
  vector<Type> R_hat(n_i);     // vector to hold latent recruits
  vector<Type> R_pred(n_i);     // vector to hold latent recruits

  //vector to hold unscaled random coefficients
  // vector<Type> beta_e_2(n_beta_e);//=beta_e;
  // //random covariate effect probabililties 
  // for( int i=0; i<n_beta_e; i++){
  //   jnll_comp(1)-=dnorm(beta_e(i) ,
  //             Type(0.0) ,
  //             //            proc_sigma,
  //             Type(1.0),
  //             true );
  //   beta_e_2(i)=beta_e(i)* sigma_e(beta_e_i(i))+hyper_beta_e(beta_e_i(i));//
  // }
  // ADREPORT(beta_e_2);
  // REPORT(beta_e_2);
  // covariate effects on errors
  vector<Type> cov_e = X * beta_e; // design matrix * coefficients
  
  //penalty on coefficients
  //jnll_comp(0)-=dnorm(beta_e,Type(0),Type(1),true).sum();
  //jnll_comp(0)-=(dexp(sigma_e,Type(1),true).sum()+log_sigma_e.sum());

Type d = 10.0;
for( int i=0; i<n_i; i++){      // loop over all streams, life-histories, and years
  //latent recruits based on Beverton-Holt model with lognormal errors

 if(mod(sl_i(i))==1){
//depensatory Beverton holt I
R_hat(i)=  (alpha(sl_i(i))* pow(S_hat(st_i(i)),fifty(sl_i(i)))) /
  (Type(1.0)+   alpha(sl_i(i))* pow(S_hat(st_i(i)),fifty(sl_i(i))) / R_max(sl_i(i)));
 }
   else { if(mod(sl_i(i))==2){
// depensatory Beverton holt II
R_hat(i)= (S_hat(st_i(i))/(S_hat(st_i(i))+fifty(sl_i(i))))  *   (alpha(sl_i(i))*S_hat(st_i(i))) /
(Type(1.0)/alpha(sl_i(i))+   S_hat(st_i(i)) / R_max(sl_i(i)));


// R_pred(i)= (Type(1)-exp((log(Type(0.5))/fifty(sl_i(i)))* S_hat(st_i(i))))  *   (alpha(sl_i(i))*S_hat(st_i(i))) /
// (Type(1.0)/alpha(sl_i(i))+   S_hat(st_i(i)) / R_max(sl_i(i)));
     
 } else{if(mod(sl_i(i))==3){
// depensatory smooth hockey stock
R_hat(i)= (S_hat(st_i(i))/(S_hat(st_i(i))+fifty(sl_i(i))))  * alpha(sl_i(i))*d*R_max(sl_i(i))*(1.0+exp(-1.0/d))*(S_hat(st_i(i))/(d*R_max(sl_i(i))) -
log((1+exp((S_hat(st_i(i))-R_max(sl_i(i)))/(d*R_max(sl_i(i)))))/(1.0+exp(-1.0/d))));
 } else{if(mod(sl_i(i))==4){
 //depensatory Ricker
 R_hat(i)= (S_hat(st_i(i))/(S_hat(st_i(i))+fifty(sl_i(i)))) * alpha(sl_i(i))* S_hat(st_i(i)) *exp(-1*R_max(sl_i(i))*S_hat(st_i(i)));

 } else{if(mod(sl_i(i))==5){
   //depensatory linear
   R_hat(i)=  (S_hat(st_i(i))/(S_hat(st_i(i))+fifty(sl_i(i)))) * (alpha(sl_i(i))* S_hat(st_i(i)));

 } else{if(mod(sl_i(i))==6){
   //Beverton Holt
   R_hat(i)=  (alpha(sl_i(i))* S_hat(st_i(i))) /
     (Type(1.0)+   alpha(sl_i(i))*S_hat(st_i(i)) / R_max(sl_i(i)));

 }else{if(mod(sl_i(i))==7){
   // smooth hockey stock
   R_hat(i)=  alpha(sl_i(i))*d*R_max(sl_i(i))*(1.0+exp(-1.0/d))*(S_hat(st_i(i))/(d*R_max(sl_i(i))) -
     log((1+exp((S_hat(st_i(i))-R_max(sl_i(i)))/(d*R_max(sl_i(i)))))/(1.0+exp(-1.0/d))));

 } else{if(mod(sl_i(i))==8){
   //Ricker
   R_hat(i)=  alpha(sl_i(i))* S_hat(st_i(i)) *exp(-1*R_max(sl_i(i))*S_hat(st_i(i)));

 } else{if(mod(sl_i(i))==9){
   //linear
   R_hat(i)=  (alpha(sl_i(i))* S_hat(st_i(i)));

 } else{if(mod(sl_i(i))==10){
   //Power function
   R_hat(i)=  alpha(sl_i(i))* pow(S_hat(st_i(i)),fifty(sl_i(i)));
 } else{if(mod(sl_i(i))==11){
   //Power function
   R_hat(i)=  (R_max(sl_i(i))*pow(S_hat(st_i(i)),fifty(sl_i(i))))/( (1/alpha(sl_i(i)))+
     pow(S_hat(st_i(i)),fifty(sl_i(i)))); //Type III functional response

 }

 }
 }
 }
 }
 }
 }
 }
 }
 }
}

    matrix<Type> sim_out(int(100),n_i);
   SIMULATE{
     for(int sim=0; sim<int(100); sim++){
      sim_out(sim,i)= R_hat(i) * exp(gamma_sim(i,sim)*proc_sigma(sl_i(i))+
        theta_sim(st_i(i),sim)*sigma_theta(s_i(i)));
     }
     REPORT(sim_out);
   }

//multiplicitive process errors
 R_pred(i) = R_hat(i)* exp(cov_e(i)+gamma(i)*proc_sigma(sl_i(i))+theta(st_i(i))*sigma_theta(s_i(i)));//+kappa(lt_i(i))))+delta(t_i(i));


  jnll_comp(1)-=dnorm(gamma(i) ,
            Type(0.0) ,
   //            proc_sigma,
   Type(1.0),
            true );



  //probability of latent recruits conditional on random effects
  //the probability is based on the estimation model for the number of juvenile emigrants
  jnll_comp(2) -= dnorm(R_obs(i) ,log(R_pred(i)) ,R_obs_sd(i) , true );
  //probability of eggs conditional on random effects

}

  //the probability is based on the estimation model for the number of juvenile emigrants
  jnll_comp(2) -= dnorm(log_S_obs ,log_S_hat  ,S_obs_CV , true ).sum();
//year random effect probs
 //  for( int i=0; i<n_t; i++){
 // jnll_comp(1)-= dnorm(delta(i), Type(0.0), sigma_delta, true );
 //
 //  }
//
// //life history by year random effect probs
//   for( int i=0; i<(n_t-1); i++){
//   for (int j = 0; j<n_l; j++){
//     jnll_comp(1)-= dnorm(kappa(i*n_l+j), Type(0.0), sigma_kappa(j), true );
//   }
//   }
//
//stream by year rand effect probs
 
    jnll_comp(1)-= dnorm(theta, Type(0.0), Type(1.0), true ).sum();


  //probability of random effects (process errors in recruitment)

  // UNSTRUCTURED_CORR_t<Type> mvnorn_proc(proc_theta);
  // REPORT(mvnorn_proc.cov());
  //   for( int i=0; i<n_t; i++){
  // jnll_comp(1)+=VECSCALE(mvnorn_proc,proc_sigma)(vector<Type>(gamma.row(i)));
  //   }

  vector<Type> ICC_sl(n_sl);
  vector<Type> ICC_s(n_s);
  for( int i=0; i<n_s; i++){
    for( int j=0; j<n_l; j++){
      ICC_sl(i*n_l+j)=pow(sigma_theta(i),2)/(pow(sigma_theta(i),2)+pow(proc_sigma(i*n_l+j),2));
    }
    ICC_s(i)=ICC_sl.segment(i*n_l,n_l).mean();
    }
ADREPORT(logit(ICC_sl));
ADREPORT(logit(ICC_s));
ADREPORT(ICC_sl);
ADREPORT(ICC_s);  
  
//return objective function 
REPORT(R_pred);
REPORT(jnll_comp);
Type obj_fun = jnll_comp.sum();
ADREPORT(log(R_pred));   //get standard deviations for unobserved emigrants
ADREPORT(log_S_hat);
ADREPORT(sigma_theta);
ADREPORT(proc_sigma);
return(obj_fun); 
}
