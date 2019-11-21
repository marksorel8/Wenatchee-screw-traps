#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_MATRIX(log_R_obs);                  // Observed recruits (posterior from juvenile modle)
  DATA_IVECTOR(brood_year);
  DATA_IVECTOR(LH);
  DATA_VECTOR(wild_S_obs);                  // Observed spawners (from redd counts)
  DATA_VECTOR(hatch_S_obs);                  // Observed spawners (from redd counts)
  vector<Type> S_obs= wild_S_obs+hatch_S_obs;
  DATA_VECTOR(hatch_carcs);
  DATA_VECTOR(total_carcs);
  DATA_VECTOR(brood_rem);
  DATA_IVECTOR(Model);
  //---------------------------------------------------------------------
  PARAMETER(log_S_obs_cv);             //log spawner observation SD
  Type S_obs_cv = exp(log_S_obs_cv);        
  REPORT(S_obs_cv);
  PARAMETER_MATRIX(log_R_hat);    // latent true number of juveniles

  PARAMETER_MATRIX(log_R_obs_sd); // annual observation error of juveniles

  PARAMETER_VECTOR(log_W_ret_init);    // latent true number of spawners
  vector<Type> W_ret_init = exp(log_W_ret_init);
  REPORT(W_ret_init);

  
  PARAMETER_VECTOR(log_alpha);           // intrinsic productivity
  vector<Type> alpha =exp(log_alpha);
  REPORT(alpha);
  PARAMETER_VECTOR(log_R_max);           // Asymptotic maximum recruitment
  vector<Type> R_max =exp(log_R_max);
  REPORT(R_max);
  PARAMETER_VECTOR(log_d);
  vector<Type> d = exp(log_d);
  REPORT(d); 
  
  PARAMETER_VECTOR(log_proc_sigma);      // process error variances
  vector<Type> proc_sigma = exp(log_proc_sigma);
  REPORT(proc_sigma);
  
  PARAMETER_VECTOR(proc_er_corr);         // process error correlations
  REPORT(proc_er_corr)  
  
  PARAMETER_VECTOR(logit_pHOS);
  vector<Type> pHOS = invlogit(logit_pHOS);
  
  DATA_MATRIX(prop_age);
  
  PARAMETER_MATRIX(logit_surv);
  REPORT(logit_surv);
  //---------------------------------------------------------------------


 
 vector<Type>wild_return(wild_S_obs.size()); //total wild spawners in brood year +3
 wild_return.fill(0);

 vector<Type>S_hat(wild_S_obs.size());      //total spawners in brood year + 3
 vector<Type> hatch_S_hat(wild_S_obs.size());

 vector<Type> wild_S_hat(wild_S_obs.size());

 
 matrix<Type> ad_LH(wild_S_obs.size(),4); //wild adult returns by LH & BY
 
 matrix<Type> ad_LH_age(wild_S_obs.size(),prop_age.cols()); //wild adult returns by age, LH & BY

 matrix<Type> R_pred(wild_S_obs.size(),4); //predicted juvenile emigrants by BY
 
 
 //
 //process population model
 
 for (int Iyear=0; Iyear<wild_S_obs.size();Iyear++){
 
    
    Type wild_ret_yr = 0;
    for(int Ilh=0;Ilh<4; Ilh++){   
      ad_LH(Iyear,Ilh)= exp(log_R_hat(Iyear,Ilh))*
                            invlogit(logit_surv(Ilh,Iyear)); //wild adult returns by LH
    
    for (int Iage=0;Iage<3; Iage++ ){
      ad_LH_age(Iyear,(Ilh*3)+Iage) = prop_age(Iyear,(Ilh*3)+Iage)  * 
                                               ad_LH(Iyear,Ilh); //wild adult returns by LH & age
      if(Iyear>=5){
        wild_ret_yr+= ad_LH_age(Iyear-Iage-3,(Ilh*3)+Iage); //accumilate wild adult returns by LH and age
       } 
      }
     }
    wild_return(Iyear) = wild_ret_yr;
    
   if(Iyear<5){
     wild_return(Iyear) = W_ret_init(Iyear);
   }else{
     wild_return(Iyear) = wild_ret_yr;
   }
   
   wild_S_hat(Iyear)=   (wild_return(Iyear) - brood_rem(Iyear));
     
   hatch_S_hat(Iyear)= wild_S_hat(Iyear)*pHOS(Iyear)/
     (1-pHOS(Iyear));
  
  S_hat(Iyear) =  wild_S_hat(Iyear)+hatch_S_hat(Iyear); //estimate total spawners
   
   
 for(int Ilh=0;Ilh<4; Ilh++){  //predict juveniles from transiton models. 

 if(Model(Ilh)==1){
   R_pred(Iyear,Ilh)= R_max(Ilh)*(1-exp(-(pow(S_hat(Iyear)/alpha(Ilh),d(Ilh))))); //weibull CDF

 }

 if(Model(Ilh)==2){
   R_pred(Iyear,Ilh)= (alpha(Ilh)*S_hat(Iyear))/(1+alpha(Ilh)*S_hat(Iyear)/R_max(Ilh)); //beverton-holt model
    }
   }//LH loop
  }//year loop
 REPORT(wild_return);
 REPORT(hatch_S_hat);
 REPORT(wild_S_hat);
 REPORT(ad_LH_age);
 REPORT(ad_LH);
 
 //                End of population model
 //------------------------------------------------------------------------------------
 
 // Load namespace which contains the multivariate distributions
 using namespace density;
 
 //survival from emigration to adults
 
 PARAMETER(logit_Phi); //logit survival autocorrelation
 Type Phi = invlogit(logit_Phi);
 REPORT(Phi);
 
 PARAMETER_VECTOR(log_surv_var);
 vector<Type> surv_var = exp(log_surv_var);
 REPORT(surv_var);
 
 PARAMETER_VECTOR(surv_cor);
 REPORT(surv_cor);
 
 PARAMETER_VECTOR(surv_coefs);
 DATA_VECTOR(LH_DAYS);
 REPORT(surv_coefs);
 
 
 vector<Type>  mean_logit_surv_LH(4);
 REPORT(mean_logit_surv_LH);
 array<Type> logit_surv_inov(logit_surv.rows(),
                             logit_surv.cols());
 for(int Ilh=0;Ilh<4; Ilh++){
   mean_logit_surv_LH(Ilh)= surv_coefs(0)+ LH_DAYS(Ilh)*surv_coefs(1);
   
   for (int Iyear=0; Iyear<wild_S_obs.size();Iyear++){
   logit_surv_inov(Ilh,Iyear) = logit_surv(Ilh,Iyear)-
                                           mean_logit_surv_LH(Ilh);
   }
 }
 
 Type like_surv = AR1(Phi,VECSCALE(UNSTRUCTURED_CORR(surv_cor),surv_var))
   (logit_surv_inov); // nll 
 
 
 //-----------------------------------------------------
 
 
 // age at return
 
// PARAMETER_MATRIX(alr_p_hyper_mu);
// //PARAMETER_MATRIX(alr_p_mu_by);
// 
// matrix<Type> p_hyper_mu(2,3);
// p_hyper_mu.col(2)=Type(1)/
//                   (Type(1)+
//                     exp(vector<Type>(alr_p_hyper_mu.col(0)))+
//                     exp(vector<Type>(alr_p_hyper_mu.col(1))));
// 
// p_hyper_mu.col(1)=Type(1)/
//   (exp(vector<Type>(alr_p_hyper_mu.col(1)))+
//     exp(vector<Type>(alr_p_hyper_mu.col(0)))+
//     exp(vector<Type>(alr_p_hyper_mu.col(1))));
// 
// p_hyper_mu.col(0)=Type(1)/
//   (exp(vector<Type>(alr_p_hyper_mu.col(0)))+
//     exp(vector<Type>(alr_p_hyper_mu.col(0)))+
//     exp(vector<Type>(alr_p_hyper_mu.col(1))));
// 
// 
// REPORT(p_hyper_mu);
// 
// alr_p_age_mat=matrix(2,2);
// 
// for (int Iage=0;Iage<2;Iage++){
// alr_p_age_mat.col(Iage)=log(invlogit(vector<Type>(logit_p_age_mean.col(Iage)))/
//                             invlogit(vector<Type>(logit_p_age_mean.col(2))));
// }
// 
// 
// PARAMETER_VECTOR(log_ages_var);
// vector<Type> ages_var = exp(log_ages_var);
// 
// PARAMETER_VECTOR(age_cor);
// 
// 
// Type age_proc_like=0;
// 
// UNSTRUCTURED_CORR_t<Type> mvn(age_cor);
// 
// 
// alr_
// 
// for(int Ilh=0;Ilh<4; Ilh++){
//  for(int I=0;I<wild_S_obs.size(); I++ ){
//    age_proc_like+= VECSCALE(mvn,ages_var)
//   (vector<Type>( vector <Type> (prop_age.row(Iyear)).segment(Ilh*3,2))-
//     log(vector<Type>(R_pred.row(I)))));
//  }
// }
 
 //------------------------------------------
 //pHOS
 Type phos_like =0;
 for (int Iyear=0; Iyear<wild_S_obs.size();Iyear++){
  phos_like += dbinom_robust(hatch_carcs(Iyear),
                                total_carcs(Iyear),
                                logit_pHOS(Iyear),true);
 }
 REPORT(pHOS)
 //--------------------------------------------
 //juvenile "recruitment" MVN process error likelihood 
 Type state_Like=0;
 
   UNSTRUCTURED_CORR_t<Type> mvn(proc_er_corr);

    for(int I=0;I<log_R_hat.rows(); I++ ){
      state_Like+= VECSCALE(mvn,proc_sigma)
      (vector<Type>( vector <Type> (log_R_hat.row(I))- 
        log(vector<Type>(R_pred.row(I)))));
    }
  
 // matrix<Type> resids = log(R_hat)-log(R_pred);
  //REPORT(resids);
  
   Type Spawner_obs_like= sum(dnorm(log(S_obs),
                                     log(S_hat),
                                      S_obs_cv,
                                      true));    //observation likelihood (spawners)

  
  Type Recruit_obs_like= 0;                   //Initialize observatin likelihood (recruits)
  for (int I=0; I< log_R_obs.cols() ;I++){
    Recruit_obs_like+= dnorm((vector<Type>(log_R_obs.col(I))),
                                   log_R_hat(brood_year(I),LH(I)),
                                   exp(log_R_obs_sd(brood_year(I),LH(I))),true).sum(); //Observatin likelihood (recruits)
  }
  
  
  Type obj_fun = -state_Like-like_surv+phos_like+
                  Spawner_obs_like+ Recruit_obs_like;
  REPORT(R_pred);
  REPORT(state_Like);
  REPORT(Spawner_obs_like);
  REPORT(Recruit_obs_like);
  REPORT(obj_fun);
  REPORT(log_R_hat);
  REPORT(S_hat);
  return(-obj_fun);
}


