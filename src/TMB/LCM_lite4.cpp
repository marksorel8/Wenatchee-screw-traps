#include <TMB.hpp>


//function to backtransform additive log ratios onto the simplex
//and one for th yearling age proportions
template<class Type> matrix<Type> alr_to_simplex(matrix<Type> XX){
  int nCols = XX.cols()+2;
  int nRows = XX.rows();
  
  matrix<Type> simplex(nRows,nCols);
  
  for (int I =0; I<nCols/2;I++ ){
    for(int J =0; J<nRows;J++ ){
      for ( int K = 0; K<2; K ++){
        int L=3*K;
        int M=L+I;
      Type numerator = Type(1);
        if(I<nCols/2-1) numerator = exp(XX(J,2*K+I));
      simplex(J,M)=numerator/(Type(1)+
                             exp(vector<Type>(XX.row(J)).
                             segment((nCols/2-1)*K,nCols/2-1)).sum());
    }
   }
  }
  
return(simplex);
  }


template<class Type>
Type objective_function<Type>::operator() ()
{

  //Data  
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
  DATA_IVECTOR(JLH_A);
  DATA_IVECTOR(PIT_age_BYs);
  DATA_IVECTOR(PIT_age_LHs);
  DATA_MATRIX(PIT_ages);
  DATA_MATRIX(Carc_adult_age);

  
  //---------------------------------------------------------------------
  // Parameters
  PARAMETER(log_S_obs_cv);             //log spawner observation SD
  Type S_obs_cv = exp(log_S_obs_cv);        
  REPORT(S_obs_cv);
  
  PARAMETER_MATRIX(log_R_hat);         // latent true number of juveniles
  PARAMETER_MATRIX(log_R_obs_sd);      // annual observation error of juveniles

  PARAMETER_VECTOR(log_W_ret_init);    // latent true number of spawners in years 1-5
  vector<Type> W_ret_init = exp(log_W_ret_init);
  REPORT(W_ret_init);

  
  PARAMETER_VECTOR(log_alpha);         // intrinsic productivity
  vector<Type> alpha =exp(log_alpha);
  REPORT(alpha);
  
  PARAMETER_VECTOR(log_R_max);         // Asymptotic maximum recruitment
  vector<Type> R_max =exp(log_R_max);
  REPORT(R_max);
  
  PARAMETER_VECTOR(log_d);             // Third parameter in sigmoid function
  vector<Type> d = exp(log_d);
  REPORT(d); 
  
  PARAMETER_VECTOR(log_proc_sigma);    // adult-juvenile transition process error variances
  vector<Type> proc_sigma = exp(log_proc_sigma);
  REPORT(proc_sigma);
  
  PARAMETER_VECTOR(logit_proc_er_corr);//adult-juvenile transition process error correlations
  vector<Type> proc_er_corr =
    (invlogit(logit_proc_er_corr)-0.5)*2;
  REPORT(proc_er_corr)  
  
  PARAMETER_VECTOR(logit_pHOS);
  vector<Type> pHOS = invlogit(logit_pHOS);
  
  PARAMETER_VECTOR(alr_p_hyper_mu);    //hyper mean addative log ratios of adult age proportions
  REPORT(alr_p_hyper_mu);
  PARAMETER_VECTOR(log_alr_p_hyper_sigma);  // variance of alr age proportions
  vector<Type> alr_p_hyper_sigma=exp(log_alr_p_hyper_sigma);
  REPORT(alr_p_hyper_sigma);
  PARAMETER_VECTOR(logit_alr_p_hyper_cor);   //correlation of alr ages
  vector<Type> alr_p_hyper_cor= 
    (invlogit(logit_alr_p_hyper_cor)-.5)*2;
  REPORT(alr_p_hyper_cor);
  
  PARAMETER_MATRIX(alr_p_age);         //alr(age proportions)annual & juvenile life history age specific
  matrix<Type> prop_age = alr_to_simplex(alr_p_age);
  REPORT(prop_age);
  
  PARAMETER_MATRIX(logit_surv);       
  REPORT(logit_surv);
  //---------------------------------------------------------------------
  // Variables
 int Nyears =wild_S_obs.size();
 
 vector<Type> wild_return(Nyears); //total wild returns in brood year
 vector<Type> S_hat(Nyears);      //total spawners in brood year
 vector<Type> hatch_S_hat(Nyears);// wild spawners
 vector<Type> wild_S_hat(Nyears);// hatchery spawners
 
 matrix<Type> ad_LH(Nyears,4); //wild adult returns by LH & BY
 matrix<Type> ad_LH_age(Nyears,4*3); //wild adult returns by age, LH & BY

 matrix<Type> R_pred(Nyears,4); //predicted juvenile emigrants by BY
 
 matrix<Type> w_ad_age(Nyears,3);
 w_ad_age.fill(0);
 
 //---------------------------------------------------------------------
 // population model
 
 for (int Iyear=0; Iyear<Nyears;Iyear++){//brood year loop
 
     for(int Ilh=0;Ilh<4; Ilh++){ //juvenile life history loop  
      ad_LH(Iyear,Ilh)= exp(log_R_hat(Iyear,Ilh))*
                            invlogit(logit_surv(Ilh,Iyear)); //wild adult returns by LH
    
    for (int Iage=0;Iage<3; Iage++ ){//adult return age loop
      ad_LH_age(Iyear,(Ilh*3)+Iage) = prop_age(Iyear,(JLH_A(Ilh)*3)+Iage)  * 
                                               ad_LH(Iyear,Ilh); //wild adult returns by LH & age
      if(Iyear>=5){
         w_ad_age(Iyear,Iage) +=ad_LH_age(Iyear-Iage-3,(Ilh*3)+Iage);
       } 
      }
     }
    
   if(Iyear<5){
     wild_return(Iyear) = W_ret_init(Iyear);
   }else{
     wild_return(Iyear) = w_ad_age.row(Iyear).sum();
   }
   
   wild_S_hat(Iyear)=   (wild_return(Iyear) - brood_rem(Iyear));
     
   hatch_S_hat(Iyear)= wild_S_hat(Iyear)*pHOS(Iyear)/
     (1-pHOS(Iyear));
  
  S_hat(Iyear) =  wild_S_hat(Iyear)+hatch_S_hat(Iyear); //estimate total spawners
   

  //predict juveniles from transiton models. 
 for(int Ilh=0;Ilh<4; Ilh++){ //juvenile life history loop
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
 REPORT(w_ad_age);
 REPORT(ad_LH_age);
 REPORT(ad_LH);
 
 //                End of population model
 //------------------------------------------------------------------------------------
 //                Beggining of data models (likelihood)
 
 
 // Load namespace which contains the multivariate distributions
 using namespace density;
 
 //survival from emigration to adults
 PARAMETER(logit_Phi); //logit survival autocorrelation (same for all LHs)
 Type Phi = invlogit(logit_Phi);
 REPORT(Phi);
 
 PARAMETER_VECTOR(log_surv_var); // variances (unique for LHs)
 vector<Type> surv_var = exp(log_surv_var);
 REPORT(surv_var);
 
 PARAMETER_VECTOR(logit_surv_cor); //logit(correlation/2+.5) transform to bound (-1,1)
 vector<Type> surv_cor= (invlogit(logit_surv_cor)-0.5)*2;
 REPORT(surv_cor);
 
 PARAMETER(surv_alpha); //intercept 
 PARAMETER(surv_beta); //and coefficient for hyper-mean survival
 DATA_VECTOR(LH_DAYS); //average day counted for each juvenile LH, covariate for survival hyper mean linear model
 REPORT(surv_alpha);
 REPORT(surv_beta);
 
 vector<Type>  mean_logit_surv_LH(4); //hyper mean logit survival
 
 array<Type> logit_surv_demean(logit_surv.rows(), //demeaned survival
                             logit_surv.cols());
 

 
 
 for(int Ilh=0;Ilh<4; Ilh++){
   mean_logit_surv_LH(Ilh)= surv_alpha+ 
                            LH_DAYS(Ilh)*
                            surv_beta; //calculate hyper mean survival
   
   for (int Iyear=0; Iyear<Nyears;Iyear++){
   logit_surv_demean(Ilh,Iyear) = logit_surv(Ilh,Iyear)-
                                           mean_logit_surv_LH(Ilh); //calculate demeaned survival
   }
 }
 REPORT(mean_logit_surv_LH);
 REPORT(logit_surv_demean);
 
 //MAR1 likelihood
 Type like_surv = AR1(Phi, 
                      VECSCALE(UNSTRUCTURED_CORR(surv_cor),surv_var)) 
   (logit_surv_demean); // nll 
 
 //penalized complexity priors
 PARAMETER(pen_com_surv_log_sigma);
 REPORT(pen_com_surv_log_sigma);
 
   like_surv -= dnorm(surv_beta,Type(0),
                    exp(pen_com_surv_log_sigma),true);
   like_surv -= dexp(exp(pen_com_surv_log_sigma),Type(1),true);               

 
 //-----------------------------------------------------
 
 
 // age at return likelihood
 Type like_rand_age =0;   //initialize likelihood of random alr(p_age) vectors
 
 //MVN distribution with unstructured VCV matrix
 UNSTRUCTURED_CORR_t<Type> rand_age_nll(alr_p_hyper_cor);
 
 vector<Type> alp_p_age_error(4); //declare annual errors of alr age proportions

  for (int I=0; I<Nyears; I++){
   
   alp_p_age_error =
     vector<Type>(alr_p_age.row(I))-
     alr_p_hyper_mu;
   
   like_rand_age += VECSCALE(rand_age_nll,alr_p_hyper_sigma)
                    (alp_p_age_error);
 }

REPORT(alp_p_age_error);

 //------------------------------------------
 //pHOS
 Type phos_like =0;
 for (int Iyear=0; Iyear<Nyears;Iyear++){
  phos_like += dbinom_robust(hatch_carcs(Iyear),
                                total_carcs(Iyear),
                                logit_pHOS(Iyear),true);
 }
 REPORT(pHOS)
 //--------------------------------------------
 //juvenile "recruitment" MVN process error likelihood 
 Type state_Like=0;
 
   UNSTRUCTURED_CORR_t<Type> mvn(proc_er_corr);

    for(int I=0;I<Nyears; I++ ){
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
  
  
  //age composition data likelihood
  
  //PIT age at return
  Type PIT_age_like = 0;
  for (int I=0; I<PIT_age_BYs.size(); I ++){
   
       PIT_age_like+=dmultinom(vector<Type>(PIT_ages.row(I)),
                             vector<Type>( vector<Type>(prop_age.
                                              row(PIT_age_BYs(I))).
                                              segment(PIT_age_LHs(I)*3,3)),
                                              true);
  }
  
  //carcass recoveries
  
  Type carcass_age_like = 0;
  for (int Iyear=5; Iyear<Nyears; Iyear ++){
    vector<Type> prop_age_tot=vector<Type>(w_ad_age.row(Iyear))/
      Type(wild_return(Iyear));
    carcass_age_like+=dmultinom(vector<Type>(Carc_adult_age.row(Iyear)),
                                prop_age_tot,
                                              true);
  }

  
  //LL
   Type obj_fun = -state_Like-like_surv-like_rand_age
    +phos_like  +Spawner_obs_like+ Recruit_obs_like+
      PIT_age_like+carcass_age_like;
  REPORT(phos_like);
  REPORT(PIT_age_like);
  REPORT(carcass_age_like);
  REPORT(R_pred);
  REPORT(log_R_obs_sd);
  REPORT(state_Like);
  REPORT(Spawner_obs_like);
  REPORT(Recruit_obs_like);
  REPORT(obj_fun);
  REPORT(log_R_hat);
  REPORT(S_hat);
  REPORT(like_surv);
  REPORT(like_rand_age);
  return(-obj_fun);
}


