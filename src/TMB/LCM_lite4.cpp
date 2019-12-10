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
  DATA_MATRIX(surv_dat);
  DATA_INTEGER(Nproj);
  
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
  
  
  PARAMETER_VECTOR(logit_pHOS);
  vector<Type> pHOS = invlogit(logit_pHOS);
  
  PARAMETER_VECTOR(alr_p_hyper_mu);    //hyper mean addative log ratios of adult age proportions
  REPORT(alr_p_hyper_mu);
  PARAMETER_VECTOR(log_alr_p_hyper_sigma);  // variance of alr age proportions
  vector<Type> alr_p_hyper_sigma=exp(log_alr_p_hyper_sigma);
  REPORT(alr_p_hyper_sigma);
  PARAMETER_VECTOR(logit_alr_p_hyper_cor);   //correlation of alr ages
  
  PARAMETER_MATRIX(alr_p_age);         //alr(age proportions)annual & juvenile life history age specific
  matrix<Type> prop_age = alr_to_simplex(alr_p_age);
  REPORT(prop_age);
  REPORT(alr_p_age);
  
  PARAMETER_VECTOR(logit_surv);       
  REPORT(logit_surv);
  

  
  PARAMETER(surv_alpha); //intercept 
  PARAMETER(surv_beta); //and coefficient for hyper-mean survival
  DATA_VECTOR(LH_DAYS); //average day counted for each juvenile LH, covariate for survival hyper mean linear model
  REPORT(surv_alpha);
  REPORT(surv_beta);
  
    //survival from emigration to adults
  PARAMETER(logit_Phi); //logit survival autocorrelation (same for all LHs)
  Type Phi = invlogit(logit_Phi);
  REPORT(Phi);
  
  PARAMETER(log_surv_var); // variances (unique for LHs)
  Type surv_var = exp(log_surv_var);
  REPORT(surv_var);
  
  
  
 
  
  //---------------------------------------------------------------------
  // Variables
 int Nyears =wild_S_obs.size();
 vector<Type>  mean_logit_surv_LH(4); //hyper mean logit survival
  
  for(int Ilh=0;Ilh<4; Ilh++){
    mean_logit_surv_LH(Ilh)= surv_alpha+ 
      LH_DAYS(Ilh)*
      surv_beta;
  }
  
  matrix<Type> logit_surv_lh_yr(4,logit_surv.size());
  
  // array<Type> surv_var_ar(logit_surv.rows(), //demeaned survival
  //                          logit_surv.cols());
  
  
  for(int Ilh=0;Ilh<4; Ilh++){
    for (int Iyear=0; Iyear<Nyears-3;Iyear++){
      logit_surv_lh_yr(Ilh,Iyear) = logit_surv(Iyear)+
        mean_logit_surv_LH(Ilh); //calculate demeaned survival
      //surv_var_ar(Ilh,Iyear) =surv_var(Ilh);
    }
  }
  REPORT(mean_logit_surv_LH);
  REPORT(logit_surv_lh_yr);
  
  
 vector<Type> wild_return(Nyears+Nproj); //total wild returns in brood year
 vector<Type> S_hat(Nyears+Nproj);      //total spawners in brood year
 vector<Type> hatch_S_hat(Nyears+Nproj);// wild spawners
 vector<Type> wild_S_hat(Nyears+Nproj);// hatchery spawners
 
 matrix<Type> ad_LH(Nyears+Nproj,4); //wild adult returns by LH & BY
 matrix<Type> ad_LH_age(Nyears+Nproj,4*3); //wild adult returns by age, LH & BY

 matrix<Type> R_pred(Nyears+Nproj,4); //predicted juvenile emigrants by BY
 
 matrix<Type> w_ad_age(Nyears+Nproj,3);
 w_ad_age.fill(0);
 
 //---------------------------------------------------------------------
 // population model
 
 //calculate adult survivals from eachjuvenile life history 
 for (int Iyear=0; Iyear<(Nyears-3);Iyear++){//brood year loop
 
     
     for(int Ilh=0;Ilh<4; Ilh++){ //juvenile life history loop  
      ad_LH(Iyear,Ilh)= exp(log_R_hat(Iyear,Ilh))*
                            invlogit(logit_surv_lh_yr(Ilh,Iyear)); //wild adult returns by LH
    
    for (int Iage=0;Iage<3; Iage++ ){//adult return age loop
      ad_LH_age(Iyear,(Ilh*3)+Iage) = prop_age(Iyear,(JLH_A(Ilh)*3)+Iage)  * 
                                               ad_LH(Iyear,Ilh); //wild adult returns by LH & age
      
   
      if(Iyear>=2){
         w_ad_age(Iyear+3,Iage) +=ad_LH_age(Iyear-Iage,(Ilh*3)+Iage);
       } 
      }
     }
 }
 

//calculate total annual spawners (wild and hatchery) 


for (int Iyear=0; Iyear<Nyears;Iyear++){//brood year loop
   if(Iyear<5){
     wild_return(Iyear) = W_ret_init(Iyear);
   }else{
     wild_return(Iyear) = w_ad_age.row(Iyear).sum();
   }
   
   wild_S_hat(Iyear)=   (wild_return(Iyear) - brood_rem(Iyear));
     
   hatch_S_hat(Iyear)= wild_S_hat(Iyear)*pHOS(Iyear)/
     (1-pHOS(Iyear));
  
  S_hat(Iyear) =  wild_S_hat(Iyear)+hatch_S_hat(Iyear); //estimate total spawners


  //predict juveniles from transition models. 
  
 for(int Ilh=0;Ilh<4; Ilh++){ //juvenile life history loop
 if(Model(Ilh)==1){
   R_pred(Iyear,Ilh)= R_max(Ilh)*(1-exp(-(pow(S_hat(Iyear)/alpha(Ilh),d(Ilh))))); //weibull CDF

 }

 if(Model(Ilh)==2){
   R_pred(Iyear,Ilh)= (alpha(Ilh)*S_hat(Iyear))/(1+alpha(Ilh)*S_hat(Iyear)/R_max(Ilh)); //beverton-holt model
    }
   }//LH loop
  }//year loop
 
 
 //                End of population model
 //------------------------------------------------------------------------------------
 //                Beggining of data models (likelihood)
 
 
 // Load namespace which contains the multivariate distributions
 using namespace density;
 
 
 
 //PARAMETER_VECTOR(logit_surv_cor); //logit(correlation/2+.5) transform to bound (-1,1)
 //vector<Type> surv_cor =invlogit(logit_surv_cor)*2 -1;
 //vector<Type> unconstrained_params(logit_surv_cor);  // Dummy parameterization of correlation matrix
 //matrix<Type> Sigma = UNSTRUCTURED_CORR(logit_surv_cor).cov();
 // matrix<Type> Sigma (4,4);
 //  int cnter=0;
 //  for(int I=0;I<4;I++){
 //   Sigma(I,I)=1;
 //   if(I>0){
 //   for(int J=0;J<I;++J){
 //    Sigma(I,J)=surv_cor(cnter);
 //     Sigma(J,I)=Sigma(I,J);
 //     cnter++;
 //   }}
 // }
 
 
 //vector<Type> surv_cor= (invlogit(logit_surv_cor)-0.5)*2;
 ///REPORT(surv_cor);
 Type like_surv = SCALE(AR1(Phi),surv_var)(logit_surv);

 //MAR1 likelihood
 // MVNORM_t<Type> surv_MVNORM(Sigma);
 // Type like_surv = AR1(Phi, VECSCALE(surv_MVNORM,surv_var))
 //   (logit_surv_demean); // nll 
 
 
 //REPORT(Sigma);   
 //penalized complexity priors
 // PARAMETER(pen_com_surv_log_sigma);
 // REPORT(pen_com_surv_log_sigma);
 // 
 //   like_surv -= dnorm(surv_beta,Type(0),
 //                    exp(pen_com_surv_log_sigma),true);
 //   like_surv -= dexp(exp(pen_com_surv_log_sigma),Type(1),true);               
 // 
 
 ///Data likelihood
 // DATA_IVECTOR(surv_yr)
 // for (int I = 0; I<surv_dat.rows();I++){
 //   Type logit_surv_i = Type(surv_alpha)+
 //     Type(Type(surv_dat(I,1))*
 //     Type(surv_beta))+Type(logit_surv(surv_yr(I)));
 //   like_surv -= dbinom_robust(surv_dat(I,0),Type(1),logit_surv_i,true);
 // }
 
 
 
 //-----------------------------------------------------
 
 
 // age at return likelihood
 Type like_rand_age =0;   //initialize likelihood of random alr(p_age) vectors
 vector <Type> alr_p_age_error(4);
 //MVN distribution with unstructured VCV matrix
UNSTRUCTURED_CORR_t<Type> rand_age_nll(logit_alr_p_hyper_cor);

   //for ( int Iage = 0; Iage<2;Iage++){
 //vector<Type> alr_p_age_error(4); //declare annual errors of alr age proportions
 // matrix<Type> cov(2,2);
 // cov(0,0)= alr_p_hyper_sigma(Iage*2);  //cov param
 // cov(1,1)= alr_p_hyper_sigma(Iage*2+1);
 // cov(0,1)=invlogit(logit_alr_p_hyper_cor(Iage))*2-1;
  
  for (int I= 0; I< (Nyears-3); I++){
    alr_p_age_error =
      vector<Type>( vector<Type>(alr_p_age.row(I))-
     vector<Type>(alr_p_hyper_mu));

    //MVNORM_t<Type> rand_age_nll(cov);

   like_rand_age += //rand_age_nll(alr_p_age_error);
     VECSCALE(rand_age_nll,
                              vector<Type>(alr_p_hyper_sigma*alr_p_hyper_sigma))
                    (alr_p_age_error);
 }

 //}
 
matrix<Type> alr_cov =rand_age_nll.cov();
REPORT(alr_cov);
 
  
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
 
   UNSTRUCTURED_CORR_t<Type> mvn(logit_proc_er_corr);

    for(int I=0;I<Nyears; I++ ){
      state_Like+= VECSCALE(mvn,
                            vector<Type>(proc_sigma*proc_sigma))
      (vector<Type>( vector <Type> (log_R_hat.row(I))- 
        log(vector<Type>(R_pred.row(I)))));
    }
    
  matrix<Type> proc_er_corr =mvn.cov();
  REPORT(proc_er_corr)  
 // matrix<Type> resids = log(R_hat)-log(R_pred);
  //REPORT(resids);
  
   Type Spawner_obs_like= sum(dnorm(log(vector<Type>(S_obs)),//.segment(5,(Nyears-5)))),
                                     log(vector<Type>(S_hat.segment(0,S_obs.size()))),//.segment(5,(Nyears-5)))),
                                      S_obs_cv,
                                      true));    //observation likelihood (spawners)

  //prior of S_obs_CV from Murdoch et al 2019
  Spawner_obs_like+=dnorm(S_obs_cv,Type(0.055),Type(0.02),true);
  
  Type Recruit_obs_like= 0;                   //Initialize observation likelihood (recruits)
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

  


  
  SIMULATE{
    vector<Type> surv_proj(Nproj);
    //survival year effects
    surv_proj(0)=logit_surv(Nyears-4)*Phi + sqrt(1-pow(Phi,2))*rnorm(Type(0),Type(surv_var));
    for ( int  I=1; I<Nproj; I++)
      surv_proj(I)=surv_proj(I-1)*Phi + sqrt(1-pow(Phi,2))*rnorm(Type(0),Type(surv_var));
    
    REPORT(surv_proj);
    
    
    //age proportions
    matrix <Type> alr_age_proj(Nproj,4);//matrix to hold props
    
    MVNORM_t<Type> rand_age_mvn(alr_cov);//mvn disctribution to simmulate from
      for ( int  I=0; I<Nproj; I++)
        alr_age_proj.row(I)=rand_age_mvn.simulate()+alr_p_hyper_mu;
      matrix<Type> prop_age_proj = alr_to_simplex(alr_age_proj);
    REPORT(prop_age_proj);
    
  
  
  //spawner->juvenile process error
  matrix <Type> juv_proc_er_proj(Nproj,4);//matrix to hold props
  
  MVNORM_t<Type> juv_proc_er_mvn(mvn);//mvn disctribution to simmulate from
  for ( int  I=0; I<Nproj; I++)
    juv_proc_er_proj.row(I)=juv_proc_er_mvn.simulate();
   REPORT(juv_proc_er_proj);
  
  //vector to hold projection pHOS
  vector<Type> pHOS_proj(Nproj);
  pHOS_proj.fill(0);
  //vector to hold projection brood_rem
  vector<Type> brood_rem_proj(Nproj);
  brood_rem_proj.fill(0);
  //calculate adult contributions from Brood years 2015-2017 with juveniles from fitted model in those years and projected survival and age comp
  //then calculate spawners in brood year 2018-2020
  // is it ok to use fitted value for smolts from 2017 brood even though no data?
  for(int Iyear=(Nyears-3); Iyear < int( Nyears); Iyear++){//brood year loop

    for(int Ilh=0;Ilh<4; Ilh++){ //juvenile life history loop

      ad_LH(Iyear,Ilh)= exp(log_R_hat(Iyear,Ilh))*
        invlogit(surv_proj(Iyear-Nyears+3)+mean_logit_surv_LH(Ilh)); //wild adult returns by LH

      for (int Iage=0;Iage<3; Iage++ ){//adult return age loop
        ad_LH_age(Iyear,(Ilh*3)+Iage) = prop_age_proj(Iyear-Nyears+3,(JLH_A(Ilh)*3)+Iage)  *
          ad_LH(Iyear,Ilh); //wild adult returns by LH & age

      }//end of adult return age loop
    }//end of life history loop
  
  //calculate spawners
  wild_return(Iyear) = w_ad_age.row(Iyear).sum();
  wild_S_hat(Iyear)=   (wild_return(Iyear) - brood_rem_proj(Iyear-Nyears));//make brood rem
  hatch_S_hat(Iyear)= wild_S_hat(Iyear)*pHOS_proj(Iyear-Nyears);//make pHOS_proj
      (1-pHOS(Iyear));
  S_hat(Iyear) =  wild_S_hat(Iyear)+hatch_S_hat(Iyear);
  }
  
  
  //loop through years Nyears+1 - Nyears+1+Nproj-3
  //calc juveniles and spawners
  
  //loop through projection years and calculate juveniles then adults for 2018:2042
  for ( int Iyear = Nyears; Iyear<Nyears+Nproj-3; Iyear++){ //start in brood year 2018

    for(int Ilh=0;Ilh<4; Ilh++){ //juvenile life history loop
          
          
          //predict juveniles
          if(Model(Ilh)==1){
            R_pred(Iyear,Ilh)= R_max(Ilh)*(1-exp(-(pow(S_hat(Iyear)/alpha(Ilh),d(Ilh))))); //weibull CDF

          }

          if(Model(Ilh)==2){
            R_pred(Iyear,Ilh)= (alpha(Ilh)*S_hat(Iyear))/(1+alpha(Ilh)*S_hat(Iyear)/R_max(Ilh)); //beverton-holt model
          }
        
        //add process error
        log_R_hat(Iyear,Ilh)=R_pred(Iyear,Ilh)+juv_proc_er_proj(Iyear-Nyears,Ilh);
        
    
        ad_LH(Iyear,Ilh)= exp(log_R_hat(Iyear,Ilh))*
          invlogit(surv_proj(Iyear-Nyears+3)+mean_logit_surv_LH(Ilh)); //wild adult returns by LH
          
          for (int Iage=0;Iage<3; Iage++ ){//adult return age loop
            ad_LH_age(Iyear,(Ilh*3)+Iage) = prop_age_proj(Iyear-Nyears+3,(JLH_A(Ilh)*3)+Iage)  *
              ad_LH(Iyear,Ilh); //wild adult returns by LH & age
            
          }//end of adult return age loop
        
        
        }//end of LH loop

    wild_return(Iyear) = w_ad_age.row(Iyear).sum();
    wild_S_hat(Iyear)=   (wild_return(Iyear) - brood_rem_proj(Iyear-Nyears));//make brood rem
    hatch_S_hat(Iyear)= wild_S_hat(Iyear)*pHOS_proj(Iyear-Nyears);//make pHOS_proj
      (1-pHOS(Iyear));
    S_hat(Iyear) =  wild_S_hat(Iyear)+hatch_S_hat(Iyear);
    
  }
  
  for ( int Iyear = Nyears+Nproj-3; Iyear<Nyears+Nproj; Iyear++){ //start in brood year 2018
    wild_return(Iyear) = w_ad_age.row(Iyear).sum();
    wild_S_hat(Iyear)=   (wild_return(Iyear) - brood_rem_proj(Iyear-Nyears));//make brood rem
    hatch_S_hat(Iyear)= wild_S_hat(Iyear)*pHOS_proj(Iyear-Nyears);//make pHOS_proj
      (1-pHOS(Iyear));
    S_hat(Iyear) =  wild_S_hat(Iyear)+hatch_S_hat(Iyear);
    
  }
  
  
  
  }
  
  
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
  REPORT(wild_return);
  REPORT(hatch_S_hat);
  REPORT(wild_S_hat);
  REPORT(w_ad_age);
  REPORT(ad_LH_age);
  REPORT(ad_LH);
  
  
  return(-obj_fun);
}


