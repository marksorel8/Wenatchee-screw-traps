#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
//Data

  //p-model
    DATA_INTEGER(Ntrials);     // number of efficiency trials
    DATA_VECTOR(rel);          // number of fish released in each efficiency trial
    DATA_VECTOR(rec);          // number of fish recaptured in each efficiency trial
    DATA_MATRIX(pDesign);      // design matrix for efficiency model fixed effects
    //DATA_MATRIX(randDesign);   // design matrix for efficiency model random effects
  
  //Migrants-model
    DATA_INTEGER(N_trap);      // numer of daily trap catch observations
    DATA_INTEGER (max_day);    // maximum day of trap catch
    DATA_INTEGER (Nyears);    // number of time series
    
    DATA_IVECTOR(seriesFac);    // indicator variable for which time series (i.e. year) a catch is from
    DATA_VECTOR(Catch);       // daily trap catch observations
    DATA_IVECTOR(catch_DOY);   // day of daily trap catch observations
    DATA_MATRIX(pDat);         // fixed effects design matrix for predicting daily efficiency
    //DATA_MATRIX(pDatRand);     // random effects design matrix for predicting daily efficiency
  
  
//Parameters

    PARAMETER_VECTOR(pCoefs);     // coefficients in efficiency model
    PARAMETER_VECTOR(mu_M);              // intercept of AR(1) process for log-mean daily outmigrants
    PARAMETER_VECTOR(phi_M);             // AR(1) coefficient for log-mean daily outmigrants (to be bounded from -1 to 1)
    PARAMETER_VECTOR(log_sigma_M);       // log AR(1) process error SD for log-mean daily outmigrants
    PARAMETER_MATRIX(log_M_hat_z);// log-mean of daily outmigrants (z-scores)
  
  
//Likelihood

    //efficiency model
    
      Type LikeFixP = 0;            //declate log likelihood of efficiency-trial data
  
      // Binomial logistic regression for trap efficiency model
      vector<Type> fixed = pDesign * pCoefs;
      vector<Type> logitPred = fixed ;//+eta;
      
      // efficiency data likelihood
      LikeFixP =dbinom_robust(rec,rel,logitPred,true).sum();
      REPORT(LikeFixP);

      
    //AR(1) for Z-scored daily migrants
      
      vector<Type> sigma_M=exp(log_sigma_M);   //transform process SD's
      matrix<Type> log_M_hat(Nyears,max_day);  //Declare matrix of log(expected migrants) un-z-scored
      
      // Load namespace which contains the multivariate distributions
      using namespace density;
      
      Type likeMHatZ=0;                       //declare likelihood component for AR1 model
      for( int Iyear=0; Iyear<Nyears; Iyear++){
        likeMHatZ+= AR1(phi_M(Iyear))(log_M_hat_z.row(Iyear)); //AR1 likelihood
        
        log_M_hat.row(Iyear)=(vector<Type> (log_M_hat_z.row(Iyear)) *
          Type (sigma_M(Iyear))) +
          Type(mu_M(Iyear)) ; //un-z-score (multiply by SD and add mean)
       
        
      }
      
       REPORT(likeMHatZ);
        REPORT(log_M_hat_z);
        REPORT(log_M_hat);
        matrix<Type> M_hat = exp(log_M_hat.array());          //exponentiate to get expected migrant estimate
        
      
    //Observation likelihood   
  
      //predict daily trap efficiency for all days with catch
      vector<Type>  logit_p  =  pDat * pCoefs;    // fixed component of logit trap efficiency for each catch day   

    
      //calculate expected catches( # migrants * trap efficiency)
      vector<Type> m_c_day(N_trap);     //vector of expected migrants on days with catch observations 
      for( int Iday=0; Iday<N_trap; Iday++){
         m_c_day(Iday)=M_hat(seriesFac(Iday)-1,catch_DOY(Iday)-1);
          }
      vector<Type> C_hat = m_c_day * invlogit(logit_p);    //expected catch
       
      Type LikeCatch =dpois(Catch,C_hat,true).sum();      //catch data likelihood
     
     
    //objective function
      Type obj_fun = LikeFixP + LikeCatch -likeMHatZ ;//
         
     
     
    //Derived quanitities 
      //vector<Type> yrlngs = log_M_hat.tail(115);  //daily yearling migrant vecotr
      //Type subTot = sum(exp(log_M_hat));               //total subyearling migrants
      //Type yrlngTot = sum(exp(yrlngs));           //total yearling migrants
    
      //process variance of daily numbr of migrants AR(1)
      // vector<Type> var_m=sigma_M*sigma_M;
      // 
      // 
      // 
      // //fry
      // matrix<Type> fry_block = M_hat.block(0,0,Nyears,90);
      // //REPORT(fry_block); // block starting (0,0) taking 1 row & 1 col
      // vector<Type> fryTot = fry_block.rowwise().sum();
      // 
      // //summer parr
      // 
      //  matrix<Type> sumParr_block = M_hat.block(0,90,Nyears,125);
      // // REPORT(sumParr_block); // block starting (0,0) taklength(87:)ing 1 row & 1 col
      // vector<Type> sumParrTot = sumParr_block.rowwise().sum();     
      // //fall parr
      //  matrix<Type> fallParr_block = M_hat.block(0,215,Nyears,130);
      // // REPORT(fallParr_block); // block starting (0,0) taking 1 row & 1 col
      //  vector<Type> fallParrTot = fallParr_block.rowwise().sum();
      // 
      // //smolt
      // matrix<Type> smolt_block = M_hat.block(0,345,Nyears,121);
      // // REPORT(fallParr_block); // block starting (0,0) taking 1 row & 1 col
      // vector<Type> smoltTot = smolt_block.rowwise().sum();
      // 
      // REPORT(fryTot);
      // REPORT(sumParrTot);
      // REPORT(fallParrTot);
      // REPORT(smoltTot);
      
 //Reporting
     // REPORT(LikeFixP);
     //REPORT(LikeRanP);
     //ADREPORT(yrlngTot);
     // REPORT(LikeCatch);
     // REPORT(likeMHatZ);
     // REPORT(obj_fun);
      REPORT(mu_M);
      REPORT(phi_M);
     //REPORT(var_m);
    //REPORT(pCoefs);
     // ADREPORT(logitPred);
      REPORT(logit_p);
     // ADREPORT(mu);
     REPORT(C_hat);
     REPORT(log_M_hat);
     ADREPORT(M_hat);
     return(-obj_fun);
}

