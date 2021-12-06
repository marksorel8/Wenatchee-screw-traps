#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  
//Data
    DATA_INTEGER(subyearlings); //flag indicating whether this is a model of subyearlings of yearlings
  
  //p-model
    DATA_VECTOR(rel);          // number of fish released in each efficiency trial
    DATA_VECTOR(rec);          // number of fish recaptured in each efficiency trial
    DATA_IVECTOR(efficIndex);  // index of catch data that each efficiency trial corresponds with 

  
  //Migrants-model
    DATA_INTEGER (N_trap);    // numer of daily trap catch observations
    DATA_INTEGER (N_day);     // maximum day of trap catch
    DATA_INTEGER (Nyears);    // number of time series (years)

    DATA_IVECTOR(seriesFac);  // indicator variable for which time series (i.e. year) a catch is from
    DATA_VECTOR(Catch);       // daily trap catch observations
    DATA_IVECTOR(catch_DOY);  // day of year of daily trap catch observations
    DATA_INTEGER(first_DOY);  // first DOY of emigrant estimates, used for calculating standardized migration mindows accross streams
    DATA_MATRIX(pDat);        // fixed effects design matrix for daily trap capture efficiency
    
    DATA_INTEGER(Use_NB);     // flag indicating whether to use Negative Binomial distribution, as opposed to Poisson, for the daily catch
 
    DATA_IVECTOR(breaks)      // day of year breaks for aggregating emigrants accross days within "windows" 
                              // based on post-hoc analysis of average timing with mixture diistribution.
  
    DATA_SCALAR(stream_length);      // lengh (km) of natal stream surveyed for spawners, for standardizing emigrant estimates
  
//Parameters

   
    PARAMETER_VECTOR(mu_M);   // intercept of log-mean daily outmigrants
    // PARAMETER(logit_phi_d);   // logit AR(1) coefficient for across-year errors in log-mean daily emigrants
    PARAMETER(logit_phi_e);   // logit AR(1) coefficient for year-specific erros in log-mean daily outmigrants
    PARAMETER(ln_tau_d);      // log AR(1) process error precission for across year erros in log-mean daily emigrants
    PARAMETER(ln_tau_e);      // log AR(1) process error precission for year-specific errors in log-mean daily emmigrants 
    PARAMETER_VECTOR(pCoefs); // coefficients in trap capture efficiency model
    PARAMETER(log_phi_NB);    // log -transformed overdispersion paramater for negative binomial 

   
//Random Effects   
    PARAMETER_VECTOR(delta);   // accross-year errors of log-means of daily emigrants 
    PARAMETER_MATRIX(epsilon); // year-specific errors of log-means of daily emigrants
     
     
    
//objective function
using namespace density;       // load "package" with AR1 distribtion
vector<Type> jnll_comp(3);     // declare vector of 3 likelihood componenets (AR1 errors, efficency trials, and catch data)
jnll_comp.setZero();           // set likelihood componenets to zero

//Variables
// Type phi_d=invlogit(logit_phi_d);// transform correlation coefficient for across-year daily outmigrant errors to "0-1" space. I don't thing a negative correlation would ever be the case here, and must be >1 to be stationar
Type phi_e=invlogit(logit_phi_e);// transform  correlation coefficient for yearly daily outmigrant errors to "0-1" space. I don't thing a negative correlation would ever be the case here, and must be >1 to be stationary
Type sigma_d = 1/exp(ln_tau_d);  // transorm log precision to standard deviation 
Type sigma_e = 1/exp(ln_tau_e);  // transorm log precision to standard deviation 
Type phi_NB=exp(log_phi_NB);     // transform negative binomial "phi" paramater 

//-----------------------------------------------------------------
//-----------------------------------------------------------------
//emigrant process
//-----------------------------------------------------------------

//likelihood of random effects (errors)

//random walk for across-year daily erros
////first time step (assuming value at time 1 is 0)
jnll_comp(0)-= dnorm(Type(delta(0)),Type(0),sigma_d,true);
////subsequent time steps
for(int i =1; i<(N_day-1); i++){
  jnll_comp(0)-= dnorm(Type(delta(i)-delta(i-1)),Type(0),sigma_d,true);
  }

vector<Type> delta2(N_day);
delta2 << Type(0),delta;

//AR1 likelihood for across-year daily errors (representing the daily deviation from the mean in a "average year")
// delta[1] ~ N(0,sigma_e)
// delta[t] ~ phi_d*delta[t-1] + sqrt(1-phi_d^2)*e[t], e[t] ~ N(0,sigma_e)
// jnll_comp(0)+=SCALE(AR1(phi_d),
          // sigma_d)(delta); //uses AR1 likelihood from density namespace of TMB, which is fast. Returns negative log likelihood.
                           //SCALE function allows us to set the standard deviation of the AR1 marginal distribution

//AR1 likelihood for year-specific daily errors (representing daily deviations on top the average within a given year)
// epsilon[1] ~ N(0,sigma_e)
// epsilon[t] ~ phi_e*epsilon[t-1] + sqrt(1-phi_e[y]^2)*e[t,y], e[t,y] ~ N(0,sigma_e[y])
for( int Iyear=0; Iyear<(Nyears); Iyear++){     //loop through years
   jnll_comp(0)+= SCALE(AR1(phi_e),
             sigma_e)(epsilon.col(Iyear)); 
  }


//calculated expected emigrants abundance (M_hat)

//M_hat[t,y] = exp(mu[y]+ delta[t] + epsilon[t,y])
matrix<Type> M_hat(N_day,Nyears);                //declare matrix to hold expected emigrant abundances
for( int Iyear=0; Iyear<Nyears; Iyear++){        //loop through years
  M_hat.col(Iyear)=
    vector<Type>(exp(
        Type(mu_M(Iyear)) +                      // annual mean 
          delta2 +                                // across-year daily error
          vector<Type> (epsilon.col(Iyear))));   // year-specific daily error
  
}


// jnll_comp(0)-=(dexp(sigma_d,Type(1),true)+ln_tau_d);
// jnll_comp(0)-=(dexp(sigma_e,Type(1),true)+ln_tau_e);

//-----------------------------------------------------------------
//-----------------------------------------------------------------
//Trapping process
//-----------------------------------------------------------------

//capture efficiency (p) on all days with catch  data
vector<Type>  logit_p  =  pDat * pCoefs;    //logit trap efficiency for each catch day  
                                            //(design matrix %*% coefficient vector) 
 
//likelihood of "trap capture efficiency" data
// recaps[t,y] ~ binomial( releases[t,y], p[t,y])
for (int I =0; I<rel.size(); I++){          //loop over all release-recapture experiments
  jnll_comp(1)-= dbinom_robust(rec(I),      
            rel(I),
            logit_p(efficIndex(I)),true);   
}


// likleihood of catch data 

// if using negative binomial distribution
// observed_catch[t,y] ~ negbinom(mu = expected_catch[t,y], 
//                                var = expected_catch[t,y]+(expected_catch[t,y]^2/p_NB))
// expected_catch[t,y]= M_hat[t,y] * p[t,y]

// if using Poisson distribution
// observed catch[t,y] ~ Poisson(expected_catch[t,y])
for( int Iday=0; Iday<N_trap; Iday++){ // loop over all days with catch data
  Type expected_catch = Type(M_hat(catch_DOY(Iday),seriesFac(Iday))*
    invlogit(logit_p(Iday)));
  if(Use_NB){                          // if using negative binomial distribution
    jnll_comp(2) -= dnbinom2 (Catch(Iday), expected_catch,
              expected_catch+(pow(expected_catch,2)/phi_NB), true);
    
  }else{                               // otherwise use Poisson observation likelihood
    jnll_comp(2) -= dpois(Catch(Iday),expected_catch, true);
  }
}


//objective function
      Type obj_fun = sum(jnll_comp);         //sum likelihood components

//Derived quantities

//Calculate sums of emigrants over days within each "migration window" for a given life history
      int columns = 1;                             // For yearlings, number of life histories = 1 (smolts) 
      if(subyearlings) columns = 3;                // For subyearling migrants, number of life histories =  3 (fry, summer parr, and fall parr) 
      matrix<Type> LH_sums(Nyears,columns);        // initialize matrix to hold sums of daily migrants over periods
      vector<Type> temp(N_day);                    //initialize vector to hold daily counts for a given year

               for( int Iyear=0; Iyear<Nyears; Iyear++){ //loop through years
           temp=M_hat.col(Iyear);                  // get daily counts for a given year
           if(subyearlings){               
             LH_sums (Iyear,0)=temp.segment(0,(breaks[0]-first_DOY)).sum();             //sum over fry migration window
             LH_sums (Iyear,1)=temp.segment((breaks[0]-first_DOY),(breaks[1]-breaks[0])).sum();       //sum over summer parr migration window;
             LH_sums (Iyear,2)=temp.segment((breaks[1]-first_DOY),(N_day-(breaks[1]-first_DOY))).sum();} //sum over fall parr migration window
           else{
             LH_sums(Iyear,0)=temp.sum();}                           //sum over all days if a yearling model
         }
  
         LH_sums=log(LH_sums.array()/stream_length);                               //take log of sums for defining lognormal posterior distribution using the delta method
         
//average log(emigrants) per day accross years         
         matrix<Type> log_M_hat= log(M_hat.array());                 //convert daily emigrants to log(emigrants)
         vector<Type> mean_day_log_M = log_M_hat.rowwise().mean();   //calculate mean log(emigrants) per day accross years 
         
         
//Reporting
     REPORT(M_hat);                // daily emigrants[t,y]
     ADREPORT(LH_sums);            //log sums of emigrants within each emigrant life history type and year
     ADREPORT(mean_day_log_M);     //average log daily emigrants
      return(obj_fun);
}
