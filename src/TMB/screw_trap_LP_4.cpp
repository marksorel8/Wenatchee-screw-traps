#include <TMB.hpp>


//  This model estimates the daily number of juvenile emigrants past screw traps
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

// Eric Buhle contributed substantially to the conception of this model 
// and inspiration was also drawn from BTSPAS. All code was written by Mark Sorel
// except for code chunks that were copied directly from the glmmTMB source code  
// <https://cran.r-project.org/web/packages/glmmTMB/index.html> 
// with permission form some glmmTMB coauthors. 


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Functions copied from glmmTMB to calculate random effects probabilities. 
// On the r-side, I specify a mixed effects model using lme4 notation and
// use glmmTMB functions to generate the design matrices and covariance structures,
// which I then feed into the TMB model for optimization.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//valid covariance structures
enum valid_covStruct {
  diag_covstruct = 0,
  us_covstruct   = 1,
  cs_covstruct   = 2
};

//defines elements of list data structure
template <class Type>
struct per_term_info {
  // Input from R
  int blockCode;     // Code that defines structure
  int blockSize;     // Size of one block
  int blockReps;     // Repeat block number of times
  int blockNumTheta; // Parameter count per block
  matrix<Type> dist;
  vector<Type> times;// For ar1 case
  // Report output
  matrix<Type> corr;
  vector<Type> sd;
};


//translates r list data structure to C/TMB list data structure. 
//Returns a vector of list, where each  list is for one random effect "component" i.e. inside one parenthesid e.g.(LH|Year)
template <class Type>
struct terms_t : vector<per_term_info<Type> > {
  terms_t(SEXP x){
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP y = VECTOR_ELT(x, i);    // y = x[[i]]
      int blockCode = (int) REAL(getListElement(y, "blockCode", &isNumericScalar))[0];
      int blockSize = (int) REAL(getListElement(y, "blockSize", &isNumericScalar))[0];
      int blockReps = (int) REAL(getListElement(y, "blockReps", &isNumericScalar))[0];
      int blockNumTheta = (int) REAL(getListElement(y, "blockNumTheta", &isNumericScalar))[0];
      (*this)(i).blockCode = blockCode;
      (*this)(i).blockSize = blockSize;
      (*this)(i).blockReps = blockReps;
      (*this)(i).blockNumTheta = blockNumTheta;
      // Optionally, pass time vector:
      SEXP t = getListElement(y, "times");
      if(!isNull(t)){
        RObjectTestExpectedType(t, &isNumeric, "times");
        (*this)(i).times = asVector<Type>(t);
      }
      // Optionally, pass distance matrix:
      SEXP d = getListElement(y, "dist");
      if(!isNull(d)){
        RObjectTestExpectedType(d, &isMatrix, "dist");
        (*this)(i).dist = asMatrix<Type>(d);
      }
    }
  }
};


//function that calculates the probability of random effects for many different random effects structures.
//Returns negative log prob of random effects for a given random effect compnenet e.g. (LH|year)
template <class Type>
Type termwise_nll(array<Type> &U, vector<Type> theta, per_term_info<Type>& term, bool do_simulate = false) {
  Type ans = 0;
  if (term.blockCode == diag_covstruct){
    // case: diag_covstruct
    vector<Type> sd = exp(theta);
    
    for(int i = 0; i < term.blockReps; i++){
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), sd, true).sum();
      if (do_simulate) {
        U.col(i) = rnorm(Type(0), sd);
      }
    }
    term.sd = sd; // For report
  }
  
  else if (term.blockCode == us_covstruct){
    // case: us_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    vector<Type> corr_transf = theta.tail(theta.size() - n);
    vector<Type> sd = exp(logsd);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == cs_covstruct){
    // case: cs_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    Type corr_transf = theta(n);
    vector<Type> sd = exp(logsd);
    Type a = Type(1) / (Type(n) - Type(1));
    Type rho = invlogit(corr_transf) * (Type(1) + a) - a;
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        corr(i,j) = (i==j ? Type(1) : rho);
    density::MVNORM_t<Type> nldens(corr);
    density::VECSCALE_t<density::MVNORM_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else error("covStruct not implemented!");
  return ans;
}


//function that creates the structures and call termwise_nll for all random effects. 
//Returns negative log prob of random effects.
template <class Type>
Type allterms_nll(vector<Type> &u, vector<Type> theta,
                  vector<per_term_info<Type> >& terms,
                  bool do_simulate = false) {
  Type ans = 0;
  int upointer = 0;
  int tpointer = 0;
  int nr, np = 0, offset;
  for(int i=0; i < terms.size(); i++){
    nr = terms(i).blockSize * terms(i).blockReps;
    // Note: 'blockNumTheta=0' ==> Same parameters as previous term.
    bool emptyTheta = ( terms(i).blockNumTheta == 0 );
    offset = ( emptyTheta ? -np : 0 );
    np     = ( emptyTheta ?  np : terms(i).blockNumTheta );
    vector<int> dim(2);
    dim << terms(i).blockSize, terms(i).blockReps;
    array<Type> useg( &u(upointer), dim);
    vector<Type> tseg = theta.segment(tpointer + offset, np);
    ans += termwise_nll(useg, tseg, terms(i), do_simulate);
    upointer += nr;
    tpointer += terms(i).blockNumTheta;
  }
  return ans;
}

//End of glmmTMB functions

//Beggining of objective function. Returns the joint negative log likelihood
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  //Data
    ////Migrants-model
  DATA_INTEGER (N_trap);    // numer of daily trap catch observations
  DATA_INTEGER (N_day);     // maximum day of trap catch
  DATA_INTEGER (Nyears);    // number of time series (years)
  
  DATA_IVECTOR(seriesFac);  // index for which time series (i.e. year) a catch is from
  DATA_VECTOR(Catch);       // daily trap catch observations
  DATA_IVECTOR(catch_DOY);  // day of year of daily trap catch observations
  DATA_INTEGER(Use_NB);     // flag indicating whether to use Negative Binomial distribution, as opposed to Poisson, for the daily catch
  
  ////p-model
  DATA_VECTOR(rel);          // number of fish released in each efficiency trial
  DATA_VECTOR(rec);          // number of fish recaptured in each efficiency trial
  DATA_IVECTOR(efficIndex);  // index of catch data that each efficiency trial corresponds with 
  DATA_MATRIX(X);           // fixed effects design matrix for daily trap capture efficiency
  DATA_SPARSE_MATRIX(Z);    // random effects design matrix for daily trap capture efficiency
  DATA_STRUCT(terms_p,
              terms_t);     // Covariance structure for random effects
  
  
  //Parameters
  ////Migrants-model
  PARAMETER_VECTOR(mu_M);   // intercept of log-mean daily outmigrants
  PARAMETER(logit_phi_e);   // logit AR(1) coefficient for year-specific erros in log-mean daily outmigrants
  PARAMETER(ln_tau_d);      // log AR(1) process error precission for across year erros in log-mean daily emigrants
  PARAMETER(ln_tau_e);      // log AR(1) process error precission for year-specific errors in log-mean daily emmigrants 
  PARAMETER_VECTOR(delta);   // accross-year errors of log-means of daily emigrants (random effects)
  PARAMETER_MATRIX(epsilon); // year-specific errors of log-means of daily emigrants  (random effects)  
  PARAMETER(log_phi_NB);    // log -transformed overdispersion paramater for negative binomial 
  
  ////p-model
  PARAMETER_VECTOR(beta); // coefficients in trap capture efficiency model
  PARAMETER_VECTOR(b);      // random effects in trap capture efficiency model 
  PARAMETER_VECTOR(theta);  // random effect covariance parameters in trap capture efficiency model
  
  
  
  
  //Variables
  using namespace density;       // load "package" with AR1 distribtion
  vector<Type> jnll_comp(3);     // declare vector of 3 likelihood componenets (random effect probs, mark-recapture (efficiency trials) data, and catch data)
  jnll_comp.setZero();           // set likelihood components to zero
  
  
  Type phi_e=invlogit(logit_phi_e);// transform  correlation coefficient for yearly daily outmigrant errors to "0-1" space. I don't thing a negative correlation would ever be the case here, and must be >1 to be stationary
  Type sigma_d = 1/exp(ln_tau_d);  // transform log precision to standard deviation 
  Type sigma_e = 1/exp(ln_tau_e);  // transform log precision to standard deviation 
  Type phi_NB=exp(log_phi_NB);     // transform negative binomial "phi" parameter 
  
  
  ////latent emigrants abundance (M_hat)
  vector<Type> delta2(N_day);
  delta2 << Type(0),delta; //assume random walk for average daily errors starts at 0
  
  
  matrix<Type> M_hat(N_day,Nyears);                //declare matrix to hold expected emigrant abundances
  for( int Iyear=0; Iyear<Nyears; Iyear++){        //loop through years
    M_hat.col(Iyear)=
      vector<Type>(exp(
          Type(mu_M(Iyear)) +                      // annual mean 
            delta2 +                                // across-year daily error
            vector<Type> (epsilon.col(Iyear))));   // year-specific daily error
    
  }
  
  
  ////capture efficiency (p) on all days with catch  data
  vector<Type> logit_p  =  X * beta; // fixed effects
  logit_p+= Z * b ;                  // random effects
  
  
  
  //Probability of random effects
  
  //random walk for across-year daily erros
  ////first time step (assuming value at time 1 is 0)
  jnll_comp(0)-= dnorm(Type(delta(0)),Type(0),sigma_d,true);
  ////subsequent time steps
  for(int i =1; i<(N_day-1); i++){
    jnll_comp(0)-= dnorm(Type(delta(i)-delta(i-1)),Type(0),sigma_d,true);
  }
  
  
  ////AR1 likelihood for year-specific daily errors (representing daily deviations on top the average within a given year)
  for( int Iyear=0; Iyear<(Nyears); Iyear++){     //loop through years
    jnll_comp(0)+= SCALE(AR1(phi_e),
              sigma_e)(epsilon.col(Iyear)); 
  }
  
  ////prob of random effects in capture efficiency model
  jnll_comp(0) += allterms_nll(b, theta, terms_p, this->do_simulate);//phi
  
  
  //likelihood
  ////"trap capture efficiency" data
  vector<Type> sim_rec(rel.size());
  for (int I =0; I<rel.size(); I++){          //loop over all release-recapture experiments
    jnll_comp(1)-= dbinom_robust(rec(I),      
              rel(I),
              logit_p(efficIndex(I)),true);   
    
    SIMULATE { //simulate data for model fit checking
      sim_rec(I) = rbinom(rel(I),
              Type(invlogit(logit_p(efficIndex(I)))));
    }
  }
  REPORT(sim_rec);
  
  
  
  //// Likleihood of catch data 
  vector<Type> sim_catch(N_trap);
  
  for( int Iday=0; Iday<N_trap; Iday++){ // loop over all days with catch data
    Type expected_catch = Type(M_hat(catch_DOY(Iday),seriesFac(Iday))*
      invlogit(logit_p(Iday)));
    if(Use_NB){                          // if using negative binomial distribution
      jnll_comp(2) -= dnbinom2 (Catch(Iday), expected_catch,
                expected_catch+(pow(expected_catch,2)/phi_NB), true);
      
      SIMULATE { //simulate data for model fit checking
        sim_catch(Iday) = rnbinom2(expected_catch,
                  expected_catch+(pow(expected_catch,2)/phi_NB));
      }
      
    }else{                               // otherwise use Poisson observation likelihood
      jnll_comp(2) -= dpois(Catch(Iday),expected_catch, true);
    }
  }
  REPORT(sim_catch);
  
  //Objective function
  Type obj_fun = sum(jnll_comp);         //sum likelihood components
  
  //Reporting
  REPORT(M_hat);                // daily emigrants[t,y]
  return(obj_fun);
}
