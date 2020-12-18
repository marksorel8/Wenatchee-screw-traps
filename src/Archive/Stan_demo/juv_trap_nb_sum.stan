data {
  int<lower=1> N_MR;                  // number of weekly mark-recapture (efficiency) observations
  int<lower=1> MR_week[N_MR];         // week of each mark-recapture observation
  int<lower=1> mark[N_MR];            // number of marked fish released
  int<lower=0> recap[N_MR];           // number of marked fish recaptured
  int<lower=1> N_trap;                // number of daily trap catch observations
  int<lower=1> trap_day[N_trap];      // day of each catch observation
  int<lower=1> trap_week[N_trap];     // week of each catch observation
  int<lower=1> NX_p;                  // number of efficiency covariates
  matrix[max(trap_week),NX_p] X_p;    // design matrix of efficiency covariates (first column is 1) 
  int<lower=1> NX_M;                  // number of covariates for daily outmigrants
  matrix[max(trap_day),NX_M] X_M;     // design matrix of outmigrant covariates (first column is 1)
  int<lower=0> C[N_trap];             // daily trap catch observations
  int<lower=1> elapsed_time[N_trap];  // time (days) of sampling for each trap catch obs
  int<lower=0, upper=1> Use_NB;                //Flag for whether to use the the Negative Binomial distribution in the observation model
}

transformed data {
  int<lower=1> max_week;  // maximum trapping week
  int<lower=1> max_day;   // maximum trapping day
  
  max_week = max(trap_week);
  max_day = max(trap_day);
}

parameters {
  vector[NX_p] beta_p;               // regression coefs for capture probability (first is intercept)
  real<lower=0> sigma_p;             // hyper-SD logit capture probability
  vector[max_week] logit_p_z;        // weekly logit capture probability (z-scores)
  vector[NX_M] beta_M;               // regression coefs for log-mean daily outmigrants (first is intercept)
  real<lower=-1,upper=1> phi_M;      // AR(1) coefficient for log-mean daily outmigrants
  real<lower=0> sigma_M;             // AR(1) process error SD for log-mean daily outmigrants
  vector[max_day] log_M_hat_z;       // log-means of daily outmigrants (z-scores)
  real<lower=0> inv_sqrt_phi_obs[Use_NB];     // inverse square root dispersion paramater for the Negative Binomial observation model for catch. 
}

transformed parameters {
  vector<lower=0,upper=1>[max_week] p;   // weekly capture probability
  vector[max_day] mu_M;                  // intercept of AR(1) process for log-mean daily outmigrants
  vector[max_day] log_M_hat;             // log-means of daily outmigrants
  vector<lower=0>[max_day] M_hat;        // means of daily outmigrants
  vector[N_trap] C_hat;                  // expected catches
  real<lower=0> phi_obs[Use_NB];         // dispersion paramater for the Negative Binomial observation model for catch. 
  
  if(Use_NB)
  phi_obs[1]= (1/inv_sqrt_phi_obs[1])^2;

  
  // Hierarchical noncentering for weekly capture probability
  // (aka "Matt trick"; see Stan 2.17.0 manual Ch. 28.6)
  p = inv_logit(X_p*beta_p + sigma_p*logit_p_z);
  
  // Hierarchical noncentering of AR(1) process for log-means of daily outmigrants
  // Prior on initial state is the stationary distribution
  mu_M = X_M * beta_M;
  log_M_hat[1] = mu_M[1] + (sigma_M/sqrt(1 - phi_M^2))*log_M_hat_z[1];
  for(t in 2:max_day)
    log_M_hat[t] = mu_M[t] + phi_M*(log_M_hat[t-1] - mu_M[t-1]) + sigma_M*log_M_hat_z[t];
  M_hat = exp(log_M_hat);
  
  // Expected catches
  // Note that the Poisson distribution is closed under addition
  C_hat = M_hat .* p[trap_week];
}

model {
  //----------------
  // Priors
  //----------------
  
  // log Jacobian of logit transform for capture probability intercept
  // implies mean(p) ~ Unif(0,1) given all covariates are at their sample means
  target += log_inv_logit(beta_p[1]) + log1m_inv_logit(beta_p[1]);
  if(NX_p > 1)
    beta_p[2:NX_p] ~ normal(0,5);
  sigma_p ~ normal(0,5);        // implicitly truncated to [0,Inf)
  logit_p_z ~ normal(0,1);      // implies logit(p) ~ N(logit(mu_p), sigma_p)
  beta_M ~ normal(0,5); 
  // phi_M ~ uniform(-1,1) implicit
  sigma_M ~ normal(0,10);       // implicitly truncated to [0,Inf)
  log_M_hat_z ~ normal(0,1);    // implies log(M_hat[t]) ~ AR1(mu_M, phi_M, sigma_M)
  
  if(Use_NB)
   inv_sqrt_phi_obs ~ normal(0,1);
  //----------------
  // Likelihood
  //----------------
  
  // Mark-recapture observations
  recap ~ binomial(mark, p[MR_week]); 
  
  // Trap catch observations
 if(Use_NB){
   C ~ neg_binomial_2(C_hat,phi_obs[1]);
   }
 else{
  // Note that a Poisson RV thinned by binomial sampling is Poisson
 C ~ poisson(C_hat);}
}

generated quantities {
  int M[max_day];  // daily outmigrants
  int M_tot;       // total outmigrants

 for(t in 1:max_day){

   if(Use_NB){
   M[t] = neg_binomial_2_rng(M_hat[t],phi_obs[1]);}
   else{
   M[t] = poisson_log_rng(log_M_hat[t]);}
 }

  M_tot = sum(M);
}
