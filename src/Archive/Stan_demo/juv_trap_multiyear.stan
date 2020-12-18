data {
  int<lower=1> N_MR;                  // number of weekly mark-recapture (efficiency) observations
  int<lower=1> MR_year[N_MR];         // brood year of each mark-recapture observation
  int<lower=1> MR_week[N_MR];         // week of each mark-recapture observation
  int<lower=1> mark[N_MR];            // number of marked fish released
  int<lower=0> recap[N_MR];           // number of marked fish recaptured
  int<lower=1> N_trap;                // number of daily trap catch observations
  int<lower=1> trap_year[N_trap];     // brood year of each catch observation
  int<lower=1> trap_week[N_trap];     // week of each catch observation
  int<lower=1> trap_day[N_trap];      // day of each catch observation
  int<lower=1> NX_p;                  // number of efficiency covariates
  matrix[max(trap_week)*max(trap_year),NX_p] X_p; // design matrix of efficiency covariates (first column is 1)
  // int<lower=1> NX_M;                  // number of covariates for daily outmigrants
  // matrix[max(trap_day)*max(trap_year),NX_M] X_M; // design matrix of outmigrant covariates (first column is 1)
  int<lower=0> C[N_trap];             // daily trap catch observations
  int<lower=1> elapsed_time[N_trap];  // time (days) of sampling for each trap catch obs
}

transformed data {
  int<lower=1> N_year;  // maximum brood year
  int<lower=1> N_week;  // maximum trapping week
  int<lower=1> N_day;   // maximum trapping day
  
  N_year = max(trap_year);
  N_week = max(trap_week);
  N_day = max(trap_day);
}

parameters {
  vector[NX_p] beta_p;               // regression coefs for capture probability (first is intercept)
  real<lower=0> sigma_p_year;        // interannual hyper-SD of logit capture probability
  real<lower=0> sigma_p_week;        // weekly hyper-SD of logit capture probability
  vector[N_year] logit_p_year_z;     // logit annual capture probability (z-scores)
  row_vector[N_week] logit_p_week_z; // logit weekly capture probability (z-scores)
  // vector[NX_M] beta_M;               // regression coefs for log-mean daily outmigrants (first is intercept)
  vector[N_year] mu_M;               // annual means of log-mean outmigrants
  vector<lower=-1,upper=1>[N_year] phi_M; // diagonal MAR(1) coefs of daily log-mean outmigrants
  vector<lower=0>[N_year] sigma_M;   // annual process error SDs of log-mean outmigrants
  cholesky_factor_corr[N_year] L_M;  // Cholesky factor of correlation matrix for log-mean outmigrants
  matrix[N_year,N_day] log_M_hat_z;  // annual log-means of daily outmigrants (z-scores)
}

transformed parameters {
  matrix<lower=0,upper=1>[N_year,N_week] p; // capture probability
  // vector[N_day] mu_M;                      // intercept of MAR(1) for log-mean daily outmigrants
  matrix[N_year,N_day] log_M_hat;           // annual log-means of daily outmigrants
  matrix<lower=0>[N_year,N_day] M_hat;      // annual medians of daily outmigrants
  vector<lower=0>[N_trap] M_hat_cumsum;     // daily means summed over days that trap is fishing
  vector<lower=0>[N_trap] C_hat;            // expected catches
  
  // Hierarchical noncentering for weekly capture probability
  // (aka "Matt trick"; see Stan 2.17.0 manual Ch. 28.6)
  p = inv_logit(to_matrix(X_p*beta_p, N_year, N_week) + 
                  rep_matrix(sigma_p_year*logit_p_year_z, N_week) + 
                  rep_matrix(sigma_p_week*logit_p_week_z, N_year));
  
  // Hierarchical noncentering of MAR(1) process for log-means of daily outmigrants
  // (Multivariate Matt trick)
  // Prior on each year's initial state is the stationary distribution
  // (NOTE: This is not correct as written; MAR(1) stationary distn is needed)
  log_M_hat[,1] = mu_M + (sigma_M ./ sqrt(1 - phi_M .* phi_M)) .* log_M_hat_z[,1];
  M_hat[,1] = exp(log_M_hat[,1]);
  for(t in 2:N_day)
  {
    log_M_hat[,t] = mu_M + phi_M .* (log_M_hat[,t-1] - mu_M) + sigma_M .* (L_M * log_M_hat_z[,t]);
    M_hat[,t] = exp(log_M_hat[,t]);
  }
  
  // Expected catches
  // Note that the Poisson distribution is closed under addition
  {
    matrix[N_day,N_year] M_hat_t = M_hat'; // local variable; transposed for efficiency
    for(i in 1:N_trap)
    {
      M_hat_cumsum[i] = sum(M_hat_t[(trap_day[i] - elapsed_time[i] + 1):trap_day[i], trap_year[i]]);
      C_hat[i] = M_hat_cumsum[i] * p[trap_year[i],trap_week[i]];
    }
  }
}

model {
  vector[N_MR] p_MR;  // capture probabilities corresponding to MR observations
  
  //----------------
  // Priors
  //----------------
  
  // log Jacobian of logit transform for capture probability intercept
  // implies mean(p) ~ Unif(0,1) given all covariates are at their sample means
  target += log_inv_logit(beta_p[1]) + log1m_inv_logit(beta_p[1]);
  if(NX_p > 1)
    beta_p[2:NX_p] ~ normal(0,3);
  sigma_p_year ~ normal(0,5);   // implicitly truncated to [0,Inf)
  sigma_p_week ~ normal(0,5);   // implicitly truncated to [0,Inf)
  // Next two lines imply logit(p) = mu_p + logit_p_year + logit_p_week,
  // logit_p_year ~ N(0,sigma_p_year), logit_p_week ~ N(0,sigma_p_week)
  logit_p_year_z ~ normal(0,1); 
  logit_p_week_z ~ normal(0,1); 
  // beta_M ~ normal(0,5); 
  mu_M ~ normal(0,5); 
  // phi_M ~ uniform(-1,1) implicit
  sigma_M ~ normal(0,10);       // implicitly truncated to [0,Inf)
  L_M ~ lkj_corr_cholesky(1);   // see Stan manual for LKJ distribution
  to_vector(log_M_hat_z) ~ normal(0,1);  // log(M_hat[,t]) ~ MAR1(mu_M, diag(phi_M), diag(sigma_M.^2))
  
  //----------------
  // Likelihood
  //----------------
  
  // Mark-recapture observations
  for(i in 1:N_MR)
    p_MR[i] = p[MR_year[i],MR_week[i]];
  recap ~ binomial(mark, p_MR); 
  
  // Trap catch observations
  // Note that a Poisson RV thinned by binomial sampling is Poisson
  C ~ poisson(C_hat);
}

generated quantities {
  matrix[N_year,N_year] Q_M;  // MAR(1) daily process error covariance matrix
  matrix[N_year,N_day] M;     // daily outmigrants
  vector[N_year] M_tot;       // total outmigrants
  
  Q_M = multiply_lower_tri_self_transpose(L_M);
  for(t in 1:N_day)
    for(y in 1:N_year)
      M[y,t] = poisson_rng(M_hat[y,t]);
  M_tot = M * rep_vector(1,N_day);  // sum across days (columns) to get annual totals
}
