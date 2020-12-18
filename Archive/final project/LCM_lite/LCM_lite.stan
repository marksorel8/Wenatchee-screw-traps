

data {
  vector[23] wild_spawners;         //wild spawner data
  vector[23] hatchery_spawners;     //HOR spawner data
  //vector[N] wild_broodstock;      // wild broodstok removals
  vector[4] em_obs_log_means[22];
  vector[4] em_obs_log_sds[22];
  
  //matrix[22,2000] emigrants[4];    // array of juvenile emigrant abundance estimates
}


parameters {
  
  vector[23] log_S_hat; 
  real <lower=0> S_hat_CV;         //spawner observation error
  
  vector[4] log_em_hat[22];
  vector<lower=0>[4] em_hat_CV[22];   //emigrant observation error
 
  vector<lower=0>[4] alpha;
  vector<lower=0>[4] R_max;
  vector<lower=0>[3] d;

  cov_matrix[4] em_proc_cov;              // MVN process error for juvenile emigrants
}

transformed parameters {
 
 vector <lower=0>[4] em_pred[22];
 
 for ( i in 1:4){
    if(i<4){
      for (j in 2:23){
      em_pred[(j-1),i]= R_max[i]*(1-exp(-((exp(log_S_hat[j])/ alpha[i] )^d[i])));
    //Weibull CDF
      }
}else
   for (j in 2:23)
   em_pred[(j-1),i]= (alpha[i]*exp(log_S_hat[j-1]))/(1+alpha[i]*exp(log_S_hat[j-1])/R_max[i]); //Bev-holt
 }
 
 }



model {
 //priors
 S_hat_CV ~ normal(.05,.01);
log_S_hat~ normal(0,5);

  for ( i in 1:22){
em_hat_CV[i] ~ normal(0,5);
}
for ( i in 1:4){
em_proc_cov[i] ~ normal(0,5);
}

alpha ~ normal(500,500);
R_max ~ normal(5000,1000);
 d ~ normal(0,5);



  for ( i in 1:22){
    for ( j in 1:4){
 log_em_hat[i,j] ~ normal(em_obs_log_means[i,j],em_obs_log_sds[i,j]);   // emigrant observatin likelihood
    }
}
 
 //likelihood 
 for ( i in 1:22){
   
 log_em_hat[i] ~ multi_normal(log(em_pred[i]),em_proc_cov);   // emigrant process likelihood
}
// for ( i in 1:22){                                  //emigrant observation likelihood
  // for (j in 1:4){
  // log(emigrants[j,i]) ~  normal(log_em_hat[i,j],em_hat_CV[i,j]);
  // }
 //}

log(wild_spawners+hatchery_spawners) ~ normal(log_S_hat,S_hat_CV);
}

