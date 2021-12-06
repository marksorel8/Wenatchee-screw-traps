# 4 helper functions for fitting spawner-to-juveniile models


##-------------------------------------------------------------------------------------------------------------------
##function to create data input for models. returns a list of data inputs to feed to TMB model "Stock_recruit_LVM". 
make_dat_func<-function(streams=c( 0:2), #streams to model 0 = chiwawa, 1 = nason, 2 = white
                        LHs=c(1:4),     #life histories to model 1=fry, 2=sumemr, 3=fall, 4=smolt
                        n_f=1,           #number of latent variable factors
                        log_J_max_prior=log(15000)){
  
  # concatenate log mean of observed emigrant distributions,
  # log standard devations of observed emigrant distributions,
  # and observed spawners for all streams, life histories, and years
  log_m_means<-log_m_sds<-numeric(0) #empty vectors
  sl_i<-st_i<-t_i<-s_i<-l_i<-BY<-redd_yr<-integer(0)
  S_obs<-matrix(nrow=0,ncol=4)
  X<-matrix(nrow = 0,ncol=24)         #design matrix
  
  
  #construction of vectors
  for (i in streams){ #loop through streams
    for(j in LHs){      #loop though life histories
      
      years<-sort(unique(all_emigrants_estimates[[i*2+1]]$dat$years))#years (of observations)
      BY<-c(BY,(years-ifelse(j==4,2,1)))#brood years corresponding with observation years for given life history
      n_year<-length(years)#number of years
      
      # likelihood penalty mean and standard deviations for juvenile emigrant abundances
      if(j<4){
        means<-(all_emigrants_estimates[[i*2+1]]$boot$boot_log_means)[,j] #subyearlings
        sds<-(all_emigrants_estimates[[i*2+1]]$boot$boot_log_sds)[,j] #subyearlings
        
      }else{
        means<-all_emigrants_estimates[[i*2+2]]$boot$boot_log_means[,1] # yearlings
        sds<- all_emigrants_estimates[[i*2+2]]$boot$boot_log_sds[,1] # yearlings
      }
      # log means- concatenate current vector with all life histories
      log_m_means<-c(log_m_means,                                 #current vector
                     means)   
      # log standard deviations- concatenate current vector with all life histories
      log_m_sds<-c(log_m_sds,                                   #current vector
                   sds)   
      
      # observed redds
      redds_sub<-filter(redds,stream== c("Chiwawa","Nason","White")[(i+1)])
      redd_yr<-c(redd_yr,redds_sub$redds[match((years-ifelse(j==4,2,1)),redds_sub$Year)]/trib_lengths[(i+1)])
      
      # design matrix of environmental covariates
      X_s<-matrix(NA,nrow=n_year,ncol=24)
      #$ winter discharge brood year (incubation)
      X_s[(1:n_year),j]<-
        scale(flow_covs$winter_high[match((years-ifelse(j==4,2,1)), #Year of incubation
                                          flow_covs$winter_high$Year),(i+2)])

      #summer discharge brood year +1 (summer rearing) Z-scored
      if(j>1){
        X_s[(1:n_year),(j+3)] <-
          scale(flow_covs$summer_low[match((years-ifelse(j==4,1,0)), #Year of incubation
                                           flow_covs$summer_low$Year),(i+2)])
      }

      #winter discharge brood year +1 (overwintering) Z-scored
      if(j==4){
        X_s[(1:n_year),8] <- 
          scale(flow_covs$winter_high[match((years-1), #Year of overwintering as parr 
                                            flow_covs$winter_high$Year),(i+2)])
      }
      # add rows to design matrix
      X<-rbind(X,X_s)
      
      
      # stream by life history index 
      sl_i<- c(sl_i,
               rep(i*4+j,n_year))
      
      #stream by year index
      st_i<-c(st_i,
              (ifelse(j==4,1,2):(n_year+ifelse(j==4,0,1)))+(i*100))
      
      #brood year index
      t_i<-c(t_i,years-ifelse(j==4,2,1))
      
      #life history index (for reference/plotting only)
      l_i<- c(l_i,
              rep(j-1,n_year))
      
      #stream index (for reference/plotting only)
      s_i <- c(s_i,
               rep(i,n_year))
      
    }
  }
  
  

  
  #get rid of columns not used in envornmental cpovariate design matrix
  X<-X[,which(apply(X,2,function(x)sum(!is.na(x)))>0),drop=FALSE]
  X[is.na(X)]<-0
  
  #construct data list
  dat<-list(n_sl = length(unique(sl_i)),           #n umber of distinct life history by stream combinations
            n_i = length(sl_i),                    # number of years x streams x life histories
            n_t = length(unique(t_i)),             # number of years
            n_l = length(LHs),             # number of life histories
            n_s = length(streams),             # number of life histories
            n_f = n_f,                             # number of latent variable factors
            J_obs = log_m_means,                   # Observed log recruits (posterior mean from juvenile model)
            J_obs_sd = log_m_sds,                  # Observed log recruits standard error (posterior sd from juvenile model)
            log_S_obs=log((redd_yr[!duplicated(st_i)] 
                           [order(st_i[!duplicated(st_i)])])), # log observed redds
            S_obs_CV = 0.1,                       # observed redd CV
            sl_i = factor(sl_i),                   # stream by life history index for observations
            st_i = factor(st_i),                   # stream by year index for observations
            t_i = factor(t_i),                     # year index of each observation
            X = as(X,"dgTMatrix"),                 # sparse design matrix of environmental covariates 
            s_i = s_i,                             # stream index for reference (not used in model)
            l_i = l_i,                             # LHP index for reference (not used in model)
            BY=BY,                                  # brood year of each observation, for reference (t_i is index used in model)
            log_J_max_prior=log_J_max_prior
  )
            return(dat)
}   

##-------------------------------------------------------------------------------------------------------------------

#function to construct parameters list. returns a list of initial values for paramaters to feed to TMB model "Stock_recruit_LVM.
make_params_func<-function(dat,rate){
  
  params<-list( 
                #fixed effects
                 beta_alpha=rnorm(dat$n_l,log(10),1),   # log alpha interceptrs (by life-history)
              
                log_sigma_alpha=rep(0,dat$n_l),                   # log alpha random effect SD
                log_sigma_Jmax=rep(0,dat$n_l),                   # log Jmax random effect SD
                log_sigma_gamma=rep(0,dat$n_l),                  # log gamma random effect SD
                log_sigma_Bgamma=rep(0,dat$n_l),                 # regularizing penalty standard deviation
                log_sigma_BJmax=rep(0,dat$n_l),                  # regularizing penalty standard deviation
                
                beta_e=rep(0,ncol(dat$X)),                     # environmental covariate coefficients for process error
                Loadings_vec=rnorm(dat$n_f*(dat$n_sl)-dat$n_f*(dat$n_f-1)/2,0,.1), # latent variable factor loadings
                log_sigma_eta=log(abs(rnorm(dat$n_sl,-1,.1))),# idiosyncratic process error log SDs
                rate=rnorm(4,log(c(1,.5,.5,1)),.5),                                     #penalty rates of regularizing penalties/priors
                
                #random effects
                beta_Jmax= rnorm(dat$n_l,log(15000),1), # log asymptotic maximum recruitment (Jmax) coefficients (by life-history) intercept
                beta_gamma= rnorm(dat$n_l,0,.2),        # log gamma coefficients (by life history)
                
                eps_alpha=rnorm(dat$n_sl,0,.1),                         # log alpha random effects
                eps_Jmax=rnorm(dat$n_sl,0,.1),                       # log Jmax random effects
                eps_gamma= rnorm(dat$n_sl,0,.1),                        # log gamma random effects
                log_S_hat = rnorm(length(dat$log_S_obs), dat$log_S_obs,.1), # latent spawners random effects
                Omega_xf=matrix(0,dat$n_t,dat$n_f),            #  latenct variable factors
                eta=numeric(dat$n_i)                         #  idisyncratic process error random effects

                )                       
  
  return(params)
}

##-------------------------------------------------------------------------------------------------------------------


#function to make map (controls which parameters are fit and which are held fixed at initials).
make_map<-function(mod_dat,fit_env=TRUE){

  map<-list()
  #fix environmental covaraites at zero if fit_env= FALSE
  if(!fit_env){
    map$beta_e<-factor(rep(NA,ncol(mod_dat$X))) # fix environmental covariates at 0.0
  }
  
  
  if(mod_dat$n_f>mod_dat$n_sl){stop("number of factors is > number of stream x LH combos")}
  return(map)
}

##-------------------------------------------------------------------------------------------------------------------

## function to attempt model fitting a specified number of times. Returns list of bets fit model (based on BIC), model and "fit", objects, and a vector of the BICs for all fitting attempts. Plus a "report" object and input data.
fit_mod_func<-function(streams, LHs, n_f, fit_env, fit_attempts,log_J_max_prior){
  dat<-make_dat_func(streams,LHs,n_f,log_J_max_prior) # make model data
  mod_map<-make_map(mod_dat=dat,fit_env = fit_env)                 # create map
  fit<-mod<-report<-NA # placeholders for converged model and "fit" objects
  BIC_vec<-rep(Inf,fit_attempts)
  for ( i in 1:fit_attempts){
    params<-make_params_func(dat,rate)        # make initial parameter values
        mod_i<-TMB::MakeADFun(dat,params,random=c("log_S_hat","Omega_xf","eps_alpha","eps_gamma","eps_Jmax","beta_gamma","beta_Jmax","eta"),DLL="Stock_recruit_LVM",map=mod_map,silent = TRUE)
    fit_i<-NA # clear previous fit object 
    try(fit_i<-TMBhelper::fit_tmb(mod_i, newtonsteps = 1,getJointPrecision = TRUE)) # optimize
    BIC_mod<- NA #clear previous BIC
    try(BIC_mod<-(2*fit_i$objective+log(dat$n_t)*fit_i$number_of_coefficients[2])) # calculate  BIC
    try(if(min(BIC_vec,na.rm = TRUE)>BIC_mod){ # save model and "fit" object if BIC is lower than previous
      fit<-fit_i
      mod<-mod_i
      try(report<-mod$report())
    })
    try(BIC_vec[i]<-BIC_mod) # save BIC
  }
  return(list(mod=mod,fit=fit,BIC_vec=BIC_vec,report=report,dat=dat))
}

#function to try a few more iterations if didn't converge in initial tries
fit_mod_iter<-function( streams, LHs, n_f, fit_env, fit_attempts, additional_attempts,log_J_max_prior){
  x<-fit_mod_func( streams, LHs, n_f, fit_env, fit_attempts,log_J_max_prior)
    i<-0
    while(i<=additional_attempts & (min(x$BIC_vec,na.rm=TRUE)==Inf)){
      x<-fit_mod_func( streams, LHs, n_f, fit_env, fit_attempts=1,log_J_max_prior)
      i<-i+1
    }
    
    return(x)
}



#function for simulation from a multivariate normal, taken fron J. Thorson's Fish Utils. Package
rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
  set.seed( random_seed )
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L = Matrix::Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  return(mu + z)
}

