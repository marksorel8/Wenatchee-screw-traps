#The purpose of this script is to fit models of transition between spawners and juvenile emigrants expressing different juvenle life histories

#libraries
library(here)
library(tidyverse)
library(viridis)
library(TMB)
library(TMBhelper)
library(readxl)


##   Load data ##


#read and combine spawner observation data from Hatchery program annual report (Johnson et al 2019)
spawners<-read.csv(here("data","Chiwawa","ChiwEscapement.csv")) %>% mutate(stream="Chiwawa") %>% 
  bind_rows(read.csv(here("data","Nason and White","Nas_Escapement.csv")) %>% mutate(stream="Nason")) %>% 
  bind_rows(read.csv(here("data","Nason and White","White_Escapement.csv")) %>% mutate(stream="White")) %>%
  rowwise() %>%  mutate(total_spawners=sum(NOS,HOS,na.rm=T))

#read egg data
Eggs <- read_csv(here("data","Wenatchee spring Chinook egg retention (12-12-19)_Sorel_mod.csv"))

#read emigrant abundance posteriors
load(here("results","Rdata","emigrant_estimates Oct 01 2020.Rdata"))

#load stream flow covariate values(flow_covs)
source(here("src","covariates.r"))

# length of habitat surveyed for spawners in each tributary (based on Montoring annual report)
trib_lengths<-c(32.5,15.4,16.1)
names(trib_lengths)<-c("Chiwawa","Nason","White")


#snowtel data
snow<-read_excel(here("data","snowpack_Apr_01.xlsx"),sheet="SNOTEL")


##   Load functions ##

#function to create data input for models
make_dat_func<-function(streams=c( 0:2), #streams to model 0 = chiwawa, 1 = nason, 2 = white
                        LHs=c(1:4)){    #life histories to model 1=fry, 2=sumemr, 3=fall, 4=smolt

  # concatenate log mean of observed emigrant distributions,
# log standard devations of observed emigrant distributions,
# and observed spawners for all streams, life histories, and years
log_m_means<-log_m_sds<-numeric(0) #empty vectors
sl_i<-st_i<-lt_i<-l_i<-t_i<-s_i<-BY<-integer(0)
S_obs<-matrix(nrow=0,ncol=4)
X<-matrix(nrow = 0,ncol=24)         #design matrix


#construction of vectors
for (i in streams){ #loop through streams
for(j in LHs){      #loop though life histories

years<-as.numeric(rownames(all_emigrants_estimates[[i*2+1]]$LH_sums))#years (of observations)
BY<-c(BY,(years-ifelse(j==4,2,1)))#brood years corresponding with observation years for given life history
n_year<-length(years)#number of years

# likelihood penalty mean and standard deviations for juvenile emigrant abundances
if(j<4){
  means<-all_emigrants_estimates[[i*2+1]]$LH_sums[,j] #subyearlings
  sds<-all_emigrants_estimates[[i*2+1]]$LH_sums_sd[,j] #subyearlings

}else{
  means<-all_emigrants_estimates[[i*2+2]]$LH_sums # yearlings
  sds<- all_emigrants_estimates[[i*2+2]]$LH_sums_sd # yearlings
}
# log means- concatenate current vector with all life histories
log_m_means<-c(log_m_means,                                 #current vector
              means)   
# log standard deviations- concatenate current vector with all life histories
log_m_sds<-c(log_m_sds,                                   #current vector
             sds)   
               

# observed eggs
egg_mat<-select(filter(Eggs,
                Stream== c("Chiwawa","Nason","White")[(i+1)]&
                  Year %in%
                  as.character(years-ifelse(j==4,2,1))),Wild:PHOS)


egg_mat[,1:3]<-egg_mat[1:3]/trib_lengths[(i+1)] # scale by spawning stream length

S_obs<-rbind(S_obs,egg_mat) #bidn with previous streams and life histories


# design matrix
X_s<-matrix(NA,nrow=n_year,ncol=24)
#winter discharge brood year (incubation)
X_s[(1:n_year),j]<-#(j+1)] <-##j+i*4] <-
  scale(flow_covs$winter_high[match((years-ifelse(j==4,2,1)), #Year of incubation 
                                                             flow_covs$winter_high$Year),(i+2)]) 

#summer discharge brood year +1 (summer rearing)
if(j>1){
  X_s[(1:n_year),(j+3)] <-##(j+3)]<-# (j+i*3)+11] <-
  scale(flow_covs$summer_low[match((years-ifelse(j==4,1,0)), #Year of incubation 
                                                            flow_covs$summer_low$Year),(i+2)]) #subyearling emigrants
}

#winter discharge brood year +1 (overwintering)
if(j==4){
  X_s[(1:n_year),8] <- ##(i)+22] <-
  scale(flow_covs$winter_high[match((years-1), #Year of overwintering as parr 
                                                                  flow_covs$winter_high$Year),(i+2)]) #yearling emigrants
}

# X_s[(1:n_year),8+j]<-#(j+1)] <-##j+i*4] <-
#   scale(snow$`Stevens Pass`[match((years-ifelse(j==4,2,1)), #Year of incubation
#                                     snow$year)])

#X_s[,9]<-egg_mat$PHOS


X<-rbind(X,X_s)

beta_e_i<-c(rep(LHs,length(streams)),rep(LHs[LHs!=1],length(streams))+3,rep(LHs[LHs=4],length(streams))+4)

# stream by life history index 
sl_i<- c(sl_i,
         rep(i*4+j,n_year))

#stream by year index
st_i<-c(st_i,
       (ifelse(j==4,1,2):(n_year+ifelse(j==4,0,1)))+(i*100))

# life history by year index
lt_i<-c(lt_i,
        years+(j*100))

 #life history index 
l_i<- c(l_i,
         rep(j-1,n_year))

#stream index 
s_i <- c(s_i,
         rep(i,n_year))
#brood years index
t_i<-c(t_i,years-ifelse(j==4,2,1))
}
}
#get rid of columns not used
X<-X[,which(apply(X,2,function(x)sum(!is.na(x)))>0),drop=FALSE]
#scale
X<-scale(X)
#fill in NAs with zeros
X[is.na(X)]<-0
#check scaling
apply(X,2,sd,na.rm=T)
apply(X,2,mean,na.rm=T)



#read capacity estimates
#capacity_estimates<-read_excel(here("data","Wenatchee_capacity_estimates.xlsx"),sheet=1,skip=2)[,1:2] %>% filter(Tributary%in%c("Chiwawa River","Nason Creek","White River"))

#construct data list
dat<-list(n_sl = length(unique(sl_i)),           #number of distinct life history by stream combinations
          n_i = length(sl_i),      # number of years x streams x life histories
          n_t = length(unique(t_i)), # number of years
          n_l = length(unique(l_i)),                   # number of life histories
          n_s = length(unique(s_i)),
          n_st = length(unique(st_i)), # number of stream by year combinations
          n_beta_e = length(beta_e_i), #number of covariate coefficients
          mod= rep(1,length(unique(sl_i))),                      # functional form of process model
          R_obs = log_m_means,    # Observed log recruits (posterior mean from juvenile model)
          R_obs_sd = log_m_sds, # Observed log recruits standard error (posterior sd from juvenile model)
          log_S_obs=log((S_obs$Total[!duplicated(st_i)][order(st_i[!duplicated(st_i)])]/4600)),
          S_obs_CV = 0.05,
          l_i= factor(l_i),         # life history index
          sl_i = factor(sl_i),  # stream by life history index
          st_i = factor(st_i),  # stream by year index for observations
          lt_i = factor(lt_i),  # life history by year index for observations
          s_theta_i=factor(rep(1:length(unique(s_i)),times=unlist(lapply(tapply(t_i,s_i,unique),length)))), #stream index of each stream*year random effect
          beta_e_i=factor(beta_e_i), # index of random covariate coefficients
          t_i=factor(t_i),#stream/ life history index of each observation
          s_i=factor(s_i),       # stream intex of each observation
          X=as(X,"dgTMatrix"),#design matrix for annual errors
          BY=BY,            #brood year for reference, not used in model
          S_obs=S_obs,      # egg info for reference
          num=100)
return(dat)
}     
    
#function to construct parameters list
make_params_func<-function(dat){

params<-list(beta_a=rep(0,dat$n_l),            #intrinsic productivity intercepts (by life-history)
             beta_rm=rep(0,dat$n_l), #Asymptotic maximum recruitment coefficient and intercept
             beta_f=rep(0,dat$n_l),
             beta_e=rep(0,ncol(dat$X)),
             #hyper_beta_e=rep(0,length(levels(dat$beta_e_i))), #hyper means of covariate effects across streams
             #zeta_a=rep(0,3),
            # eta_a=rep(0,4),
            # zeta_rm=rep(0,3),
            # eta_rm=rep(0,4),
            # zeta_f=rep(0,3),
            # eta_f=rep(0,4),
             #log_sigma_a=rnorm(dat$n_l,-2,1),              #standard deviation of stream productivity errors (relative to life-stage expectation)
             #log_sigma_rm=rnorm(dat$n_l,-2,1),             #standard deviation of stream max recruitment errors (relative to life-stage expectation)
             #log_sigma_e=log(rnorm(length(levels(dat$beta_e_i)),.5,.05)), #standard deviations of stream by covariate coefficients
            #log_sigma_f=rnorm(dat$n_l,-2,1),
             #log_sigma_e=rnorm(ncol(dat$X),-2,1),
             #log_sigma_zeta_a=0,
             #log_sigma_eta_a=0,
             #log_sigma_zeta_rm=0,
             #log_sigma_eta_rm=0,
             #log_sigma_zeta_f=0,
             #log_sigma_eta_f=0,
             #par_theta = rep(0,3),
             log_proc_sigma=(rnorm(dat$n_sl,.1,.015)),    #process error standard deviations
            #log_proc_sigma=2,    #process error standard deviations
            #log_sigma_delta = rnorm(1,1,.15),            #not used standard deviation of global erros
            log_sigma_theta = (rnorm(dat$n_s,.1,.015)),    # sd of stream specific errors
            #log_sigma_kappa = rnorm(dat$n_l,-2,1),
            #log_sigma_log_fifty = rnorm(dat$n_sl,-1,.2),
            #log_sigma_log_rm = rnorm(dat$n_sl,-1,.2),
             #proc_theta=rep(0,(3*3-3)/2),                   #process errors correlations
             #random effects
            log_S_hat = rnorm(dat$n_st, dat$log_S_obs,.1),        #log spawner abundance
             eps_a= rnorm(dat$n_sl,log(50),.5),                  #productivity stream specific
             eps_rm=rnorm(dat$n_sl,log(500),.5),                  #max recruitment  stream specific
             eps_f=rnorm(dat$n_sl,log(.5),.1),                     # depensation stream specific
            #kappa=rep(0,length(unique(lt_i))),  # random life history by year error
            theta=rep(0,dat$n_st),              # random stream by year error
            delta=rep(0,dat$n_t),               # random year errors
#             gamma=matrix(rnorm(12*22),nrow=dat$n_t,ncol=dat$n_sl))             #process errors in juvenile recruitment 
gamma=numeric(dat$n_i))                        # random life history by year error 
return(params)
}


 #function to make map for models
 make_map<-function(dat){
   
   vec1<-factor(seq(1,dat$n_sl))
   vec2<-factor(seq(1,dat$n_sl))
   #vec1.1<-vec2.2<-1:4
   for ( i in 1:dat$n_sl){
     if(dat$mod[i]%in%c(5,9:10)){ vec1[i]<-factor(NA) #eps_rm map
     #vec1.1[i]<-factor(NA)
     }
     if(dat$mod[i]%in%(6:9)){ vec2[i]<-factor(NA)  #eps_f map
     #vec2.2[i]<-factor(NA)
     }
   }   
   map=list(beta_f=factor(rep(NA,dat$n_l)),#factor(vec2.2),
            beta_rm=factor(rep(NA,dat$n_l)),#factor(vec1.1),
            beta_a =factor(rep(NA,dat$n_l)),
            eps_f=factor(vec2),
            #log_sigma_a=factor(rep(NA,4)),            
            eps_rm=factor(vec1))
           # log_sigma_rm=factor(rep(NA,4)),#factor(vec1.1),
            #log_sigma_f=factor(rep(NA,4)))#factor(vec2.2)))
   #,
            #log_sigma_log_fifty = factor(vec2))
   
if(dat$n_s<2){
  map$log_sigma_theta<-factor(NA)
  map$theta<-factor(rep(NA,length(params$theta))) 
}
   return(map)
 }
 

 
 
#funciton to make bounds for model
 make_bounds_func<-function(){
   L=c(beta_e =rep(-Inf,length(params$beta_e)-sum(is.na(map$beta_e))),
       #hyper_beta_e=rep(-Inf,length(mod$env$parameters$hyper_beta_e)),
       #log_sigma_e=rep(.0001,length(mod$env$parameters$hyper_beta_e)),
       log_proc_sigma=rep(.0001,length(params$log_proc_sigma)),
       #log_sigma_delta=.0001,
       log_sigma_theta=rep(.0001,length(params$log_sigma_theta)),
       eps_a=rep(-Inf,length(params$eps_a)),
       eps_rm=rep(-Inf,length(params$eps_rm)-sum(is.na(map$eps_rm))),
       eps_f=rep(-15,length(params$eps_a)-sum(is.na(map$eps_f))))
   
   U=c(beta_e =rep(Inf,length(params$beta_e)-sum(is.na(map$beta_e))),
       #hyper_beta_e=rep(Inf,length(mod$env$parameters$hyper_beta_e)),  
       #log_sigma_e=rep(Inf,length(mod$env$parameters$hyper_beta_e)),
       log_proc_sigma=rep(Inf,length(params$log_proc_sigma)),
     #log_sigma_delta=Inf,
       log_sigma_theta=rep(Inf,length(params$log_sigma_theta)),
       eps_a=rep(7.5,length(params$eps_a)),
       eps_rm=rep(Inf,length(params$eps_rm)-sum(is.na(map$eps_rm))),
       eps_f=rep(15,length(params$eps_a)-sum(is.na(map$eps_f))))
   return(list(U=U,L=L))
 }
 
 
 # source plotting functions
 source(here("src","SR_plotting_funcs.R"))
 
 
 ##   Fit models  ##
 
# begin dredge where all life histories same funcitonal form
 
 #load TMB model
 setwd(here("src","TMB"))
 TMB::compile("Stock_recruit.cpp")
 dyn.load(dynlib("Stock_recruit"))
 
 
#models to try all combinations of     
mods<-c(9, #linear
        10,#power
        6, #Beverton-holt
        1) #depensatory B-H I

#all combinations, where all streams within a given juvenile lfie histories have same functional form
mod_mat<-expand.grid(fry=mods,sum=mods,fall=mods,spring=mods)

#arrray to hold outputs
dredge_mat<-array(NA,dim=c(nrow(mod_mat),10,2))

start<-Sys.time()
for(i in 1:nrow(mod_mat)){
  for(j in 1:2){
    if(!is.na(dredge_mat[i,5,1])){next}
  params<-make_params_func(dat)
  mod<-NA
  fit<-NA
  
  dat$mod<- as.numeric(rep(mod_mat[i,],times=4))
  map<-make_map(dat)
  bounds<-make_bounds_func()        # make bounds
  
  #map$beta_a<-factor(rep(NA,4))
  #map$log_sigma_a<-factor(rep(NA,4))
  #map<-list()
  #map$eps_rm<-factor(rep(NA,12))
  mod<-TMB::MakeADFun(dat,params,random=c("log_S_hat","theta","gamma"),DLL="Stock_recruit",map=map,silent = TRUE)

  try(fit<-TMBhelper::fit_tmb(mod, newtonsteps = 1,lower = bounds$L,upper=bounds$U))
  try(if(is.na(fit)) next)
 # try(if(fit$Convergence_check=="The model is likely not converged") next )
  try(dredge_mat[i,1:4,j]<-as.numeric(mod_mat[i,]))
  try(dredge_mat[i,5,j]<-fit$AIC)
  try(dredge_mat[i,6,j]<-fit$max_gradient)
  try(dredge_mat[i,7:9,j]<-fit$number_of_coefficients)
  }
}


#replace failed first attempts at fitting with second attempt when succesfull
dredge_mat[,,1][which(is.na(dredge_mat[,,1]))]<-dredge_mat[,,2][which(is.na(dredge_mat[,,1]))]
elapsed<-Sys.time()-start
sum(!is.na(dredge_mat[,5,])) #number of models that converged
table(apply(dredge_mat[,5,],1,function(x)sum(is.na(x))))
summary(apply(dredge_mat[,5,],1,function(x)diff((x))),na.rm=T)
apply(dredge_mat[,5,],1,mean,na.rm=T)

#best 20 models
top_20_dredge<-as.data.frame(head(dredge_mat[order(dredge_mat[,5,1]),1:5,1],20))
colnames(top_20_dredge)<-c("Fry","Summer","Fall","Smolts","AIC")  
top_20_dredge$Delta_AIC<-round(top_20_dredge$AIC-min(top_20_dredge$AIC),3)
exp.delta<- exp(-0.5*(top_20_dredge$AIC - min(top_20_dredge$AIC,na.rm=T)))
weights <- exp.delta/sum(exp.delta,na.rm=T)
top_20_dredge$weights<-round(weights,4)
for (i in 1:4)top_20_dredge[,i]<-c("BH II",NA,NA,NA,NA,"BH",NA,NA,"linear","Power")[top_20_dredge[,i]]
top_20_dredge[,5]<-round(as.numeric(top_20_dredge[,5]),3)
apply(dredge_mat[order(apply(dredge_mat[,5,],1,mean,na.rm=T)),5,],1,function(x)sum(is.na(x)))
dredge_mat[order(apply(dredge_mat[,5,],1,mean,na.rm=T)),6,]
aic<-c(apply(dredge_mat[,5,],1,mean,na.rm=T),fit$AIC)
exp.delta <- exp(-0.5*(aic - min(aic,na.rm=T)))
weights <- exp.delta/sum(exp.delta,na.rm=T)
head(sort(weights,decreasing = TRUE),10)
dredge_mat[1:10,,]

save(top_20_dredge,file=here("results","dredge2.1.Rdata"))
load(file=here("results","Rdata","dredge2.Rdata"))

write.csv(top_20_dredge,here("results","dredge2.1.csv"))
# begin dredge where all life histories unique funcitonal form

#models to try all combinations of where each individual stream and life history can have its own functional form    
mods<-c(9, #linear
        5,# power
        10, # power
        6, #Beverton-holt
        2,
        1)#1) #depensatory B-H I

mods<-1:10
dredge_indivudal<-array(NA,dim=c(length(mods),3,4))


for ( i in 0:2){
  
  for( j in 1:4){
    
    dat<-make_dat_func(i,j)
    params<-make_params_func(dat)
   

    for (k in 1:length(mods)){
    
    dat$mod<- rep(mods[k],3)
    map<-make_map(dat)
    fit<-NA
    bounds<-make_bounds_func()
    mod<-TMB::MakeADFun(dat,params,random=c("log_S_hat","theta","gamma","beta_e"),DLL="Stock_recruit",map=map,silent = TRUE,lower = bounds$L,upper=bounds$U)
    try(fit<-TMBhelper::fit_tmb(mod, newtonsteps = 1))  
    try(dredge_indivudal[k,i+1,j]<-fit$AIC)
    }
  }
}
#best model for each
apply(dredge_indivudal,c(3:2),which.min)

# begin fitting best models

#load TMB model
setwd(here("src","TMB"))
TMB::compile("Stock_recruit.cpp")
dyn.load(dynlib("Stock_recruit"))

#make data and parameters
dat<-make_dat_func()
params<-make_params_func(dat)
#str(dat)
#str(params)


## begin model averaging ##

#matrix of four top models 
mod_mat<-rbind(c(10,10,6,6),
               c(1,10,6,6),
               c(10,1,6,6),
               c(1,1,6,6))

#fit and save each model
mod_outs<-list()
for ( i in 1:4){
  dat$mod<-rep(mod_mat[i,],times=3) # OR specify functional forms by life history
  map<-make_map(dat)                #make map
  mod<-TMB::MakeADFun(dat,params,random=c("log_S_hat","theta","gamma"),DLL="Stock_recruit",map=map,silent = TRUE)
  bounds<-make_bounds_func()        # make bounds
  fit<-TMBhelper::fit_tmb(mod, lower = bounds$L, upper=bounds$U, newtonsteps = 1, getJointPrecision = TRUE, getReportCovariance = TRUE)
  report<-mod$report()
  all_dat_plot<-plot_fun_pan()
  mod_outs[[i]] <- list(aic=fit$AIC,
                      report=report,
                      mod=mod,
                      fit=fit,
                      preds=all_dat_plot$pred_mat) 
  
}

#AIC weights
aic_vec<-unlist(lapply(mod_outs,function(x)x$aic))
exp.delta <- exp(-0.5*(aic_vec - min(aic_vec,na.rm=T)))
weights <- exp.delta/sum(exp.delta,na.rm=T)


#predicted spawners 
log_s_pred_pos<-which(names(mod_outs[[1]]$fit$SD$value)=="log_S_hat")
log_s_pred<-matrix(unlist(lapply(mod_outs,function(x)x$fit$SD$value[log_s_pred_pos])),ncol=4)
log_s_sd_pred<-matrix(unlist(lapply(mod_outs,function(x)x$fit$SD$sd[log_s_pred_pos])),ncol=4)%>% `colnames<-`(paste("sd",1:4))
#model averaged spawners adn SEs
MA_s<-log_s_pred %>% as.tibble() %>% rowwise() %>% mutate(MA_S=sum(c(V1,V2,V3,V4)*weights)) %>% bind_cols(as.tibble(log_s_sd_pred)) %>%   rowwise() %>% mutate(MA_S_sd=sum(weights*sqrt(c(`sd 1`,`sd 2`,`sd 3`,`sd 4`)^2+(c(V1,V2,V3,V4)-MA_S)^2)))
#predicted jvueniles
log_r_pred_pos<-which(names(mod_outs[[1]]$fit$SD$value)=="log(R_pred)")
log_r_pred<-matrix(unlist(lapply(mod_outs,function(x)x$fit$SD$value[log_r_pred_pos])),ncol=4)
log_r_sd_pred<-matrix(unlist(lapply(mod_outs,function(x)x$fit$SD$sd[log_r_pred_pos])),ncol=4) %>% `colnames<-`(paste("sd",1:4))
#model averaged juveniles and SEs and combine wiht spawners 
MA_r<- log_r_pred %>% as_tibble() %>% rowwise() %>% mutate(MA_R=sum(c(V1,V2,V3,V4)*weights)) %>% bind_cols(as.tibble(log_r_sd_pred)) %>%   rowwise() %>% mutate(MA_R_sd=sum(weights*sqrt(c(`sd 1`,`sd 2`,`sd 3`,`sd 4`)^2+(c(V1,V2,V3,V4)-MA_R)^2))) %>% ungroup() %>%  mutate(stream=c("Chiwawa","Nason","White")[as.numeric(dat$s_i)],LH=c("Fry","Summer","Fall","Smolt")[as.numeric(dat$l_i)]) %>% bind_cols(MA_s[as.numeric(dat$st_i),]) %>% mutate(LH=fct_relevel(LH,"Fry","Summer","Fall","Smolt"))

#predictions of expected juveniles over a range of spawner abundances
preds<-full_join(mod_outs[[1]]$preds,mod_outs[[2]]$preds,by=c("spawners","LH","stream")) %>% 
  full_join(mod_outs[[3]]$preds,by=c("spawners","LH","stream")) %>% 
  full_join(mod_outs[[4]]$preds, by=c("spawners","LH","stream")) %>% drop_na() %>% 
  rowwise() %>%  mutate(MA =exp(sum(log(c(juveniles.x,juveniles.y,juveniles.x.x,juveniles.y.y))*weights))) %>% 
  mutate(MA_pred_SE_synch =sum(weights*sqrt(c(theta.x,theta.y,theta.x.x,theta.y.y)^2+(log(c(juveniles.x,juveniles.y,juveniles.x.x,juveniles.y.y))-log(MA))^2)))%>% mutate(MA_pred_SE_total =sum(weights*sqrt(c(eps.x,eps.y,eps.x.x,eps.y.y)^2+(log(c(juveniles.x,juveniles.y,juveniles.x.x,juveniles.y.y))-log(MA))^2))) %>%  mutate(LH=fct_relevel(LH,"Fry","Summer","Fall","Smolt")) %>% replace_na(list(MA_pred_SE_synch=0,MA_pred_SE_total=0))

#Begin plot of spawners vs juveniles 
scale_juv<-100
#facet_wrap plot
SR_plot<-ggplot(data=preds,aes(x= spawners,y=MA/scale_juv))+facet_wrap(~LH+stream,scales = "free",nrow=4) + geom_ribbon(aes(ymin=MA*exp(-1.96*MA_pred_SE_total)/scale_juv,
                                                                                                               ymax=MA*exp(1.96*MA_pred_SE_total)/scale_juv,fill=rgb(.7,.1,.1,.2)))+
  geom_ribbon(aes(ymin=MA*exp(-1.96*MA_pred_SE_synch)/scale_juv,
                  ymax=MA*exp(1.96*MA_pred_SE_synch)/scale_juv,fill=rgb(.7,.1,.1,.21)))+geom_line(color="firebrick3",size=1.25)+
  geom_point(data=MA_r,aes(x=exp(MA_S),y=exp(MA_R)/scale_juv))+ 
  geom_linerange(data=MA_r,aes(x=exp(MA_S),y=exp(MA_R)/scale_juv,ymin=exp(qnorm(.025,MA_R,MA_R_sd))/scale_juv,ymax=exp(qnorm(.975,MA_R,MA_R_sd))/scale_juv))+
  geom_linerange(data=MA_r,aes(x=exp(MA_S),y=exp(MA_R)/scale_juv,xmin=exp(qnorm(.025,MA_S,MA_S_sd)),xmax=exp(qnorm(.975,MA_S,MA_S_sd)))) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.margin = margin(25.5, 5.5, 5.5, 5.5, "pt"),
    panel.spacing.y=unit(-.25,"cm"),
    legend.box.spacing=unit(1.5,"cm"))+
  xlab("Spawners")+ylab("Emigrants x100/ km")+
  scale_fill_manual(values=c(rgb(.7,.1,.1,.2),rgb(.7,.1,.1,.3)),labels=c("Idiosyncratic","Synchronous"),name="Prediction interval")+guides(fill = guide_legend(override.aes= list(alpha = c(0.2,.3))))


#grid_wrap plot to fet facet labels grobs from
p_grid<-ggplot(data=preds,aes(x= spawners,y=MA))+facet_grid(LH~stream,scales = "free")+geom_ribbon(aes(ymin=MA*exp(-1.96*MA_pred_SE_total),
                                                                                                       ymax=MA*exp(1.96*MA_pred_SE_total)),fill=rgb(.7,.1,.1,.2)) + 
  theme(
    plot.margin = margin(2, 2, 2, 2, "cm")
  ) 


# funciton to selectively remove certain grobs
gtable_filter_remove <- function (x, name, trim = TRUE){
  matches <- !(x$layout$name %in% name)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  if (trim) 
    x <- gtable_trim(x)
  x
}

#convert facet_wrap plot to Grob list
p_grid_tab<-ggplotGrob(p_grid)
#convert facet_grid plot to Grob list
SR_plot_tab<-ggplotGrob(SR_plot)

# remove bottom axes from all but smolt facets
SR_plot_filtered<-gtable_filter_remove(SR_plot_tab,
                                        name = c(paste0("axis-b-",rep(1:3,each=3) ,"-", 1:3),
                                                 paste0("strip-t-",rep(1:3,each=3) ,"-", 2:4)),
                                        trim = FALSE)

#add facet labels for columns (Stream)
SR_plot_filtered<-gtable::gtable_add_grob(SR_plot_filtered, p_grid_tab$grobs[grep('strip-t', p_grid_tab$layout$name)],t=rep(1,3), l= c(5,9,13))

#add facet labels for rows (LHs)
SR_plot_filtered<-gtable::gtable_add_grob(SR_plot_filtered, p_grid_tab$grobs[grep('strip-r', p_grid_tab$layout$name)],t=c(8,13,18,23),l=rep(16,4) )

#adjust placement of facet labels
SR_plot_filtered$heights[7]<-unit(-2,"pt")
SR_plot_filtered$widths[15]<-unit(-.35,"cm")

#render plots
grid::grid.newpage()
grid::grid.draw(SR_plot_filtered)

png(file=here("results","plots","spawn_em_MA.png"),units="in",height=5,width=6.5,res=300)
grid::grid.newpage()
grid::grid.draw(SR_plot_filtered)
dev.off()

#cowplot::plot_grid(SR_plot_filtered,SR_plot_filtered,nrow=2)


#model average covariate effects
beta_e_pos<-which(names(mod_outs[[1]]$fit$par)=="beta_e")
beta_e_pred<-matrix(c(unlist(lapply(mod_outs,function(x)x$fit$par[beta_e_pos])),
                       unlist(lapply(mod_outs,function(x)diag(x$fit$SD$cov.fixed)[beta_e_pos]))),ncol=8) %>% as_tibble()%>% rowwise() %>% mutate(MA=sum(c(V1,V2,V3,V4)*weights),MA_se=sum(weights*sqrt(c(V5,V6,V7,V8)+(c(V1,V2,V3,V4)-MA)^2)))

# processs variance explained by covariates
error_sds_pos<-which(names(mod_outs[[1]]$fit$par) %in% c("log_proc_sigma","log_sigma_theta"))
error_sds_pred<-matrix(unlist(lapply(mod_outs,function(x)x$fit$par[error_sds_pos])),ncol=4) %>% as_tibble() %>% rowwise() %>% mutate(MA=sum(c(V1,V2,V3,V4)*weights))


## fit all mdodels with no covariates ##

#make data and parameters
dat<-make_dat_func()
params<-make_params_func(dat)
#fit all mdodels with no covariates
mod_outs_NoCov<-list(4)
for ( i in 1:4){
  dat$mod<-rep(mod_mat[i,],times=3) # OR specify functional forms by life history
  map<-make_map(dat)                #make map
  map$beta_e<-factor(rep(NA,length(params$beta_e)))
  mod<-TMB::MakeADFun(dat,params,random=c("log_S_hat","theta","gamma"),DLL="Stock_recruit",map=map,silent = TRUE)
  bounds<-make_bounds_func()        # make bounds
  fit<-TMBhelper::fit_tmb(mod, lower = bounds$L,upper=bounds$U,newtonsteps = 1, getJointPrecision = TRUE,getReportCovariance = TRUE)
  report<-mod$report()
  all_dat_plot<-plot_fun_pan()
  mod_outs_NoCov[[i]]<-list(aic=fit$AIC,
                      report=report,
                      mod=mod,
                      fit=fit,
                      preds=all_dat_plot$pred_mat) 
  
}

#AIC weights
aic_vec_NoCov<-unlist(lapply(mod_outs_NoCov,function(x)x$aic))
exp.delta_NoCov <- exp(-0.5*(aic_vec_NoCov - min(aic_vec_NoCov,na.rm=T)))
weights_NoCov <- exp.delta_NoCov/sum(exp.delta_NoCov,na.rm=T)


#model average ICC
logit_ICC_pos<-which(names(mod_outs_NoCov[[1]]$fit$SD$value)%in%c("logit(ICC_sl)","logit(ICC_s)"))
logit_ICC_pred<-matrix(c(unlist(lapply(mod_outs_NoCov,function(x)x$fit$SD$value[logit_ICC_pos])),
                       unlist(lapply(mod_outs_NoCov,function(x)x$fit$SD$sd[logit_ICC_pos]))),ncol=8) %>% as_tibble() %>% rowwise() %>%  mutate(MA=sum(c(V1,V2,V3,V4)*weights_NoCov),MA_se=sum(weights*sqrt(c(V5,V6,V7,V8)^2+(c(V1,V2,V3,V4)-MA)^2)))


# percent of variance explained by covariates
error_sds_pos_NoCov<-which(names(mod_outs_NoCov[[1]]$fit$par) %in% c("log_proc_sigma","log_sigma_theta"))
error_sds_pred_NoCov<-matrix(unlist(lapply(mod_outs_NoCov,function(x)x$fit$par[error_sds_pos_NoCov])),ncol=4) %>% as_tibble() %>% rowwise() %>% mutate(MA=sum(c(V1,V2,V3,V4)*weights_NoCov))
                      
1-error_sds_pred$MA/error_sds_pred_NoCov$MA

C_tab<-matrix(1-(error_sds_pred$MA[1:12]^2+rep(error_sds_pred$MA[13:15]^2,each=4))/
  (error_sds_pred_NoCov$MA[1:12]^2+rep(error_sds_pred_NoCov$MA[13:15]^2,each=4)),ncol=3,dimnames = list(LH=c("Fry","Summer","Fall","Smolt"),stream=c("Chiwawa","Nason","White")))
  
write.csv(C_tab,file=here("results","C_tab.csv"))




#specify functional forms
#dat$mod<- mods[c(apply(dredge_indivudal,c(3:2),which.min))] #each indiviudal best form dredge with unique life history x stream funcitonal forms

#make data and parameters
dat<-make_dat_func()
params<-make_params_func(dat)

dat$mod<-rep(c(10,1,6,6),times=3) # OR specify functional forms by life history
map<-make_map(dat)                #make map
#map$beta_e<-factor(rep(NA,length(params$beta_e)))     #optional fix environmental covariate effects at 0
#map$hyper_beta_e<-factor(rep(NA,length(params$hyper_beta_e))) 
#map$log_sigma_e<-factor(rep(NA,length(params$log_sigma_e)))

#initialize model
mod<-TMB::MakeADFun(dat,params,random=c("log_S_hat","theta","gamma"),DLL="Stock_recruit",map=map,silent = TRUE)

bounds<-make_bounds_func()        # make bounds

#mod$gr()       #check initial gradient 
#mod$fn()       #check initial function value

#fit model 
fit<-TMBhelper::fit_tmb(mod, lower = bounds$L,upper=bounds$U,newtonsteps = 1, getJointPrecision = TRUE,getReportCovariance = TRUE)
fit$AIC         # AIC
report<-mod$report() #Report
tab_SR_1<-plot_fun_pan() #plot all 
# plot functional relationships and year specific values of emigrants and spawners
png(file=here("results","plots","spawn_em.png"),units="in",height=5,width=6.5,res=300)
tab_SR_1<-plot_fun_pan()
dev.off()


#weird model averaging crap. 
thing<-thing %>% full_join(tab_SR_1$pred_mat %>% rename(juveniles_4=juveniles))
thing<-thing %>% mutate(MA_juv=juveniles*weights[1]+
                           juveniles_2*weights[2]+
                           juveniles_3*weights[3]+
                           juveniles_4*weights[4])
thing<-thing %>% mutate(LH=fct_relevel(LH,"Fry","Summer","Fall","Smolt"))

ggplot(data=thing,aes(x=spawners,y=MA_juv))+geom_line()+geom_point(data=tab_SR_1$pred_mat_2%>% mutate(LH=fct_relevel(LH,"Fry","Summer","Fall","Smolt")),aes(x=exp(spawners),y=exp(juveniles)))+facet_wrap(~stream+LH,scales = "free")


thing<-tab_SR_1$pred_mat

300*(weights[1]+weights[3])+100*(weights[2]+weights[4])

#end of weird model averaging crap

test<-mvtnorm::rmvnorm(1,fit$SD$value,fit$SD$cov)

set.seed(1)
apply(mvtnorm::rmvnorm(100000,fit$SD$par.fixed[1:5],fit$SD$cov.fixed[1:5,1:5]),2,sd)

#posterior means
report<-mod$report()
report$sigma_delta^2 #year
report$sigma_kappa # life history by year
var_theta_noCov<-report$sigma_theta^2 # stream by brood year  # to-do: compare with stream by year
var_sig_noCov<-report$proc_sigma^2  # residual

tot_var_noCov<-rep(var_theta_noCov,each=4)+var_sig_noCov

var_theta<-report$sigma_theta^2 # stream by brood year  # to-do: compare with stream by year
var_sig<-report$proc_sigma^2  # residual

tot_var<-rep(var_theta,each=4)+var_sig


matrix(1-(tot_var/tot_var_noCov),4)

matrix(report$alpha,4)
matrix(report$R_max,4)
matrix(report$fifty,4)
report$sigma_log_fifty
report$sigma_log_R_max
fit$opt$par
plot(report$beta_e)
report$beta_e_2
plot(report$beta_e_2[order(dat$beta_e_i)])

plot(report$beta_e)
## ICC 
fit$SD$value[names(fit$SD$value)=="ICC_sl"] #life history by stream specific
fit$SD$sd[names(fit$SD$value)=="ICC_sl"]  #sd
a<-fit$SD$value[names(fit$SD$value)=="logit(ICC_s)"] #across- life-history average (by stream)
b<-fit$SD$sd[names(fit$SD$value)=="logit(ICC_s)"] #sd
plogis(qnorm(.025,a,b))
plogis(qnorm(.975,a,b))


# plot functional relationships and year specific values of emigrants and spawners
png(file=here("results","plots","egg_em_dif_ax_fit.png"),units="in",height=5,width=6.5,res=300)
tab_SR_1<-plot_fun_pan()
dev.off()

#spawners and juvenile sby year
png(here("results","plots","spaw_em.png"),units="in",res=300,height=5,width=7)
plot_spawn_em_year()
dev.off()


##expected emigrants vs. spawners
png(here("results","plots","juv_LH.png"),units="in",res=300,height=4,width=5)
expec_spawn_em_ploft_func()
dev.off()

## plot ICC 
png(here("results","plots","ICC_plot.png"),units="in",res=300,height=4,width=5)
plot_ICC_func()
dev.off()

# plot environmental covariate coefficients
png(here("results","plots","coef_plot_all.png"),units="in",res=300,height=4,width=5)
cov_coef_plot_func(FALSE)
dev.off()





#bootstrapping ICC, functions below
ICC_boot<-bootstrap_ICC()
range(ICC_boot)
dim(ICC_boot)
dev.new()
par(mfrow=c(5,3))
apply(ICC_boot,2,hist)

thing<-ICC_boot %>% as.tibble() %>% `colnames<-`(paste(rep(1:5,times=3),rep(1:3,each=5),sep=".")) %>% pivot_longer(1:15) %>% mutate(stream=rep(c("Chiwawa","Nason","White"),each=5,length.out=150000),LH=rep(c("Fry","Summer","Fall","Smolt","Average"),length.out=150000)) %>% mutate(LH=fct_relevel(LH,"Fry","Summer","Fall","Smolt","Average"))

ICC_plot_boot<-ggplot(data=thing, aes(x=stream,y=value,fill=LH))+geom_boxplot(position = position_dodge(.8),width=.5,outlier.colour = NA)+scale_fill_manual(values=cols) + xlab("") +ylab("ICC")



which(names(fit$par)=="log_proc_sigma")
which(names(fit$par)=="log_sigma_theta")

fit$SD$cov.fixed



#Function to bootstrap model parameters 
rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
  set.seed( random_seed )
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L = Matrix::Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  return(mu + z)
}

#function to boostrap ICC
bootstrap_ICC<-function(dat, last_best=mod$env$last.par.best , joint_precis=fit$SD$jointPrecision ,n_sim=10000, seed=1234){
  test_sim<-rmvnorm_prec(last_best ,
                         joint_precis ,
                         n_sim,seed)
  test_sim[9:23,][test_sim[9:23,]<.0001]<- .0001
  
  row_proc_sigma<-which(names(last_best)=="log_proc_sigma")
  row_sigma_theta<-which(names(last_best)=="log_sigma_theta")
  
  out<-matrix(NA,n_sim,15)
  for ( i in 1:3){
    for(j in 1:4){
      out[,(i-1)*5+j]<-(test_sim[row_sigma_theta[i],])^2/
        ((test_sim[row_sigma_theta[i],])^2+(test_sim[row_proc_sigma[(i-1)*4+j],])^2)
    }
    out[,(i-1)*5+5]<-apply( out[,((i-1)*5+1):((i-1)*5+4)],1,mean)
  }
  
 return(out)
  
}



cov_var



#smooth hockey stick asymptotic max recruitment
report$alpha*10*report$R_max*(1+exp(-.1))*(.1+log(1+exp(-.1)))


## Look at post hoc correlations in process erros##

##correlation idiosyncratic
cor_test<-matrix(NA,nrow=dat$n_t,ncol=dat$n_sl)
for ( i in 1:dat$n_i){
  cor_test[as.numeric(dat$t_i)[i],as.numeric(dat$sl_i)[i]]<-report$gamma[i]
}

cor_test_2<-cor(cor_test,use="pairwise.complete.obs")
lh_s_names<-paste(rep(c("Fry","Summer","Fall","Smolt"),times=3),rep(c("Chiwawa","Nason","White"),each=4),sep=".")
dimnames(cor_test_2)<-list(lh_s_names,lh_s_names)
corrplot::corrplot(cor_test_2,type="lower")


#correlation stream level
cor_test<-matrix(NA,nrow=dat$n_t,ncol=dat$n_s)
for ( i in 1:dat$n_i){
  cor_test[as.numeric(dat$t_i)[i],(as.numeric(dat$s_i)[i])]<-report$theta[dat$st_i[i]]
}

cor_test_2<-cor(cor_test,use="pairwise.complete.obs")
s_names<-c("Chiwawa","Nason","White")
dimnames(cor_test_2)<-list(s_names,s_names)
corrplot::corrplot(cor_test_2,type="lower")





View(cor_test_2)
dim(cor_test_2)
mean(cor_test_2)
mean(cor_test_2[1:4,1:4])
mean(cor_test_2[5:8,5:8])
mean(cor_test_2[9:12,9:12])


cor_test<-matrix(NA,nrow=dat$n_t,ncol=3)
for ( i in 1:dat$n_i){
  cor_test[as.numeric(dat$t_i)[i],as.numeric(dat$sl_i)[i]]<-report$gamma[i]
}


cor(cor_test,use="pairwise.complete.obs")

#calculate Coefficients of Variation
juvenile_dat<-bind_cols(Juveniles  =report$R_pred/100,stream=c("Chiwawa","Nason","White")[dat$s_i],BY=dat$BY,LH=c("Fry","Summer","Fall","Smolt")[dat$l_i]) %>% mutate(LH=fct_relevel(LH,"Fry","Summer","Fall","Smolt"))   #drop first and last brood years which have incomplete juvenile data
  
  
  sum_juv<-juvenile_dat %>% group_by(BY,LH) %>% summarise(Juveniles = sum(Juveniles)/n(),count=n()) %>% ungroup() %>% mutate(stream="Average") %>% filter(count==3)

juvenile_dat<-bind_rows(juvenile_dat,sum_juv) %>% group_by(BY,stream) %>% 
  group_by(stream) %>% mutate(max_BY=max(BY),min_BY=min(BY)) %>% ungroup() %>% mutate(min_overall=max(min_BY)) %>% filter(BY>min_overall&BY<max_BY) %>% group_by(stream,BY) %>% summarize(total=sum(Juveniles),n=n()) %>% group_by(stream) %>% summarize(mean=mean(total),sd=sd(total),cv=sd(total)/mean(total),n=n())

#spawner estimatesx by BY and stream
dat$S_obs %>% mutate(stream=c("chiw","nas","Whi")[dat$s_i],BY=dat$BY) %>% distinct()
