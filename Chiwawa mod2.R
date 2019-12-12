pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}#end of function

pkgTest("here")
pkgTest("mvtnorm")
pkgTest("TMB")


###function to make data and inits for model
mod_ins<-function(dat, mod_lifestage, start_day, end_day,  p_mod_formula,use_NB){
  
  #subset to data rows which have catch for the lifestage being modeled or efficiency trial information
  dat<-subset(dat,(!is.na(rel)|lifestage==mod_lifestage)&DOY<=end_day)
  dat$log_scale_Dis<-scale(log(dat$dis))
  #row indices for efficiency trials
  effic_indices<-which(!is.na(dat$rel))
  
  #number of efficiency trials
  Ntrials<-length(effic_indices)
  
  # #released per trial
  rel<- dat$rel[effic_indices]
  # # recaptured per trial
  rec<- dat$recap[effic_indices]
  
  #number of days to estimate migrants
  max_day<-length(start_day:end_day)
  
  
  #Catch Data
  
  #index of rows with catch observations
  catch_index<-which(!is.na(dat$Catch)&
                       dat$lifestage==mod_lifestage)
  # number of years of catch data
  Nyears<-length(unique(dat$year[catch_index]))
  
  #Vector of year index of catch data
  seriesFac<-as.numeric(dat$year[catch_index])-min(as.numeric(dat$year[catch_index]))+1
  
  # Vector of Catch of subyearlings (fry and parr combined)
  catch<-dat[catch_index,"Catch"]
  
  # day of year of catches
  catch_DOY<-dat$DOY[catch_index]-start_day+1
  
  #Total length of catch vector
  N_trap<-length(catch)
  
  #model matrix for predicting daily capture efficiency
  p_pred_mat<-model.matrix(p_mod_formula,dat)
  
  if (sum(is.na(p_pred_mat))>0) stop("missing values in capture efficiency design matrix")
  
  # list of data to be fed to model
  mod_data <- list(Ntrials=Ntrials,rec=rec,rel=rel,efficIndex=effic_indices-1, pDat=p_pred_mat,Catch=catch,catch_DOY=catch_DOY,seriesFac=seriesFac,N_trap=N_trap, max_day=max_day,Nyears=Nyears,use_NB=use_NB)
  
  
  #initial values
  #coefficients for capture efficiency model
  pCoefs <-rnorm(ncol(p_pred_mat))
  #random year effects on efficiency
  rand_year<-rnorm(Nyears)
  #random day effects on efficiency 
  rand_day<-rnorm(length(seriesFac))
  #log stnadard deviation of random year effects
  log_sigma_year<-rnorm(1,-1,.5)
  #log stnadard deviation of random day effects
  log_sigma_day<-rnorm(1,-1,.5)
  
  # intercepts of AR(1) process for log-mean daily outmigrants
  mu_M<-runif(Nyears,2,7)                                                 
  # AR(1) process error SD for log-mean daily outmigrants
  log_sigma_M<-runif(Nyears,-1,1)
  
  # AR(1) coefficient for log-mean daily outmigrants
  logit_phi_M<-rnorm(Nyears)     
  
  #latent state z-scored log(expected migrants)
  log_M_hat_z<-matrix(rnorm(Nyears*max_day),nrow=max_day,ncol=Nyears)   
  
  #variance of negative binomial
  log_phi=rnorm(ifelse(use_NB==1,Nyears,0))
  
  #p paramater of negative binomial used in the "...NB_2.cpp" script, Which used the same paramaterixation as R.  
  logit_p_NB=rnorm(ifelse(use_NB==1,Nyears,0))
  
  #paramater list
  mod_pars=list(pCoefs=pCoefs,rand_year=rand_year,rand_day=rand_day,log_sigma_year=log_sigma_year,log_sigma_day=log_sigma_day,mu_M=mu_M, logit_phi_M=logit_phi_M, log_sigma_M=log_sigma_M,log_M_hat_z=log_M_hat_z,log_phi=log_phi,logit_p_NB=logit_p_NB)
  
  
  
  #bounds
  lower=c(rep(-5,times=length(pCoefs)), -5,
          rep(-7,Nyears),rep(-5,Nyears),rep(-5,times=Nyears),rep(-5,times=Nyears))
        
          #rep(-5,Nyears),rep(-5,times=max_day*Nyears))
  
  
  
  upper=c(rep(5,times=length(pCoefs)), 5,
          rep(7,Nyears),rep(5,Nyears),rep(5,times=Nyears),rep(5,times=Nyears))
  
  #rep(-5,Nyears),rep(-5,times=max_day*Nyears))
  
  #return
  list(mod_data=mod_data,mod_pars=mod_pars,lower=lower,upper=upper)
}


###Function to bootstrap recruits

bootstrap_migrants<-function(n_Draws=2000,dat,SD_obj,redds,breaks=NA){
  library(mvtnorm)
  #Bootstrap migrants
  N_years<-dat$mod_data$Nyears
  N_days<-ncol(SD_obj$cov)/ dat$mod_data$Nyears
  out<-array(NA,dim=c(n_Draws,N_days,N_years))
  
  log_M_hat_array<-array(NA, dim=c(N_days,N_years,3))
  log_M_hat_SD<-summary(SD_obj)[]
  m_hat_SD<-log_M_hat_SD[rownames(log_M_hat_SD)=="log_M_hat",]
  
  
  for ( i in 1:N_years){
    
    draws<-rmvnorm(n_Draws , SD_obj$value[((1+(i-1)*N_days):(i*N_days))], SD_obj$cov[((1+(i-1)*N_days):(i*N_days)),((1+(i-1)*N_days):(i*N_days))] )
    
    out[,,i]<-draws
    
    log_M_hat_array[,i,1]<-(m_hat_SD[((1+(i-1)*N_days):(i*N_days)),1])
    log_M_hat_array[,i,2]<-qnorm(.025,m_hat_SD[((1+(i-1)*N_days):(i*N_days)),1],m_hat_SD[((1+(i-1)*N_days):(i*N_days)),2])
    log_M_hat_array[,i,3]<-qnorm(.975,m_hat_SD[((1+(i-1)*N_days):(i*N_days)),1],m_hat_SD[((1+(i-1)*N_days):(i*N_days)),2])
    
  }
  
  
  #break into groups
  num_groups<-length(breaks)+1-sum(is.na(breaks))
  groups<-list()
  for ( i in 1:num_groups){
    first_day<-ifelse(i==1,1,breaks[i-1]+1)
    end_day<-ifelse(i==num_groups,N_days,breaks[i])
    groups[[i]]<-apply((exp(out[,first_day:end_day,])),c(1,3),sum)
    log_means<-apply(log(groups[[i]]),2,mean)
      log_sds<-apply(log(groups[[i]]),2,sd)
    
    plot(redds,exp(log_means),main=paste(first_day,"-",end_day),ylab="recruits")
    
    segments(redds,exp(log_means+2*log_sds),redds,exp(log_means-2*log_sds))
    
  }
  return(list(groups,log_M_hat_array))
}
########End of funcion #####################


###Function to fit a spawner-emigrant transition model

fit_SR<-function(R_obs,S_obs,pHOS,Model,SD=FALSE,N_years){
  compile("Stock_recruit.cpp")
  dyn.load(dynlib("Stock_recruit"))
  
  
  SR_dat<-list(R_obs=R_obs,                 # Observed recruits (posterior from juvenile modle)
               S_obs=S_obs,               # Observed spawners (from redd counts)
   
               S_obs_cv_hyper_mu=0.05,           #spawner observation error hyper mean
               S_obs_cv_hyper_sd=0.01,        #spawner observation error hyper variance
               Model=Model)
  
  
  mean_R<-colMeans(R_obs)
  SR_pars<-list(beta=.05,log_R_hat=log(mean_R),      # latent true number of juveniles
                log_R_obs_sd=log(apply(R_obs,2,sd)),  # juvenile observation error
                log_S_hat=log(as.numeric(S_obs)),  # latent true number of spawners
                log_S_obs_cv=log(0.05),                         # Spawner observation error
                log_alpha=ifelse(Model==3,1.5,log(60)),         # intrinsic productivity
                log_R_max=ifelse(Model==3,5,log(60000)),          # Asymptotic maximum recruitment
                log_proc_sigma=log(.4),  # process error standard deviation
                log_d=0)
  if(Model==2|Model==4){
    map<-list(log_d=factor(1))
  }else{
    if(Model==5){
      map<-list(log_d=factor(NA),log_R_max=factor(NA))}else{
        map<-list(log_d=factor(NA))
      }}
  #browser()
  SR<-MakeADFun(SR_dat,SR_pars,random = c("log_R_hat","log_S_hat"),DLL="Stock_recruit", map = map,silent = T)
  lower<-c(rep(-20,length(mean_R)),rep(-50,length(mean_R)),rep(-10,length(mean_R)),-50,-10,-50,-50,-50)
  upper<-c(rep(log(max(as.numeric(S_obs)*10000)),length(mean_R)),
           rep(500,length(mean_R)),
           rep(500,length(mean_R)),5
           ,log(2500),
           log(max(as.numeric(S_obs))*5000),
           50,
           50)
  
  SR_fit<-nlminb(SR$par,SR$fn,SR$gr, 
                 control=list(rel.tol=1e-12,eval.max=1000000,
                              iter.max=10000))#,lower=lower,upper=upper)
  
  Npar<-ifelse(Model==2|Model==4,3,ifelse(Model==5,1,2))
  Nsamp<-N_years
  AICc<-SR$fn()*2+2*Npar+(2*Npar*(Npar+1))/(Nsamp-Npar-1)
  
  SR_Rep<-SR$report()
  
  
  plot(SR_Rep$S_hat,SR_Rep$R_hat,col="red")
  points(S_obs,mean_R)
  points(SR_Rep$S_hat,SR_Rep$R_pred,col="blue")
  
  if(isTRUE(SD)){
    SR_sd<-sdreport(SR)
    
  }
  
  out<-list(AICc=AICc,SR_Rep=SR_Rep)
  
  if(isTRUE(SD)){
    out[[3]]<-SR_sd
  }
  
  return(out)
  
}
###################### End Function ##############################


###Function to compare Spawner-emigrant candidate mdoels

SR_compare<-function(R_dat,S_dat,pHOS,N_years){
  out<-list()
  for ( i in 1:5){
    out[[i]]<-fit_SR(R_obs=R_dat,S_obs=S_dat,pHOS=pHOS,Model=i,N_years = N_years)
    print(out[[i]][[1]])
  }
  
  return(out)
}



###  Function to plot SR

plot_SR<-function(fit_list,model,main_title,Riv_length,Y_max,y_ax){
  s_hats<-fit_list[[model]]$SR_Rep$S_hat
  r_hats<-fit_list[[model]]$SR_Rep$R_hat
  Riv_length<-Riv_length*1000
  plot(s_hats,r_hats/Riv_length,xlim=c(0,max(s_hats)*1.1),ylim=c(0,Y_max/1000),main=main_title,xlab="",ylab="",yaxt="n",font=2,cex=1.2,type="n")#max(r_hats*2)/Riv_length)
  axis(2,labels=y_ax,font=2,cex=1.2)
  
  alpha<-fit_list[[model]]$SR_Rep$alpha
  R_max<-fit_list[[model]]$SR_Rep$R_max
  d<-fit_list[[model]]$SR_Rep$d
  predS<-seq(0,max(s_hats)*1.3,by=5)
  
  if(model==1) R_pred<-(alpha*predS)/(1+alpha*predS/R_max)
  
  if(model==2) R_pred<-(R_max*predS^d)/(alpha+predS^d)
  
  if(model==3) R_pred<- R_max*predS^alpha
  
  if(model==4) R_pred<-R_max*(1-exp(-((predS/alpha)^d)))
  
  if(model==5) R_pred<-predS*alpha
  
  proc_sigma<-fit_list[[model]]$SR_Rep$proc_sigma
  upper_R_pred<-exp(qnorm(.975,log(R_pred),proc_sigma))/Riv_length
  lower_R_pred<-exp(qnorm(.025,log(R_pred),proc_sigma))/Riv_length
  polygon(x=c(predS,rev(predS)),y=c(upper_R_pred,rev(lower_R_pred)),border =    FALSE,col="grey")
  
  points(predS,R_pred/Riv_length,type="l",lwd=2)
  points(s_hats,r_hats/Riv_length,pch=19)
  
  R_obs_sd<-fit_list[[model]]$SR_Rep$R_obs_sd
  S_obs_sd<-fit_list[[model]]$SR_Rep$S_obs_sd
  segments(exp(qnorm(.975,log(s_hats),S_obs_sd)),r_hats/Riv_length,exp(qnorm(.025,log(s_hats),S_obs_sd)),r_hats/Riv_length)
  segments(s_hats,exp(qnorm(.975,log(r_hats),R_obs_sd))/Riv_length,s_hats,exp(qnorm(.025,log(r_hats),R_obs_sd))/Riv_length)
  box()
}


###--------------------------------------------------------------------
#                       End of Functions loading
###--------------------------------------------------------------------


###Load data
if(file.exists(here("src","Chiw_dat_2.Rdata"))){
  
  load(here("src","Chiw_dat_2.Rdata"))
  
}else{
  source(here("src","Chiwawa Data Proc 2.R"))
  chiw_dat<-Chiw_dat_Proc()
  save(chiw_dat,file=here("src","Chiw_dat_2.Rdata"))
}

###load model

require(TMB)
#TMB model
setwd(here("src","TMB"))
compile("multi_year_independant_1NB_2_rand_p.cpp")
dyn.load(dynlib("multi_year_independant_1NB_2_rand_p"))

###Process data for models


#range of days of year in data
range_of_days<-range(chiw_dat$DOY)

#subyearlings
chiw_SBC_ins<-mod_ins(chiw_dat,"SBC",range_of_days[1],range_of_days[2],as.formula(~lifestage+log_scale_Dis*Position2-1),1)

#Yearlings
chiw_YCW_ins<-mod_ins(chiw_dat,"YCW",range_of_days[1],181,as.formula(~lifestage+log_scale_Dis*Position2-1),1)


### specifiy and optimize models
#specify
#subyearlings
model_ch_subs_NB2<- MakeADFun(chiw_SBC_ins$mod_data, chiw_SBC_ins$mod_pars,  
                              random=c("log_M_hat_z","rand_year"),
                              DLL="multi_year_independant_1NB_2_rand_p",silent=T)

#yearlings
model_ch_yrlngs_NB2<- MakeADFun(chiw_YCW_ins$mod_data, chiw_YCW_ins$mod_pars,  
                              random=c("log_M_hat_z","rand_year"),
                              DLL="multi_year_independant_1NB_2_rand_p",silent=T)


#optimize
#subyearlings
model_ch_subs_fit_NB2<- nlminb(model_ch_subs_NB2$par, model_ch_subs_NB2$fn, model_ch_subs_NB2$gr, 
                               control=list(rel.tol=1e-12,eval.max=1000000,
                                            iter.max=10000),upper=chiw_SBC_ins$upper,lower=chiw_SBC_ins$lower)

model_ch_subs_fit_NB2$objective
model_ch_subs_NB2$fn()

#quick check on results
rep_ch_subs_NB2<-model_ch_subs_NB2$report()
log_m_hat_ch_subsNB2<-rep_ch_subs_NB2$log_M_hat
plot(rowMeans(exp(log_m_hat_ch_subsNB2)),type="l",ylab="",xlab="",xaxt="n")
abline(v=87,col="red")
abline(v=215,col="red")
days<-c(91,182,274,366,(365+91))-50
labs<-c("Apr","Jul","Oct","Jan","Apr")
axis(1,at=days,labels=labs)

#Yearlings
model_ch_yrlngs_fit_NB2<- nlminb(model_ch_yrlngs_NB2$par, model_ch_yrlngs_NB2$fn, model_ch_yrlngs_NB2$gr, 
                               control=list(rel.tol=1e-12,eval.max=1000000,
                                            iter.max=10000),upper=chiw_YCW_ins$upper,lower=chiw_YCW_ins$lower)

model_ch_yrlngs_fit_NB2$objective
model_ch_yrlngs_NB2$fn()


#quick check on results
rep_ch_yrlngs_NB2<-model_ch_yrlngs_NB2$report()
log_m_hat_ch_yrlngsNB2<-rep_ch_yrlngs_NB2$log_M_hat
plot(rowMeans(exp(log_m_hat_ch_yrlngsNB2)),type="l",ylab="",xlab="",xaxt="n")
days<-c(91,182,274,366,(365+91))-50
labs<-c("Apr","Jul","Oct","Jan","Apr")
axis(1,at=days,labels=labs)


### Generate distributions of total migrants per period by drawing samples of log daily migrants from a multivariate normal, exponenetiating and summing.

#load chiwawa spawners data for plotting vs juvenile migrants
Chiw_redds<-read.csv(here("data","Chiwawa","ChiwEscapement.csv"))
Chiw_sub_spawners<-rowSums(Chiw_redds[8:29,2:3])
Chiw_yrlng_spawners<-rowSums(Chiw_redds[7:28,2:3])

#calculate covariance matrix by inverting the hessian and via the delta methods
ch_subs_SD_NB2<-sdreport(model_ch_subs_NB2)
ch_yrlngs_SD_NB2<-sdreport(model_ch_yrlngs_NB2)

#generate distributions and plot
chiw_subs<-bootstrap_migrants(n_Draws=2000,dat=chiw_SBC_ins,SD_obj=ch_subs_SD_NB2,redds=Chiw_sub_spawners,breaks=c(87,214))

chiw_yrlngs<-bootstrap_migrants(n_Draws=2000,dat=chiw_YCW_ins,SD_obj=ch_yrlngs_SD_NB2,redds=Chiw_yrlng_spawners)


all_boots_trans<-array(unlist(list(t(chiw_subs[[1]][[1]]),t(chiw_subs[[1]][[2]]),t(chiw_subs[[1]][[3]]),t(chiw_yrlngs[[1]][[1]]))),dim=c(22,2000,4))

hist((all_boots_trans[15,,1]))

lognorm_pars<-list(log_means=apply(log(all_boots_trans),c(1,3),mean),
log_sds=apply(log(all_boots_trans),c(1,3),sd))
save(lognorm_pars,file="C:/Users/Mark Sorel/Documents/558/final project/LCM_lite/lognorm_pars.Rdata")      


dim(all_boots_trans)
save(all_boots_trans,file="C:/Users/Mark Sorel/Documents/558/final project/LCM_lite/emigrants.Rdata")                                                          
mean(all_boots_trans[2,,4])                                                                     

mean(chiw_yrlngs[[1]][[1]][,2])
all_LHs<-array(c(unlist(chiw_subs[[1]]),unlist(chiw_yrlngs[[1]])),dim=c(2000,22,4))


mean(all_LHs[,2,4])


all_boots_trans<-matrix(unlist(list(t(chiw_subs[[1]][[1]]),t(chiw_subs[[1]][[2]]),t(chiw_subs[[1]][[3]]),t(chiw_yrlngs[[1]][[1]]))),dim=c(22,2000,4))



### compare different spawner->emigrant transtions for different life histories
Chiw_fry_SR<-SR_compare(chiw_subs[[1]][[1]],Chiw_sub_spawners,pHOS=NULL,N_years=22)

Chiw_summer_SR<-SR_compare(chiw_subs[[1]][[2]],Chiw_sub_spawners,pHOS=NULL,N_years=22)

Chiw_fall_SR<-SR_compare(chiw_subs[[1]][[3]],Chiw_sub_spawners,pHOS=NULL,N_years=22)

Chiw_smolt_SR<-SR_compare(chiw_yrlngs[[1]][[1]],Chiw_yrlng_spawners,pHOS=NULL,N_years=22)


###examine residuals of best models

#fry: model 4-Weibull CDF
fry_weib_resid<-Chiw_fry_SR[[4]]$SR_Rep$resids
plot(fry_weib_resid)
acf(fry_weib_resid)


#summer: model 4-Weibull CDF
summer_weib_resid<-Chiw_summer_SR[[4]]$SR_Rep$resids
plot(summer_weib_resid)
acf(summer_weib_resid)

#fall: model 4-Weibull CDF
fall_weib_resid<-Chiw_fall_SR[[4]]$SR_Rep$resids
plot(fall_weib_resid)
acf(fall_weib_resid)

#smolt: model 1- B-H
smolt_weib_resid<-Chiw_smolt_SR[[1]]$SR_Rep$resids
plot(smolt_weib_resid)
acf(smolt_weib_resid)

resid_cor<-cor(data.frame(fry_weib_resid,summer_weib_resid,fall_weib_resid,smolt_weib_resid))



trans_mod_par_func<-function(thing){
list(alpha=thing$alpha,
     R_max=thing$R_max,
     d=thing$d,
     R_obs_sd=thing$R_obs_sd,
     S_obs_cv=thing$S_obs_cv,
     proc_sigma=thing$proc_sigma,
     R_hat=thing$R_hat,
     S_hat=thing$S_hat)
}


outs_list<-

par_outs<-lapply(list(Chiw_fry_SR[[4]]$SR_Rep,Chiw_fall_SR[[4]]$SR_Rep,Chiw_summer_SR[[4]]$SR_Rep,Chiw_smolt_SR[[1]]$SR_Rep),trans_mod_par_func)
save(par_outs,file="C:/Users/Mark Sorel/Documents/558/final project/LCM_lite/par_outs.Rdata")

function(){
     inits_list<-list(S_hat,
                 S_hat_CV=.05,
                 em_hat=,
                 em_hat_CV,
                 em_proc_cov=,
                 alpha=,
                 R_max=,
                 d=
                   )
  } 
  for ( i in 1:4){
    
    inits_list$em_hat<-cbind
  
  
}


save()




str(par_outs)

trans_mod_par_func(Chiw_smolt_SR[[1]]$SR_Rep)



all_boots_trans<-matrix(unlist(list((chiw_subs[[1]][[1]]),(chiw_subs[[1]][[2]]),(chiw_subs[[1]][[3]]),(chiw_yrlngs[[1]][[1]]))),nrow=2000)


apply(log(all_boots_trans),2,function(x){qqnorm(x);qqline(x,col="red")})

apply(all_boots_trans,2,function(x){
  thing<-hist(x,breaks=30);par(new=T);plot(seq(thing$breaks[1],max(thing$breaks),100), dnorm(log(seq(thing$breaks[1],max(thing$breaks),100)),mean(log(x)),sd(log(x))),axes=TRUE,type="l",col="red",xlim=range(thing$breaks))
})


save(all_boots_trans,file=here("all_boots_trans_12-1-19.Rdata"))


dim(all_boots_trans)
colMeans(log(all_boots_trans))
apply(log(all_boots_trans),2,sd)
exp(SR_4$env$last.par.best)

library(TMB)
library(here)
library(tidyverse)
setwd(here("src","TMB"))
compile("LCM_lite4.cpp")
dyn.load(dynlib("LCM_lite4"))

load(here("all_boots_trans_12-1-19.Rdata"))
Chiw_redds<-read.csv(here("data","Chiwawa","ChiwEscapement.csv"))

spawners<-rowSums(Chiw_redds[7:29,2:3])
wild_spawners<-Chiw_redds[7:29,2]
hatch_spawners<-Chiw_redds[7:29,3]
phos<-Chiw_redds[7:29,4]

#hatchery broodstock removals
hatch_BS<-read.csv(here("data","Chiwawa","hatch_broodstock.csv"))

head(hatch_BS)

wild_brood_rmoveal<-hatch_BS[7:29,2]

pNOB<-hatch_BS[7:29,2]/rowSums(hatch_BS[7:29,c(2,4)]) 
pNOB[is.na(pNOB)]<-1

#carcass origins
carc_origins<-read.csv(here("data","Chiwawa","carc_phos.csv"))

head(carc_origins)

hatch_carcs<-carc_origins[3:25,3]
total_carcs<-rowSums(carc_origins[3:25,2:3])

pNI<-pNOB/(pNOB+test$pHOS) 
plot(pNI)
mean(pNI)


#Survival data
surv_dat<-read.csv(here("data","Chiwawa","surv_dat.csv"))
surv_dat2<-surv_dat %>% mutate(Length.mm=(Length.mm-mean(Length.mm))/sd(Length.mm)) %>% mutate(surv_TUM=as.numeric(surv_TUM)) %>% subset(brood_year_guess<=2014)%>%droplevels()

#PIT tag age at return
PIT_adult_age<-read.csv(here("data","Chiwawa","LH_age_BY"))

PIT_age_LHs<-as.numeric(factor(PIT_adult_age$life_stage_guess))-1
PIT_age_BYs<- PIT_adult_age$brood_year-1995
PIT_ages<-as.matrix(PIT_adult_age[,3:5] )
pit_props<-PIT_ages/matrix(rowSums(PIT_ages),nrow(PIT_ages),ncol(PIT_ages))
plot(pit_props[,1],type="o",ylim=c(0,1))
points(pit_props[,2],type="o",col="blue")
points(pit_props[,3],type="o",col="red")


#carcass (wild only) age at return
Carc_adult_age<-read.csv(here("data","Chiwawa","carcass ages.csv"))%>%mutate("cnt3"=(X3*n),"cnt4"=(X4*n),"cnt5"=(X5*n))%>%
  slice(2:24)%>%
  select(starts_with("cnt"))%>%
  mutate_each(round)%>%
  as.matrix()
car_age_props<-Carc_adult_age/matrix(rowSums(Carc_adult_age),nrow=nrow(Carc_adult_age),ncol=ncol(Carc_adult_age))
plot(car_age_props[,1],type="o",ylim=c(0,1))
points(car_age_props[,2],type="o",col="blue")
points(car_age_props[,3],type="o",col="red")


alr_Carc_adult_age<-cbind(log((Carc_adult_age[,1]+.01/rowSums(Carc_adult_age)+.01)/
                     (Carc_adult_age[,3]+.01/rowSums(Carc_adult_age)+.01)),
                      log((Carc_adult_age[,2]+.01/rowSums(Carc_adult_age)+.01)/
                         (Carc_adult_age[,3]+.01/rowSums(Carc_adult_age)+.01)))
apply(alr_Carc_adult_age[-1,],2,sd)

plot(alr_Carc_adult_age[,1],type="o")#,ylim=c(0,1))
points(alr_Carc_adult_age[,2],type="o",col="blue")



#prop_age inits
prop_ages_init<-c(.125,.725,.15,.06,.65,.29)

alr_prop_ages_init<-log(c(prop_ages_init[1:2]/prop_ages_init[3],prop_ages_init[4:5]/prop_ages_init[6]))

alr_prop_ages_init_mat<-matrix(alr_prop_ages_init,nrow=23,ncol=4,byrow = T)


##data
SR_dat<-list(log_R_obs=log(all_boots_trans)[1:100,],                 # Observed recruits (posterior from juvenile modle)
             brood_year=c(rep(1:22,times=3),0:21),
             LH=rep(0:3,each=22),
             wild_S_obs=as.numeric(wild_spawners),  
             hatch_S_obs = as.numeric(hatch_spawners),  
             brood_rem = wild_brood_rmoveal,
             hatch_carcs=hatch_carcs,
             total_carcs=total_carcs,
             LH_DAYS=(c(36,64,80,92)-mean(surv_dat$Length.mm))/sd(surv_dat$Length.mm),
             # Observed spawners (from redd counts)
             Model=c(1,1,1,2),
             JLH_A=c(0,0,0,1), # juvenile life history ages at emigration.
             PIT_age_BYs=PIT_age_BYs,
             PIT_age_LHs=PIT_age_LHs,
             PIT_ages=PIT_ages,
             Carc_adult_age=Carc_adult_age,
             surv_dat=as.matrix(surv_dat2[,1:2]),
             surv_yr=(surv_dat2$brood_year_guess-1995),
             Nproj=25,
             rem_rate=0.33,
             NORcutoff =500,
             Hmax=200,
             rule=0
             
             )




#matrix(apply(log(all_boots_trans),2,mean),nrow=22)

#matrix(apply(log(all_boots_trans),2,sd),nrow=22)



SR_pars<-list(log_R_hat=matrix(10,nrow=23,ncol=4),      # latent true number of juveniles
              log_R_obs_sd=matrix(-1.4,nrow=23,ncol=4),  # juvenile observation error
              log_S_obs_cv=log(0.05),                         # Spawner observation error
              log_alpha=log(c(700,700,400,60)), #rep(log(60),4),         # intrinsic productivity
              log_R_max=rep(log(50000),4),          # Asymptotic maximum recruitment
              log_d=c(1,1.8,1.4),
              
              log_proc_sigma=rep(-.5,4)/2,  # process error standard deviation
              logit_proc_er_corr=rnorm(6,0,.03),#,c(.75,.086,.5,.0033,.6,.1),#qlogis((rep(.1,6)/2)+.5),
              
              log_W_ret_init=c(-1,3.7,4.1,4,4.5),
              logit_pHOS=qlogis(phos-.02),
              logit_surv=rep(0,20),
              logit_Phi=.8,
              log_surv_var=-2,
              #logit_surv_cor=qlogis((rep(.98,6))),
              surv_alpha=-5,
              surv_beta=.18,
              alr_p_hyper_mu=alr_prop_ages_init,
              log_alr_p_hyper_sigma=log(rep(.6,4)),
              logit_alr_p_hyper_cor=rnorm(6,0,.03),#c(.7,-.2,.4,-.18,.8,.8),#qlogis(c(-1,1,-1,-1,1,-1)/4+.5),
              alr_p_age=alr_prop_ages_init_mat[1:20,],
              logit_Phi_alr=qlogis(.4)
              #pen_com_surv_log_sigma=exp(1)
              )


#browser()

map<-list(log_S_obs_cv=factor(NA))
map<-list(alr_p_age=factor(rep(NA,length(alr_prop_ages_init_mat))))
log_R_obs_sd_map<-1:(23*4)
log_R_obs_sd_map[c(1,24,47,92)]<-NA
map<-list(log_R_obs_sd=factor(log_R_obs_sd_map),surv_beta=factor(NA))

,log_surv_var=factor(rep(NA,4)))
,log_S_obs_cv=factor(NA))



SR_5<-MakeADFun(SR_dat,SR_pars,random = c("log_R_hat","log_W_ret_init","logit_surv","alr_p_age"),DLL="LCM_lite4",silent = T,map=map)


,"logit_pHOS"



,"logit_alr_p_hyper_cor"
,"logit_proc_er_corr"

,"log_S_obs_cv"
,"logit_surv_cor"




SR_fit<-nlminb(SR_5$par,SR_5$fn,SR_5$gr, 
               control=list(rel.tol=1e-6,eval.max=1000000,
                            iter.max=10000))#,lower=lower,upper=upper)

SR_5$fn()
SR_fit<-nlminb(SR_5$env$last.par.best[-SR_5$env$random],SR_5$fn,SR_5$gr, 
               control=list(rel.tol=1e-12,eval.max=1000000,
                            iter.max=10000))

SR_5$fn()
sd_SR<-sdreport(SR_5) 

#mle and standard errors of s_hats
sd_SR_sum_s_hat_hist<-summary(sd_SR)[(nrow(summary(sd_SR))-(22+23)):(nrow(summary(sd_SR))-23),]
#mle and standard errors of wild_s_hats
wild_sd_SR_sum_s_hat_hist<-summary(sd_SR)[(nrow(summary(sd_SR))-22):nrow(summary(sd_SR)),]

#replace SD of first value with S_obs_CV because its wonky, likely because value is essentially 0. 
wild_sd_SR_sum_s_hat_hist[1,2]<- exp(SR_fit$par["log_S_obs_cv"])


#simulate
result<-function(){

  par(mfcol=c(3,3),mar=c(2,2,2,2),oma=c(2,3.5,0,0))
  #  scenario<-"no hatchery"
  no_hatch<-sim_func(SR_5,"no hatchery",TRUE,2.5)
  
SR_dat$rule<-1
SR_dat$Hmax<-200
SR_dat$NORcutoff<-500
SR_dat$rem_rate<-.33
SR_6<-MakeADFun(SR_dat,SR_pars,random = c("log_R_hat","log_W_ret_init","logit_surv","alr_p_age"),DLL="LCM_lite4",silent = T,map=map)

#scenario<-"Hmax=200, NORcut=500"
rule1<-sim_func(SR_6, "Hmax=200, NORcut=500",FALSE)

# SR_dat$rule<-1
# SR_dat$Hmax<-300
# SR_dat$NORcutoff<-300
# SR_dat$rem_rate<-.33
# SR_7<-MakeADFun(SR_dat,SR_pars,random = c("log_R_hat","log_W_ret_init","logit_surv","alr_p_age"),DLL="LCM_lite4",silent = T,map=map)

#rule2<-sim_func(SR_7, "Hmax=300, NORcut=300")

SR_dat$rule<-1
SR_dat$Hmax<-300
SR_dat$NORcutoff<-300
SR_dat$rem_rate<-.4
SR_7<-MakeADFun(SR_dat,SR_pars,random = c("log_R_hat","log_W_ret_init","logit_surv","alr_p_age"),DLL="LCM_lite4",silent = T,map=map)

rule3<-sim_func(SR_7, "Hmax=600, NORcut=300",FALSE)

cbind(no_hatch,rule1,rule3)

}

sim_func<-function(mod,scenario,labs,line){

set.seed(0114) 
sim_rep<-replicate(1000,{
 
  sim<-mod$simulate(SR_5$env$last.par)
  matrix(c(sim$wild_S_hat,
        sim$hatch_S_hat,
        sim$S_hat,
        sim$wild_return,
        rep(0,23),sim$brood_rem_proj),nrow=(23+25))
})

sim_rep_proj<-sim_rep[24:(23+25),4,]

f4<-rep(1/4,4)
QET<-apply(sim_rep_proj,2,function(x){
  x_lag <- stats::filter(x, f4, sides=1)
  min(x_lag,na.rm = T)<50
})

pQET<-sum(QET)/length(QET)

pQET

geo_mean<-mean(exp(apply(log(sim_rep_proj),2,mean)))

geo_mean


pHOS_mat<-sim_rep[24:(23+25),2,]/sim_rep[24:(23+25),3,]
pHOS_mean<-mean(apply(pHOS_mat,2,mean))
pHOS_mean
pNOB_mat<-sim_rep[24:(23+25),5,]/74
pNOB_mean<-mean(apply(pNOB_mat,2,mean))
pNOB_mean
PNI_mat<-pNOB_mat/(pNOB_mat+pHOS_mat)
pNI_mean<-mean(apply(PNI_mat,2,mean))
pNI_mean
 if(mod$env$data$rule){
plot(c(sim_rep[24:(23+25),4,]),c(PNI_mat),xlim=c(0,400),ylab="pNI",xlab="Wild Return",main=scenario)}else{
  plot(0,type="n",xlim=c(0,400),ylim=c(0,1),main=scenario)
  abline(h=1,lwd=4)
}
if(labs)mtext("pNI",2,line)
segments(0,0,176,0)
segments(176,0,176,.4)
segments(176,.4,207,.5)
segments(208,.5,277,.67)
segments(277,.67,372,.8)
segments(372,.8,1000,1)

out<-c(pQET,geo_mean,pHOS_mean,pNOB_mean,pNI_mean)
#statusQuo<-c(pQET,geo_mean,pHOS_mean,pNOB_mean,pNI_mean)
#dec_tab<-cbind(no_hatch,statusQuo)


s_hat_proj<-sim_rep[24:(23+25),3,]
wild_s_hat_proj<-sim_rep[24:(23+25),1,]

plot_all<-function(sim_out,sd_out,sp_dat,origin){
  
  s_hat_proj_quant<-apply(sim_out,1,
                          quantile,c(.025,.25,.5,.75,.975))
  
  s_hat_hist_quant<-exp(matrix(sd_out[,1],nrow=5,ncol=23,byrow = T)+((matrix(sd_out[,2],nrow=5,ncol=23,byrow = T))*qnorm(c(.025,.25,.5,.75,.975))))
  
  
  s_hat_all<-cbind(s_hat_hist_quant,s_hat_proj_quant)
  
  years<-1995:(2017+25)
  plot(0,type="n",ylim=range(s_hat_all),xlim=range(years),xlab="",ylab="")
  if(labs)mtext(paste(origin,"spawners"),2,line)
  polygon(c(years,rev(years)),c(s_hat_all[1,],rev(s_hat_all[5,])),border=F,col="lightgrey")
  
  polygon(c(years,rev(years)),c(s_hat_all[2,],rev(s_hat_all[4,])),border=F,col="darkgrey")
  
  points(years,s_hat_all[3,],type="l")
  
  points(years[1:23],sp_dat)
  
}# end of function

plot_all(wild_s_hat_proj,wild_sd_SR_sum_s_hat_hist,wild_spawners,"wild\n")

plot_all(s_hat_proj,sd_SR_sum_s_hat_hist,spawners,"")


out
}



result()



plot_all(wild_s_hat_proj,wild_sd_SR_sum_s_hat_hist,wild_spawners)


plot_all(s_hat_proj,sd_SR_sum_s_hat_hist,spawners)

exp(sd_SR_sum_s_hat_hist[,1]+1.96*sd_SR_sum_s_hat_hist[,2])



wild_s_hat_proj<-sim_rep[24:(23+25),1,]
wild_s_hat_proj_quant<-apply(sim_rep_sub,1,
                        quantile,c(.025,.25,.5,.75,.975))




plot(exp(sd_SR_sum_s_hat_hist[,1]+1.96*sd_SR_sum_s_hat_hist[,2]),type="l")
points(exp(sd_SR_sum_s_hat_hist[,1]-1.96*sd_SR_sum_s_hat_hist[,2]),type="l")
points(exp(sd_SR_sum_s_hat_hist[,1]),type="l")
points(spawners[],col="blue",pch=19)



plot(exp(wild_sd_SR_sum_s_hat_hist[,1]+1.96*wild_sd_SR_sum_s_hat_hist[,2]),type="l",ylim=c(0,1000))
points(exp(wild_sd_SR_sum_s_hat_hist[,1]-1.96*wild_sd_SR_sum_s_hat_hist[,2]),type="l")
points(exp(wild_sd_SR_sum_s_hat_hist[,1]),type="l")
points(wild_spawners[],col="blue",pch=19)



lower_s_hat<-apply(sim_rep_sub,1,quantile,.025)

mid_s_hat<-apply(sim_rep_sub,1,quantile,.5)


test$S_hat[1:23]
test$S_obs_cv
test$pHOS


plot(0,type="n",ylim=c(0,max))



sim_rep_sub<-sim_rep[,5,]


upper_s_hat<-apply(sim_rep_sub,1,quantile,.975)

lower_s_hat<-apply(sim_rep_sub,1,quantile,.025)

mid_s_hat<-apply(sim_rep_sub,1,quantile,.5)

plot(0,type = "n",ylim=range(c(upper_s_hat,lower_s_hat)),xlim=c(0,48))


points(upper_s_hat,type="l")
points(lower_s_hat,type="l")
points(mid_s_hat,type="l")









plot(sim$S_hat,type="o")
plot(sim$wild_S_hat,type="o")
plot(sim$surv_proj,type="o")
plot(sim$S_hat,type="o")
plot(sim$wild_S_hat,type="o")








abline(v=23)
sim$prop_age_proj
sim$juv_proc_er_proj
sim$w_ad_age

SR_5$fn()
range(SR_5$gr())
test2<-test
test3<-test
test<-SR_5$report()
test$Recruit_obs_like
test$Spawner_obs_like
test$state_Like
test$carcass_age_like
test$PIT_age_like
test$like_surv
test$like_rand_age
test$phos_like
test$obj_fun
test$alpha
test$wild_S_hat
test$hatch_S_hat
test$wild_return
test$pHOS
test$alr_p_hyper_mu
hist(test$logit_surv)
test$R_max
test$proc_sigma
test$corr_rec_er_mat
test$alr_p_hyper_cor
test$p_cor_mat
test$a
alr_prop_ages_init
hist(plogis(test$logit_surv))
test$S_hat
plot()
test$hatch_S_hat
test$pHOS
test$R_max

library(ggplot2)

p_age<-test$prop_age

plot(p_age[,1],ylim=c(0,1),type="l",col="blue",pch=19,xaxt="n",xlab="",ylab=expression(p))
points(p_age[,2],type="l",col="red")
points(p_age[,3],type="l")
rowSums(p_age[,1:3])
#plot(p_age[,4],ylim=c(0,1),type="o",col="blue")
points(p_age[,4],type="l",lty=2,col="blue")
points(p_age[,5],type="l",col="red",lty=2)
points(p_age[,6],type="l",lty=2)
axis(1,at=seq(1,20,5),labels=seq(1995,2010,5))

legend(y=1.75,x=1,col=rep(c("blue","black","red"),each=2),legend=c("1-3","2-3","1-4","2-4","1-5","3-5"),xpd=NA,ncol=3,lty=rep(1:2,each=3))
rowSums(p_age[,4:6])

sr_out<-data.frame(test$S_hat[1:23],exp(test$log_R_hat))
colnames(sr_out)<-c("Spawners","fry","summer","fall","smolts")

sr_out<-gather(sr_out,"lifestage","Juveniles",2:5)

ggplot(sr_out,aes(Spawners,Juveniles))+geom_point()+facet_wrap(~lifestage)

LH_lab<-c("fry","summer","fall","smolts")

ord<-order(test$S_hat[1:23])
par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,0,0))
for ( i in 1:3){
plot(test$S_hat[1:23],exp(test$log_R_hat[,i]),main=LH_lab[i],xlab="",ylab="",ylim=c(0,200000),pch=19)


polygon(c(test$S_hat[ord],rev(test$S_hat[ord])),
        c(exp(log(test$R_pred[,i][ord])+1.96*test$proc_sigma[i]),
          rev(exp(log(test$R_pred[,i][ord])-1.96*test$proc_sigma[i]))),border = F,col="lightgrey")
points(test$S_hat[ord], test$R_pred[,i][ord],col="blue",type="l",lwd=2)

points(test$S_hat[1:23],exp(test$log_R_hat[,i]),pch=19)
box()

}

plot(test$S_hat[1:23],exp(test$log_R_hat[,4]),main="smolt",xlab="",ylab="",ylim=c(0,200000),pch=19)

i<-4
polygon(c(test$S_hat[ord],rev(test$S_hat[ord])),
        c(exp(log(test$R_pred[,i][ord])+1.96*test$proc_sigma[i]),
          rev(exp(log(test$R_pred[,i][ord])-1.96*test$proc_sigma[i]))),border = F,col="lightgrey")
points(test$S_hat[ord], test$R_pred[,i][ord],col="blue",type="l",lwd=2)

points(test$S_hat[1:23],exp(test$log_R_hat[,i]),pch=19)


points(test$S_hat[ord], test$R_pred[,4][ord],col="blue",type="l",lwd=2)
mtext("Spawners",1,0.5,outer=T,xpd=NA)
mtext("Juveniles",2,0.5,outer=T,xpd=NA)

par(mfrow=c(1,1))
plot(1995:2017,test$pHOS,pch=19,type="o")
points(1995:2017,phos,col="blue",type="l")
mtext("pHOS",2,0.5,outer=T,xpd=NA)

plot(hatch_carcs/total_carcs,test$pHOS,main="pHOS")
mtext("pHOS_Obs",1,0.5,outer=T,xpd=NA)
mtext("pHOS_pred",2,0.5,outer=T,xpd=NA)
abline(0,1,col="blue")

plot(exp(colMeans(log(all_boots_trans))),exp(c(test$log_R_hat[-c(1,24,47,92)])),main="Juveniles")
mtext("log(J_Obs)",1,0.5,outer=T,xpd=NA)
mtext("log(J_pred)",2,0.5,outer=T,xpd=NA)
abline(0,1,col="blue")


plot(apply(log(all_boots_trans),2,sd),exp(c(test$log_R_obs_sd[-c(1,24,47,92)])),main="Juveniles")
mtext("log(J_Obs)_SD",1,0.5,outer=T,xpd=NA)
mtext("log(J_pred)_SD",2,0.5,outer=T,xpd=NA)
abline(0,1,col="blue")

plot(log(spawners),log(test$S_hat),main="spawners")
mtext("log(S_Obs)",1,0.5,outer=T,xpd=NA)
mtext("log(S_pred)",2,0.5,outer=T,xpd=NA)
abline(0,1,col="blue")

test$S_obs_cv


abline(0,1)
plot(spawners,test$S_hat,xlab="spawners",ylab="s_hat")
abline(0,1)



test$S_obs_cv


colMeans(test$prop_age)

colSums(PIT_ages[1:11,])/sum(PIT_ages[1:11,])
colSums(PIT_ages[12:23,])/sum(PIT_ages[1:23,])

surv_out<-data.frame(1995:2014,t(plogis(test$logit_surv_lh_yr)))

colnames(surv_out)<-c("year","fry","summer","fall","smolts")
surv_out<-gather(surv_out,"lifestage","survival",2:5)


ggplot(data=surv_out,aes(year,survival,col=lifestage))+geom_line(size=1.5)+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))


plot(1995:(2017+25),test$S_hat,type="l",ylim=c(0,3000))
points(1995:2017,spawners,col="blue",type="l")
mtext("log(S_Obs)",1,0.5,outer=T,xpd=NA)
mtext("Spawners",2,0.5,outer=T,xpd=NA)
legend(x="topleft",legend=c("model","observed"),col=c("black","blue"),lty=1)




