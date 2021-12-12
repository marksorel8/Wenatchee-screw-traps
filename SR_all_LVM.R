#The purpose of this script is to fit models of transition between spawners and juvenile emigrants expressing different juvenile life histories

#libraries
library(here)
library(tidyverse)
library(viridis)
library(TMB)
library(TMBhelper)
library(readxl)
library(glmmTMB)


##  Analysis of screw trap data to estimate daily juvenile emigrants ##
refit<-FALSE # whether to refit models rather that loading existing fitted model objects if they exist (takes about ~20 minutes; 5 for model fitting and 15 for parametric bootstrap)

#fit/load model
source(here("src","ST_all.R"))
ST_all<-ST_all_func(refit=refit)


# load functions to calculate geometric means of dialy emigrants and plot
source(here("src","ST_plotting_funcs.R"))

#plot average daily emigrants and LHP breaks
# png(file=here("results","plots","emigration_timing.png"),units="in",width = 6,height=4,res=300)
with(ST_all,ggplot_timing_func(all_emigrants_estimates,all_data_lists,across_stream_geomean,mix_geomean_dens,breaks))
# dev.off()


#plot average daily temperature and discharge
# png(file=here("results","plots","temp_dis.png"),units="in",width = 4.5,height=3,res=300)
with(ST_all,plot_dis_temp_func(all_data_lists))
# dev.off()

rm(list=ls()[which(ls()!="ST_all")])#remove all all the functions used for screw trap model


#----------------------------------------------------------------------------
##  Analysis of screw trap data to estimate daily juvenile emigrants ##

##   Load data ##
{
  #read redd data from Hatchery program annual report (Hillman et al 2020)
  redds<-read_csv(here("data","redd_counts.csv")) %>% pivot_longer(c("Chiwawa","Nason","White"),"stream",values_to="redds") 
  
  pHOS<-read_excel(here("data","adult data.xlsx"),2)
  
  test<-full_join(redds,pHOS %>% rename("Year"="Brood_year","stream"="Stream"))
  
  summary(lm(redds~stream+pHOS_weighted:stream-1,data = test))
  
  ##emigrant abundance estimates (log means and standard deviations from screw trap model)
  all_emigrants_estimates<-ST_all$all_emigrants_estimates; rm("ST_all")
  
  #load stream flow covariate values (flow_covs)
  source(here("src","covariates.r"))
  
  # length of habitat surveyed for spawners in each tributary (based on motoring annual report, Hillman et al. 2020)
  trib_lengths<-c(32.5,15.4,16.1)
  names(trib_lengths)<-c("Chiwawa","Nason","White")
}



##   Load functions ##

# source helper functions to make data object, make initial parameter values, and makes map (which fixes certain parameter),  all to feed to the TMB model "Stock_recruit_LVM"
source(here("src","SR_helper_functions.R"))
 
 # source plotting functions for visualizing results.
 source(here("src","SR_plotting_funcs_factor.R"))
 
 #--------------------------------------------------------------------- 
 #--------------------------------------------------------------------- 
 ##   Fit models  ##
 
 #load TMB model
 setwd(here("src","TMB"))
 TMB::compile("Stock_recruit_LVM.cpp")
 dyn.load(dynlib("Stock_recruit_LVM"))
 
#---------------------------------------------------------------------

#sometimes have to try a number of times with different starting calues to get a model to converge (takes several minutes)
#there should ideally be multiple converged model with the same BIC = 612.9517
 ## if its not working, I recommend closing r and then trying again
 set.seed(1234)
fit_mod_result<-fit_mod_iter(streams=0:2,     #streams to include (Chiwawa, Nason, White)
                             LHs=1:4,         #life histories to include (spring subs, sumemr subs, fall subs, spring yearlings)
                             n_f=1,           # number of latent variable factors to include
                             fit_env =1,      #whether to fit environmentel covariates (1=year, 0 = no)
                             fit_attempts=100, # numebr of times to attempt to fit model
                             additional_attempts = 0, #additional attempts after that until a model converges
                             log_J_max_prior=log(15000))   # prior mean on Jmax hyper-mean
 
#BIC value from each iteration running the model with different initial parameters
##There should ideally be multiple converged model with the same BIC = 612.9517
fit_mod_result$BIC_vec                      
# AIC of best fit
fit_mod_result$fit$AIC

save(fit_mod_result,file=here("results","fit_mod_result_12_10_21.Rdata"))

#list of final parameter values
par_out<-fit_mod_result$mod$env$parList(par=fit_mod_result$mod$env$last.par.best)
par_out$rate %>% exp() #penalty parameters
gc() #garbage clean
  

  #parametric bootstrap from posterior of fitted model
  sim_post<-rmvnorm_prec(fit_mod_result$mod$env$last.par.best, 
               fit_mod_result$fit$SD$jointPrecision, 50000, 623 )
    

  
  #reporting
  
  ##calculate and store functional form parameters (alpha, gamma, Jmax) for each posterior samples
  FFparams<-array(dim=c(12,50000,3))
  for (i in 1:50000){
    out<-fit_mod_result$mod$report(par=sim_post[,i])
    FFparams[,i,1]<-out$alpha
    FFparams[,i,2]<-out$gamma
    FFparams[,i,3]<-out$Jmax
  }
  ##calculate quantiles of FF parameters 
  quant_FF<-apply(FFparams,c(1,3),quantile,probs=c(0.025,.5,0.975))
  
 ## table of  functional relationship fit parameters
   FF_params<- tibble(stream=rep(c("Chiwawa","Nason","White"),each=4),
  LH=factor(rep(c("Spr-0","Sum-0","Fall-0", "Spr-1"),times=3),levels=c("Spr-0","Sum-0","Fall-0", "Spr-1")),           
  alpha=fit_mod_result$report$alpha,

  alpha_lcl=quant_FF[1,,1],
  alpha_ucl=quant_FF[3,,1],
  gamma=fit_mod_result$report$gamma,

  gamma_lcl=quant_FF[1,,2],
  gamma_ucl=quant_FF[3,,2],
  Jmax=fit_mod_result$report$Jmax,

  Jmax_lcl=quant_FF[1,,3],
  Jmax_ucl=quant_FF[3,,3],) %>% 

    arrange(LH,stream) %>% 
    mutate(across(3:11,round,2))
  
  View(FF_params)

 write.csv(FF_params,here("results","FF_params.csv"))
 
 
 #hyper means of parameters
 ##gamma
 exp(par_out$beta_gamma)
 exp(
 sim_post[
 which(names(fit_mod_result$mod$env$last.par.best)=="beta_gamma"),]) %>% 
   apply(1,quantile,probs=c(.025,.5,.975))
 ##Jmax
 exp(par_out$beta_Jmax)
 exp(
   sim_post[
     which(names(fit_mod_result$mod$env$last.par.best)=="beta_Jmax"),]) %>% 
   apply(1,quantile,probs=c(.025,.5,.975))
 
 
 

  # plot of latent and expected juveniles vs spawners 
  ##png(file=here("results","plots","spawn_em_11112021.png"),units="in",height=5,width=6.5,res=300)
  spawn_em<-ggplot_spawner_juveniles(mod_fit = fit_mod_result$fit,mod_dat =fit_mod_result$dat,mod_rep=fit_mod_result$report)
  ##dev.off()

  
#spawners and juvenile by year
# png(here("results","plots","spaw_em.png"),units="in",res=300,height=5,width=7)
   plot_spawn_em_year(mod_dat=fit_mod_result$dat,sum_out = spawn_em$sum_out)
   # dev.off() 


##expected emigrants vs. spawners
# png(here("results","plots","juv_LH.png"),units="in",res=300,height=4,width=5)
   expec_spawn_em_ploft_func(preds=spawn_em$preds)
  # dev.off()
  


# Probability environmental covariate coefficients are not equal to zero
pnorm(0,abs(fit_mod_result$fit$SD$par.fixed),sqrt(diag(fit_mod_result$fit$SD$cov.fixed)),lower.tail = TRUE)[names(fit_mod_result$fit$SD$par.fixed)=="beta_e"]*2


# plot environmental covariate coefficients
# png(here("results","plots","coef_plot_all.png"),units="in",res=300,height=4,width=5)
bootstrap_env_cov(fit_mod_result$dat,last_best=fit_mod_result$fit$par , precis=fit_mod_result$fit$SD$cov.fixed )
# dev.off()  

# quantiles of paramateric bootstrap distribution of environmental covariates
fit_mod_result$mod$env$last.par.best[which(names(fit_mod_result$mod$env$last.par.best) =="beta_e")] 
sim_post[which(names(fit_mod_result$mod$env$last.par.best) =="beta_e"),] %>% apply(1,quantile,probs=c(.025,.5,.975))




#Refit model with no environmental covariates to look at process erros covariance (&correlation) based on the loadings onto the latent variable


#sometimes have to run this a few times to get a model to converge
fit_mod_result_no_env<-fit_mod_iter(streams=0:2,LHs=1:4,n_f=1,fit_env =0, fit_attempts=50,additional_attempts = 0,log_J_max_prior=log(15000))   

#BIC value from each iteration running the model with different initial parameters
##There should ideally be multiple converged model with the same BIC = 606.946
fit_mod_result_no_env$BIC_vec                      
# AIC of best fit
fit_mod_result_no_env$fit$AIC

save(fit_mod_result_no_env,file=here("results","fit_mod_result_no_env_12_10_21.Rdata"))

# plot process error correlation (takes a minute for bootstrapping p-values)
## png(here("results","plots","correlation.png"),units="in",res=300,height=10,width=10)
corr_plot<-plot_cor(mod_rep = fit_mod_result_no_env$report, mod_fit =  fit_mod_result_no_env$fit, mod_dat = fit_mod_result_no_env$dat)
## dev.off()  

#end of analysis