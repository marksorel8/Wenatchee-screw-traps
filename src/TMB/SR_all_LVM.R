#The purpose of this script is to fit models of transition between spawners and juvenile emigrants expressing different juvenile life histories

#libraries
library(here)
library(tidyverse)
library(viridis)
library(TMB)
library(TMBhelper)
library(readxl)



##  Analysis of screw trap data to estimate daily juvenile emigrants ##

#fit/load model
source(here("src","ST_all.R"))
ST_all<-ST_all_func()

# load functions to calculate geometric means of dialy emigrants and plot
source(here("src","ST_plotting_funcs.R"))

#plot average daily emigrants and LHP breaks
# png(file=here("results","plots","temp_dis.png"),units="in",width = 6,height=4,res=300)
with(ST_all,ggplot_timing_func(all_emigrants_estimates,all_data_lists,across_stream_geomean,mix_geomean_dens,breaks))
# dev.off()


#plot average daily temperature and discharge
# png(file=here("results","plots","temp_dis.png"),units="in",width = 4.5,height=3,res=300)
with(ST_all,plot_dis_temp_func(all_data_lists))
# dev.off()

rm(list=ls()[which(ls()!="ST_all")])#remove all all teh functions used for screw trap model


#----------------------------------------------------------------------------
##  Analysis of screw trap data to estimate daily juvenile emigrants ##

##   Load data ##
{
  #read redd data from Hatchery program annual report (Hillman et al 2020)
  redds<-read_csv(here("data","redd_counts.csv")) %>% pivot_longer(c("Chiwawa","Nason","White"),"stream",values_to="redds") 
  
  ##emigrant abundance estiamtes (log means  emigra emigrascrew trapabunel
  all_emigrants_estimates<-ST_all$all_emigrants_estimates; rm("ST_all")eriors from screw trap model
  
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
 # fit each stream x LHP combination individually (exploratory)
 
 #models to try all combinations of where each individual stream and life history can have its own functional form    
 mods<-c(1, #depensatory B-H II
         2, #Beverton-holt
         3,#power
         4) #linear


 #empty array to hold BICs
 dredge_indivudal<-array(NA,dim=c(length(mods),3,4),dimnames=list(c("BH II","BH","power","Linear"),c("Chiwawa","Nason","White"),c("fry","summer","fall","smolt")))
 
 
 # note: this will throw errors or warnings about models that don't converge for one reason or another. 
 for ( i in 0:2){ # loop over streams
   
   for( j in 1:4){ # loop over LHs
     
     for (k in 1:length(mods)){ # loop over functional forms
       #attempt fitting 5 times. Throws errors when doesn't converge for one reason or another
       fit_mod_result<-fit_mod_iter(k,i,j,n_f=0,no_rand_FF = 1, fit_env = 0, fit_attempts=2,additional_attempts=10)   
       # save best BIC of model
     dredge_indivudal[k,i+1,j]<-min(fit_mod_result$BIC_vec,na.rm=TRUE)
     }
   }
   }
     
    
 #results table
dredge_indivudal
 #best model for each LHP * stream fit individually
 best_mod<-mods[apply(dredge_indivudal,c(3:2),which.min)]; print(c("BH II","BH","power","Linear")[best_mod])
 
 #write results to CSV file
 #as_tibble(dredge_indivudal) %>% t() %>% `colnames<-`(c("BH II","BH","power","Linear")) %>% as.data.frame() %>% write.csv (file=here("results","individual_BIC.csv"),)
 
 
# fit best models all together
 fit_mod_result<-fit_mod_iter(best_mod,streams=0:2,LHs=1:4,n_f=1,no_rand_FF =1, fit_env = 0, fit_attempts=10,additional_attempts = 10)   
 
 fit_mod_result$BIC_vec # BICs of fit attempts
 
 # plot of latent and expected juveniles vs spawners 
 ##png(file=here("results","plots","spawn_em_ind_14020.png"),units="in",height=5,width=6.5,res=300)
 dev.new() 
 spawn_em<-ggplot_spawner_juveniles(mod_fit = fit_mod_result$fit,mod_dat =fit_mod_result$dat,mod_rep=fit_mod_result$report)
 ## dev.off()
 
 # plot process error correlation (takes a minute for bootstrapping p-values)
 dev.new() 
 corr_plot<-plot_cor(mod_rep = fit_mod_result$report, mod_fit =  fit_mod_result$fit, mod_dat = fit_mod_result$dat)
 
 #spawners and juvenile by year
 plot_spawn_em_year(mod_dat=fit_mod_result$dat,sum_out = spawn_em$sum_out)

 ##expected emigrants vs. spawners
 dev.new() 
 expec_spawn_em_ploft_func(preds=spawn_em$preds)
 
 
 #---------------------------------------------------------------------
 ##fit each LH individually across streams, assuming a common functional form across streams (exploratory)

  ##empty array to hold BIC results
 dredge_ff<-matrix(NA,length(mods),4,dimnames=list(c("BH II","BH","power","Linear"),c("fry","summer","fall","smolt")))
 
 # note: this will throw lots of errors or warnings about models that don't converge for one reason or another. 
 for ( k in 1:length(mods)){ # loop over models
   for (j in 1:4){           # loop over life histories
    
     #attempt fitting 10 times. Throws errors when doesn't converge for one reason or another
     fit_mod_result<-fit_mod_iter(rep(k,times=3),streams =0:2,LHs=j,n_f=0,no_rand_FF = 0, fit_env = 0, fit_attempts=30, additional_attempts = 20)   
     # save best BIC of model
     dredge_ff[k,j]<-min(fit_mod_result$BIC_vec,na.rm=TRUE)
     
      }
   print(k) # print iteration in loop over model functional forms after fitting for all stream and LHs to track progress. 
 }
 
 #results table
dredge_ff
#write results to csv file
#write.csv(t(dredge_ff),file=here("results","ind_LHP_dredge.csv"))
 
#best model for each LHP * stream fit individually
 best_mod<-apply(dredge_ff,2,which.min); print(c("BH II","BH","power","Linear")[best_mod])
 #---------------------------------------------------------------------
# fit all combinations of different functional forms for different LHPs (4^4 = 256 unique combinations), assuming functional forms are common among streams *within* LHPs. This takes some time (~ 1 hour or more depending on cores), and not all models converge. 

# matrix of all combinations, where all streams within a given juvenile life histories have same functional form
mod_mat<-expand.grid(fry=mods,sum=mods,fall=mods,spring=mods)

 if(length(list.files(here("results"))[substr(list.files(here("results")),1,18)=="dredge_mat"])>0){
   
   load(file=here("results",list.files(here("results"))[substr(list.files(here("results")),1,10)=="dredge_mat"][which.max(lubridate::mdy(substr(list.files(here("results"))[substr(list.files(here("results")),1,10)=="dredge_mat"],12,22)))]))             #load previous presults
 }else{

  library(foreach)
 library(doParallel)
 
 ncores<-detectCores()-2    # number of cores to use
 
 cl <- makeCluster(ncores)  #settup clusters for parallel computing
 registerDoParallel(cl)
 
 system.time(
   dredge_mat2<-foreach(i = iter(1:nrow(mod_mat),chunksize=nrow(mod_mat)/ncores),.packages=c('TMB','tidyverse','here'),.combine='rbind',.multicombine=TRUE,.inorder=FALSE) %dopar%{
     
     
     #load TMB model
     setwd(here("src","TMB"))
     TMB::compile("Stock_recruit_LVM.cpp")
     dyn.load(dynlib("Stock_recruit_LVM"))
     
     out<-rep(NA,10)
     out[1:4]<-as.numeric(mod_mat[i,])  # save functional form combination in output array
     #attempt fitting 2 times. Throws errors when doesn't converge for one reason or another
     fit_mod_result<-fit_mod_iter(rep(as.numeric(mod_mat[i,]),3),streams =0:2,LHs=1:4,n_f=1,no_rand_FF = 0, fit_env = 0, fit_attempts=10,additional_attempts = 20)   
      
     
     # save best BIC of model
     try(out[5]<-min(fit_mod_result$BIC_vec,na.rm=TRUE))           # save BIC
     try(out[6]<- length(fit_mod_result$mod$env$last.par.best))# total coefficients
     try(out[7]<-length(fit_mod_result$mod$par))# number of fixed effects
     try(out[8]<- out[6]-out[7])# number of random effects
     return(out)
   }#end of loop
   
 )
 
  
 
  stopCluster(cl) #shut down parallel computing clusters
 unregister <- function() {
   env <- foreach:::.foreachGlobals
   rm(list=ls(name=env), pos=env)
 }
 unregister()
 
 
 #save result table as R object 
save(dredge_mat,file=here("results",paste0("dredge_mat",substr(date(),4,10),substr(date(),20,25),".Rdata")))  #save results
}



sum(dredge_mat[,5]==Inf) #number of models out of 256 that didn't converged

#best 20 models
top_20_dredge<-as.data.frame(head(dredge_mat2[order(dredge_mat2[,5]),c(1:5,7:8)],20)) # grab 20 best models and sort by BIC
colnames(top_20_dredge)<-c("Fry","Summer","Fall","Smolts","BIC","Fixed effects","Random effects") # rename columns
top_20_dredge$Delta_BIC<-round(top_20_dredge$BIC-min(top_20_dredge$BIC),3) #calculate delta BIC
for (i in 1:4)top_20_dredge[,i]<-c("BH II","BH","Power","linear")[top_20_dredge[,i]] # give functional forms meaningful names

#write results to csv file
# write.csv(top_20_dredge,here("results","dredge_1-5-2020.csv")) # write. csv of top 20 table

#---------------------------------------------------------------------

# fit best combinations of functional forms, where all streams have common function form within a given LHP
# functions forms: fry, summer, smolt = power; fall = Beverton-Holt
# fit best models all together
fit_mod_result<-fit_mod_iter(rep(c(3,3,2,3),times=3),streams=0:2,LHs=1:4,n_f=1,no_rand_FF =0, fit_env = 0, fit_attempts=10,additional_attempts = 10)   

fit_mod_result$BIC_vec # BICs of fit attempts

# look at some interesting model values
  matrix(fit_mod_result$report$alpha,4)  # alphas
  matrix(fit_mod_result$report$gamma,4)  # gammas (note, some are just the initial values if not included in a model)
  matrix(fit_mod_result$report$Jmax,4)   # Jmaxes (also some not iuncluded so jsut intial values)
  
  fit_mod_result$report$log_sigma_alpha      # alpha random effect standard deviations
  fit_mod_result$report$log_sigma_gamma      # gamma random effect standard deviations (not all meaningfl, e.g. no gamma in B-H for falls)
  fit_mod_result$report$log_sigma_Jmax       # Jmax random effect standard deviations (only meanful for summers)
    
  
  
  # plot of latent and expected juveniles vs spawners 
  ##png(file=here("results","plots","spawn_em_123020.png"),units="in",height=5,width=6.5,res=300)
  spawn_em<-ggplot_spawner_juveniles(mod_fit = fit_mod_result$fit,mod_dat =fit_mod_result$dat,mod_rep=fit_mod_result$report)
  ##dev.off()

  # plot process error correlation (takes a minute for bootstrapping p-values)
  ## png(here("results","plots","correlation.png"),units="in",res=300,height=10,width=10)
  corr_plot<-plot_cor(mod_rep = fit_mod_result$report, mod_fit =  fit_mod_result$fit, mod_dat = fit_mod_result$dat)
  ## dev.off()  
 
#spawners and juvenile by year
# png(here("results","plots","spaw_em.png"),units="in",res=300,height=5,width=7)
   plot_spawn_em_year(mod_dat=fit_mod_result$dat,sum_out = spawn_em$sum_out)
   # dev.off() 


##expected emigrants vs. spawners
# png(here("results","plots","juv_LH.png"),units="in",res=300,height=4,width=5)
   expec_spawn_em_ploft_func(preds=spawn_em$preds)
  # dev.off()
  

## fit  model with environmental covariates
   fit_mod_result<-fit_mod_iter(rep(c(3,3,2,3),times=3),streams=0:2,LHs=1:4,n_f=1,no_rand_FF =0, fit_env = 1, fit_attempts=10,additional_attempts = 10)   
   
   fit_mod_result$BIC_vec # BICs of fit attempts
   
# p-values of environmental covariate coefficients
pnorm(0,abs(fit_mod_result$fit$SD$par.fixed),sqrt(diag(fit_mod_result$fit$SD$cov.fixed)),lower.tail = TRUE)[names(fit_mod_result$fit$SD$par.fixed)=="beta_e"]*2


# plot environmental covariate coefficients
# png(here("results","plots","coef_plot_all.png"),units="in",res=300,height=4,width=5)
bootstrap_env_cov(fit_mod_result$dat,last_best=fit_mod_result$fit$par , precis=fit_mod_result$fit$SD$cov.fixed )
# dev.off()  

