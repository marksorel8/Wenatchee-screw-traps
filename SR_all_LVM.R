#The purpose of this script is to fit models of transition between spawners and juvenile emigrants expressing different juvenile life histories

#libraries
library(here)
library(tidyverse)
library(viridis)
library(TMB)
library(TMBhelper)
library(readxl)



##  Analysis of screw trap data to estimate daily juvenile emigrants ##
refit<-FALSE # whether to refit models rather that loading existing fitted model objects if they exist

#fit/load model
here::i_am("src/ST_all.R")
source(here("src","ST_all.R"))
ST_all<-ST_all_func(refit=refit)


load(here("results","emigrant_estimates May 13 2021.Rdata"))


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
# fit all combinations of different functional forms for different LHPs (4^4 = 256 unique combinations), assuming functional forms are common among streams *within* LHPs. This takes some time (~ 1 hour or more depending on cores), and not all models converge. 

 #models to try all combinations of where each individual stream and life history can have its own functional form    
 mods<-c(1, #depensatory B-H II
         2, #Beverton-holt
         3,#power
         4) #linear
 
 
# matrix of all combinations, where all streams within a given juvenile life histories have same functional form
mod_mat<-expand.grid(fry=mods,sum=mods,fall=mods,spring=mods)
refit=TRUE
 
 if(!refit & length(list.files(here("results"))[substr(list.files(here("results")),1,10)=="dredge_mat"])>0){
   
   last_file_name<-list.files(here("results"))[substr(list.files(here("results")),1,10)=="dredge_mat"][which.max(lubridate::mdy(substr(list.files(here("results"))[substr(list.files(here("results")),1,10)=="dredge_mat"],12,22)))]
   
   print(paste("Loading file:",last_file_name))
   
   load(file=here("results",list.files(here("results"))[substr(list.files(here("results")),1,10)=="dredge_mat"][which.max(lubridate::mdy(substr(list.files(here("results"))[substr(list.files(here("results")),1,10)=="dredge_mat"],12,22)))]))             #load previous presults
 }else{

  library(foreach)
 library(doParallel)
 
 ncores<-detectCores()-2    # number of cores to use
 
 cl <- makeCluster(ncores)  #settup clusters for parallel computing
 registerDoParallel(cl)
 
 
 
 
 
 
 system.time(
   dredge_mat<-foreach(i = iter(1:nrow(mod_mat),chunksize=nrow(mod_mat)/ncores),.packages=c('TMB','tidyverse','here'),.combine='rbind',.multicombine=TRUE,.inorder=FALSE) %dopar%{
     
     # if(dredge_mat2[i,5]==Inf){
     #load TMB model
     setwd(here("src","TMB"))
     TMB::compile("Stock_recruit_LVM.cpp")
     dyn.load(dynlib("Stock_recruit_LVM"))
     
     out<-rep(NA,10)
     out[1:4]<-as.numeric(mod_mat[i,])  # save functional form combination in output array
     #attempt fitting 2 times. Throws errors when doesn't converge for one reason or another
     fit_mod_result<-fit_mod_iter(rep(as.numeric(mod_mat[i,]),3),streams =0:2,LHs=1:4,n_f=1,no_rand_FF = 0, fit_env = 1, fit_attempts=10,additional_attempts = 20)   
      
     
     # save best BIC of model
     try(out[5] <- min(fit_mod_result$BIC_vec,na.rm=TRUE))           # save BIC
     try(out[6] <- length(fit_mod_result$mod$env$last.par.best))# total coefficients
     try(out[7] <- length(fit_mod_result$mod$par))# number of fixed effects
     try(out[8] <- out[6]-out[7])# number of random effects
     return(out)
     # }else{
       # return(dredge_mat2[i,])
     # }
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
# and .csv
write.csv(dredge_mat,here("results",paste0("dredge_mat",substr(date(),4,10),substr(date(),20,25),".csv"))) # write. csv of top 20 table

}

sum(dredge_mat[,5]==Inf) #number of models out of 256 that didn't converged

#best 20 models
top_20_dredge<-as.data.frame(head(dredge_mat[order(dredge_mat[,5]),c(1:5,7:8)],20)) # grab 20 best models and sort by BIC
colnames(top_20_dredge)<-c("Fry","Summer","Fall","Smolts","BIC","Fixed effects","Random effects") # rename columns
top_20_dredge$Delta_BIC<-round(top_20_dredge$BIC-min(top_20_dredge$BIC),3) #calculate delta BIC
for (i in 1:4)top_20_dredge[,i]<-c("BH II","BH","Power","linear")[top_20_dredge[,i]] # give functional forms meaningful names

#write results to csv file
# write.csv(top_20_dredge,here("results","dredge_5-22-2020.csv")) # write. csv of top 20 table

#---------------------------------------------------------------------
##for cross-validation folds
set.seed(626)
rand_ord<-sample(0:195,196)
set.seed(626) 
fold<-
  c(seq(from=1,by=14,length=14),197)
 
 #combinations of priors to try
prior_combs<-expand.grid(lamb_gam=c(2.5,5,10),lamb_j=c(1,5,10),nu_J=log(c(10000,50000,100000)))
 

#grod seach 
library(foreach)
library(doParallel)
ncores<-5
cl <- makeCluster(ncores)  #settup clusters for parallel computing
registerDoParallel(cl)

#loop through parameter sets and fit Astoria departure (migration) timing and store fits
#this takes about three hours
start<-Sys.time()
results<-foreach(
  i=iter(1:nrow(prior_combs),chunksize=nrow(prior_combs)/ncores),
  .packages=c('TMB','tidyverse','here'),.export = c(),.combine='rbind',.inorder=TRUE) %dopar%{
  
  setwd(here("src","TMB"))
  dyn.load(dynlib("Stock_recruit_LVM"))

LL_exclude<-rep(NA,196)
for (i in 1:nrow(prior_combs)){
for(j in 85:196){
  if(is.na(LL_exclude[j])){
fit_mod_result<-fit_mod_iter(rep(c(1,1,1,1),times=3),streams=0:2,LHs=1:4,n_f=1,no_rand_FF =0, fit_env =1, fit_attempts=5,additional_attempts = 4,rate=c(5,5,5,5,1,log(5000)),fold_exclude=(j-1))
                               
        gc()                       # rand_ord[fold[j]:(fold[j+1]-1)])   
try(LL_exclude[j]<-fit_mod_result$report$ll_exclude)
        print(j)
        print(LL_exclude[j])
   }
}
}
return(LL_exclude)
}#end of loop
sum(is.na(results))
end<-Sys.time()
rowSums(results) %>% sort

# write.csv(results,here("results","LL_exclude_loo.csv"))
results<-read.csv(here("results","LL_exclude_loo.csv"))[,-1]
rowSums(results) %>% `names<-`(1:27) %>%  sort()
prior_combs2<- prior_combs %>% mutate(LL_exc=rowSums(results))

prior_combs2 %>% arrange(LL_exc)
write.csv(prior_combs2,here("results","prior_combs2.csv"))




stopCluster(cl) #shut down parallel computing clusters
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister()


results2<-array(unlist(results),dim=c(dim(results[[1]]),Niter))
return(results2)
}

 for(j in 1:nrow(prior_combs)){
 
 }


 
thing2<- sum(ll_fold_exclude2)
sum(ll_fold_exclude2)
sum(ll_fold_exclude3)
sum(ll_fold_exclude4)

test_sim<-rmvnorm_prec(mu=fit_mod_result$mod$env$last.par.best ,
                       prec=fit_mod_result$fit$SD$jointPrecision ,
                       n=500,random_seed=500)


test_x<-apply(test_sim,2,function(x){fit_mod_result$mod$report(x)[["juv_like"]]})

loo::loo(t(test_x))


fit_mod_result<-fit_mod_iter(rep(c(1,1,1,1),times=3),streams=0:2,LHs=1:4,n_f=1,no_rand_FF =0, fit_env =1, fit_attempts=5,additional_attempts = 2,rate=c(5,5,5,5,1,log(5000)),fold_exclude=integer(0)) 

lower<-rep(-Inf,length(fit_mod_result$mod$env$last.par.best))
lower[which(names(fit_mod_result$mod$env$last.par.best)=="Loadings_vec")[5]]<-0
stand_test<-tmbstan::tmbstan(fit_mod_result$mod, cores=4,lower = lower,upper=rep(Inf,length(lower)),
                             init=fit_mod_result$mod$env$last.par.best)
                             
shinystan::launch_shinystan(stand_test)


thing<-as.data.frame(stand_test)
colnames(thing)

sd(exp(thing$`beta_alpha[1]`+thing$`eps_alpha[2]`*exp(thing$`log_sigma_alpha[1]`)))

fit_mod_result$BIC_vec # BICs of fit attempts

# look at some interesting model values
  matrix(fit_mod_result$report$alpha,4)  # alphas
  matrix(fit_mod_result$report$gamma,4)  # gammas (note, some are just the initial values if not included in a model)
  matrix(fit_mod_result$report$Jmax,4)   # Jmaxes (also some not iuncluded so jsut intial values)
  
  fit_mod_result$report$log_sigma_alpha      # alpha random effect standard deviations
  fit_mod_result$report$log_sigma_gamma      # gamma random effect standard deviations (not all meaningfl, e.g. no gamma in B-H for falls)
  fit_mod_result$report$log_sigma_Jmax       # Jmax random effect standard deviations (only meanful for summers)
  
  
 FF_params<- data.frame(stream=rep(c("Chiwawa","Nason","White"),each=4),
  LH=rep(c("Spr-0","Sum-0","Fall-0", "Spr-1"),times=3),           
  alpha=fit_mod_result$report$alpha,
  alpha_se=fit_mod_result$fit$SD$sd[names(fit_mod_result$fit$SD$value)=="alpha"],
  gamma=fit_mod_result$report$gamma,
  gamm_se=fit_mod_result$fit$SD$sd[names(fit_mod_result$fit$SD$value)=="gamma"],
  Jmax=fit_mod_result$report$Jmax,
  Jmax_se=fit_mod_result$fit$SD$sd[names(fit_mod_result$fit$SD$value)=="Jmax"])
 FF_params<-FF_params[rep(c(1,5,9),times=4)+rep(0:3,each=3),]
 FF_params[,c(-1:-2)]<-round(FF_params[,c(-1:-2)],2)
 write.csv(FF_params,here("results","FF_params.csv"))
 
 

  # plot of latent and expected juveniles vs spawners 
  ##png(file=here("results","plots","spawn_em_625020.png"),units="in",height=5,width=6.5,res=300)
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
   fit_mod_result<-fit_mod_iter(rep(c(4,4,2,3),times=3),streams=0:2,LHs=1:4,n_f=1,no_rand_FF =0, fit_env = 1, fit_attempts=10,additional_attempts = 2)   
   
   fit_mod_result$BIC_vec # BICs of fit attempts
   
# p-values of environmental covariate coefficients
pnorm(0,abs(fit_mod_result$fit$SD$par.fixed),sqrt(diag(fit_mod_result$fit$SD$cov.fixed)),lower.tail = TRUE)[names(fit_mod_result$fit$SD$par.fixed)=="beta_e"]*2


# plot environmental covariate coefficients
# png(here("results","plots","coef_plot_all.png"),units="in",res=300,height=4,width=5)
bootstrap_env_cov(fit_mod_result$dat,last_best=fit_mod_result$fit$par , precis=fit_mod_result$fit$SD$cov.fixed )
# dev.off()  

test<-fit_mod_result$fit$par[names(fit_mod_result$fit$par)=="beta_e"]
test<-fit_mod_result$mod$env$last.par.best[names(fit_mod_result$mod$env$last.par.best)=="beta_e"]
names(test)<-fit_mod_result$dat$X %>% colnames()
test

test5<-cbind(test,sqrt(diag(fit_mod_result$fit$SD$cov.fixed))[names(fit_mod_result$fit$par)=="beta_e"])  




test<-tibble(stream=c(rep(c("Chiwawa","Nason","White"),each=4),
                      rep(c("Chiwawa","Nason","White"),each=3),
                      rep(c("Chiwawa","Nason","White"),each=1)
                      # rep(c("Chiwawa","Nason","White"),each=4)
                      ),
                      
                        
             LH=factor(c(rep(c("Spr. Sub","Sum. Sub","Fall Sub", "Spr. Yrl."),times=3),
                  rep(c("Sum. Sub","Fall Sub", "Spr. Yrl."),times=3),
                  rep(c( "Spr. Yrl."),times=3)
                  # rep(c("Spr. Sub","Sum. Sub","Fall Sub", "Spr. Yrl."),times=3)
                  ),levels=c(c("Spr. Sub","Sum. Sub","Fall Sub", "Spr. Yrl."))
                  ),
               
              
             cov= c(substr(colnames(fit_mod_result$dat$X)[1:12],1,5),
             substr(colnames(fit_mod_result$dat$X)[13:24],nchar(colnames(fit_mod_result$dat$X)[13:24])-4,nchar(colnames(fit_mod_result$dat$X)[13:24]))),
             mu=fit_mod_result$fit$SD$value[names(fit_mod_result$fit$SD$value)=="beta_e"],
se=fit_mod_result$fit$SD$sd[names(fit_mod_result$fit$SD$value)=="beta_e"],
par=colnames(fit_mod_result$dat$X))

test2<-test %>% filter(cov!=":pHOS") %>% mutate(cov=factor(cov,levels=c("win_0","sum_0","win_1")))

ggplot(test2, aes(x=LH,y=mu,color=stream))+ geom_hline(yintercept=0,linetype=2)+geom_point(position=position_dodge(width = .5))+facet_grid(~cov,scale="free_x", space = "free_x")+ geom_linerange(aes(ymin=mu+1.96*se,ymax=mu-1.96*se),position=position_dodge(width = .5))

test3<-test %>% filter(cov==":pHOS") 

ggplot(test3, aes(x=LH,y=mu,color=stream))+ geom_hline(yintercept=0,linetype=2)+geom_point(position=position_dodge(width = .5))+facet_wrap(~cov)+ geom_linerange(aes(ymin=mu+1.96*se,ymax=mu-1.96*se),position=position_dodge(width = .5))


#----------------------------------------------------------------------------
test<-tibble(
             
             
             LH=factor(c(rep(c("Spr. Sub","Sum. Sub","Fall Sub", "Spr. Yrl."),times=1),
                         rep(c("Sum. Sub","Fall Sub", "Spr. Yrl."),times=1),
                         rep(c( "Spr. Yrl."),times=1),
                         rep(c("Spr. Sub","Sum. Sub","Fall Sub", "Spr. Yrl."),times=1)),levels=c(c("Spr. Sub","Sum. Sub","Fall Sub", "Spr. Yrl."))),
             
             
             cov= c(substr(colnames(fit_mod_result$dat$X)[1:4],1,5),
                    substr(colnames(fit_mod_result$dat$X)[5:12],nchar(colnames(fit_mod_result$dat$X)[5:12])-4,nchar(colnames(fit_mod_result$dat$X)[5:12]))),
             mu=fit_mod_result$fit$SD$value[names(fit_mod_result$fit$SD$value)=="beta_e"],
             se=fit_mod_result$fit$SD$sd[names(fit_mod_result$fit$SD$value)=="beta_e"],
             par=colnames(fit_mod_result$dat$X))

test2<-test %>% filter(cov!=":pHOS") %>% mutate(cov=factor(cov,levels=c("win_0","sum_0","win_1")))

ggplot(test2, aes(x=LH,y=mu,))+ geom_hline(yintercept=0,linetype=2)+geom_point(position=position_dodge(width = .5))+facet_grid(~cov,scale="free_x", space = "free_x")+ geom_linerange(aes(ymin=mu+1.96*se,ymax=mu-1.96*se),position=position_dodge(width = .5))

test3<-test %>% filter(cov==":pHOS") 

ggplot(test3, aes(x=LH,y=mu))+ geom_hline(yintercept=0,linetype=2)+geom_point(position=position_dodge(width = .5))+facet_wrap(~cov)+ geom_linerange(aes(ymin=mu+1.96*se,ymax=mu-1.96*se),position=position_dodge(width = .5))



#----------------------------------------------------------------------------



#--------------------------------------------------------------------- 
# fit each stream x LHP combination individually (exploratory)


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
##png(file=here("results","plots","spawn_em_ind_16020.png"),units="in",height=5,width=6.5,res=300)
# dev.new() 
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
#---------------------------------------------------------------------
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
