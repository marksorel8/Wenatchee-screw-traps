# Mark Sorel

# script estimates the number of juvenile emigrants from natal streams in the Wenatchee River basin 

#takes data collected from screw traps provided by Washington Department of Fish and Wildlife (Chiwawa River) and Yakama Nation Fisheries (Nason and White Rivers) and returns estimates of daily emigrants and log sums of emigrants within migrations windows representing juvenile life histories (Fry, summer parr, fall parr, yearling smolt), and standard errors. Also plots the time series of geometric mean number of migrants on each day of year across years.  

library(TMB)
library(here)
library(tidyverse)

#------------------------------------------------------------------------------------
# calling functions to conduct analysis. Function are below and must be loaded first.
#------------------------------------------------------------------------------------

#process data
chiw_data<-Chiw_dat_Proc() #chiwawa
nas_whi_data<-Nason_White_data_Func() #nason and white

#make data lists for all models
all_data_lists<-make_data_list_func(chiw_data = chiw_data, nas_whi_data = nas_whi_data)

#fit all models
ts<-Sys.time() 
all_emigrants_estimates<-fit_all(all_data_lists=all_data_lists)
Sys.time()-ts

#save as ".Rdata" objects
save(all_emigrants_estimates,file=here("results",paste0("emigrant_estimates",substr(date(),4,10),substr(date(),20,25),".Rdata")))

#plot geomeans of emigrants for each day of year across years
plot_all_geomean_daily_emigrants(all_emigrants_estimates)

#------------------------------------------------------------------------------------
# Function for analysis
#------------------------------------------------------------------------------------

#functions to read and process data
source(here("src","Chiwawa Data Proc 2.R")) # Chiwawa 
source(here("src","Nason White Data Proc 2.R")) # Nason and White



#function to make data lists for model
make_screw_trap_model_data<-function(data_in,stream="Chiwawa", lifestage="yrlng",Use_NB=1,subyearlings=0){
  
   #subset data to days when there was catch data (including 0 catch)
  data_in<-filter(data_in,!is.na(data_in[,paste0("count.",lifestage)]))
  
  #design matrix for trap efficiency model
  #if(stream=="Nason")
  #pdat<-model.matrix(~moved*scale(Disch.cfs) , data=data_in) else
  pdat<-model.matrix(~scale(Disch.cfs) , data=data_in)
      
  #liost of data inputs for model
  data_list<-list(subyearlings=subyearlings,                                       #flag indicating if model is of subyearlings or yearling migrants
             rel=c(na.exclude(data_in[,paste0(lifestage,"_rel")])),   # releases in efficiency trils
             rec=c(na.exclude(data_in[,paste0(lifestage,"_recap")])), # recaps
             efficIndex=which(!is.na(
               data_in[,paste0(lifestage,"_rel")])),               # index of vector of trap days that efficiency trials correspond too
             Catch=data_in[,paste0("count.",lifestage)],           # vector of catch by trap day
             catch_DOY = data_in$DOY-min(data_in$DOY),             # vector of DOY of trap day
             first_DOY = min(data_in$DOY),                         # first day of year, for reference
             seriesFac = data_in$year_factor-1,                    # vector of year index of trap day
             years=data_in$Year,                                   # years, for reference
             pDat=pdat,                                            # design matrix for model of efficiency on each trap day
             N_trap=length(data_in[,paste0("count.",lifestage)]),  # number of trap days
             N_day=diff(range(data_in$DOY))+1,                     # number of days per year to estimate emigrant abundance
             Nyears=length(unique(data_in$year_factor)),           # number of years
             Use_NB=Use_NB)                                        # flag to use negative binomial observation model.
  
  return(data_list)
  }



# function to make data lists for all model of all streams and ages (subyearling and yearling) 
make_data_list_func<-function(chiw_data=chiw_data,nas_whi_data=nas_whi_data,Use_NB=1){
out<-list(6)
index<-1
for( i in c("Chiwawa","Nason","White")){
  for (j in c("sub","yrlng")){
    
    if(i=="Chiwawa") data_in<-chiw_data$dat else
      data_in<-nas_whi_data[[i]]
    
    if(j=="yrlng"){data_in<-subset(data_in,DOY<=200) #trim time series for yearlings to spring through early summer
       }
   
    
  out[[index]] <- make_screw_trap_model_data(data_in=data_in,stream=i, lifestage=j,Use_NB=Use_NB,subyearlings=ifelse(j=="sub",1,0))

  index <- index + 1
  }
}
return(out)
}



#function to make inital parameters for a model for a given data set
make_screw_trap_model_inits<-function(data_in){
  params<-list(pCoefs=rep(-.5,ncol(data_in$pDat)),
             mu_M=log(c(tapply(data_in$Catch *3,data_in$seriesFac,mean))),
             logit_phi_d=qlogis(.9),
             logit_phi_e=qlogis(.9),
             ln_tau_d=0,
             ln_tau_e=0,
             delta=numeric(data_in$N_day),
             epsilon=matrix(0,nrow=data_in$N_day,ncol=data_in$Nyears),
             logit_phi_NB=0)

return(params)
}



#function to fit model 
fit_model<-function(data_in){
setwd(here("src","TMB"))
TMB::compile("screw_trap_LP_3.cpp") #compile TMB model
dyn.load("screw_trap_LP_3") #load TMB model
params<-make_screw_trap_model_inits(data_in) #initial parameters
str(params)
if(data_in$Use_NB) map<-list() else map=list(logit_p_NB=factor(rep(NA,1))) #if using Poisson, dont optimize NB "prob" param 
mod<-TMB::MakeADFun(data_in,params,random=c("epsilon","delta"),DLL="screw_trap_LP_3",silent=T,map=map)# construct model
fit3<-TMBhelper::fit_tmb(mod,mod$fn,mod$gr, getsd = TRUE,newtonsteps = 1) #optimize model

#extract log sums of emigrants within life history migration windows, and standard errors
LH_sums<-LH_sums_sd<-NULL
try({
 LH_sums<-matrix(fit3$SD$value[names(fit3$SD$value)=="LH_sums"],ncol=ifelse(data_in$subyearlings,3,1))
 LH_sums_sd<-matrix(fit3$SD$sd[names(fit3$SD$value)=="LH_sums"],ncol=ifelse(data_in$subyearlings,3,1))
})
 
return(list(LH_sums=LH_sums,LH_sums_sd=LH_sums_sd,fit3=fit3))
}



#function to fit models for all streams and ages
fit_all<-function(all_data_lists){
chiw_subs<-fit_model(all_data_lists[[1]])
chiw_yrlngs<-fit_model(all_data_lists[[2]])
nason_subs<-fit_model(all_data_lists[[3]])
nason_yrlngs<-fit_model(data=all_data_lists[[4]])
white_subs<-fit_model(all_data_lists[[5]])
white_yrlngs<-fit_model(all_data_lists[[6]])

return(list(chiw_subs=chiw_subs, chiw_yrlngs=chiw_yrlngs,
            nason_subs=nason_subs, nason_yrlngs=nason_yrlngs,
            white_subs=white_subs, white_yrlngs=white_yrlngs
            ))
}



#function to plot geometric mean daily emigrant time series for a given natal stream
plot_average_ts_func<-function(sub_fit,yrlng_fit,data_ins,save_plot=FALSE,river){
log_means<-sub_fit$fit3$SD$value[names(sub_fit$fit3$SD$value)!="LH_sums"]
log_mean_sds<-sub_fit$fit3$SD$sd[names(sub_fit$fit3$SD$value)!="LH_sums"]

log_means_yrlngs<-yrlng_fit$fit3$SD$value[names(yrlng_fit$fit3$SD$value)!="LH_sums"]
log_mean_sds_yrlngs<-yrlng_fit$fit3$SD$sd[names(yrlng_fit$fit3$SD$value)!="LH_sums"]

if(save_plot) png(here("results","plots","chiw_daily_ts.png"),units="in",height=3,width=4,res=300)
par(cex=.8)
ymax<-max(c(exp(log_means+1.96*log_mean_sds),exp(log_means_yrlngs+1.96*log_mean_sds_yrlngs)))
{plot(1,type="n",xlim=c(0,500),ylim=c(0,ymax),xlab="",ylab="",xaxt="n",main=river)
polygon(c(1:length(log_means),length(log_means):1),c(exp(log_means+1.96*log_mean_sds),rev(exp(log_means-1.96*log_mean_sds))),border=F,col="grey")

polygon(c(1:length(log_means_yrlngs),length(log_means_yrlngs):1)+365,c(exp(log_means_yrlngs+1.96*log_mean_sds_yrlngs),rev(exp(log_means_yrlngs-1.96*log_mean_sds_yrlngs))),border=F,col="grey")

points(exp(log_means),type="l")
points(1:length(log_means_yrlngs)+365,exp(log_means_yrlngs),type="l")

abline(v=c(139,266)-data_ins$first_DOY,col="red")
abline(v=c(365)-data_ins$first_DOY,col="red")
days<-c(91,182,274,366,(365+91))-data_ins$first_DOY
labs<-c("Apr","Jul","Oct","Jan","Apr")
axis(1,at=days,labels=labs)
mtext("Emigrants/ day",2,3)
mtext("Average brood year",1,3)
text(c(91,182,300,520)-data_ins$first_DOY,y=ymax,pos=1,labels=c("Fry","Summer\nparr","Fall\nParr","Smolts"),cex=.9)}
if(save_plot) dev.off()
}



#function to plot geometric mean daily emigrant time series for all streams
plot_all_geomean_daily_emigrants<-function(all_emigrants_estimates){
  plot_average_ts_func(all_emigrants_estimates[["chiw_subs"]],all_emigrants_estimates[["chiw_yrlngs"]],all_data_lists[[1]],river="Chiwawa")
  plot_average_ts_func(all_emigrants_estimates[["nason_subs"]],all_emigrants_estimates[["nason_yrlngs"]],all_data_lists[[3]],river="Nason")
  plot_average_ts_func(all_emigrants_estimates[["white_subs"]],all_emigrants_estimates[["white_yrlngs"]],all_data_lists[[5]],river="White")
}


#end of script
