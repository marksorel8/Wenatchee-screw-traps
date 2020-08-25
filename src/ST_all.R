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
geo_means<-plot_all_geomean_daily_emigrants(all_emigrants_estimates)

#geometric mean subyearlings across streams
across_stream_geomean<-data.frame(x=1:365,y=exp(rowMeans(geo_means$all_log_means[,1:3],na.rm=T)))

#expand to daily observation for each (rounded) fish
observation_for_each_average_fish<-rep(across_stream_geomean$x[which(!is.na(across_stream_geomean$y))],times=round(na.exclude(across_stream_geomean$y)))

#fit three-normal mixture distribution
mix_geomean<-mixtools::normalmixEM(observation_for_each_fish,lambda=c(.2,.4,.5),mu=c(100,200,300),)

#calculate density
mix_geomean_dens<-mix_geomean$lambda[1]* dnorm(seq(50,350,by=1),mix_geomean$mu[1],mix_geomean$sigma[1])+
  mix_geomean$lambda[2]*dnorm(seq(50,350,by=1),mix_geomean$mu[2],mix_geomean$sigma[2])+
  mix_geomean$lambda[3]*dnorm(seq(50,350,by=1),mix_geomean$mu[3],mix_geomean$sigma[3])

#plot
plot(mix_geomean,2,breaks=200)
#hist(observation_for_each_fish,breaks=200,freq=FALSE)
points(seq(50,350,by=1),mix_geomean_dens,type="l",lwd=2)
breaks<-numeric(2)
for ( i in 1:2) breaks[i]<-which.min(mix_geomean_dens[((round(mix_geomean$mu[i]):round(mix_geomean$mu[(i+1)]))-49)])+round(mix_geomean$mu[i])
abline(v=breaks,lwd=2)

plot_all_geomean_daily_emigrants(all_emigrants_estimates,breaks=breaks)


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
  names(out)[[index]]<-paste(i,j,sep="_")
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
plot_average_ts_func<-function(sub_fit,yrlng_fit,data_ins,save_plot=FALSE,river,breaks=NULL){
log_means<-sub_fit$fit3$SD$value[names(sub_fit$fit3$SD$value)!="LH_sums"]
log_mean_sds<-sub_fit$fit3$SD$sd[names(sub_fit$fit3$SD$value)!="LH_sums"]

log_means_yrlngs<-yrlng_fit$fit3$SD$value[names(yrlng_fit$fit3$SD$value)!="LH_sums"]
log_mean_sds_yrlngs<-yrlng_fit$fit3$SD$sd[names(yrlng_fit$fit3$SD$value)!="LH_sums"]

if(save_plot) png(here("results","plots","chiw_daily_ts.png"),units="in",height=3,width=4,res=300)
par(cex=.8)
ymax<-max(c(exp(log_means+1.96*log_mean_sds),exp(log_means_yrlngs+1.96*log_mean_sds_yrlngs)))
{plot(1,type="n",xlim=c(0,500),ylim=c(0,ymax),xlab="",ylab="",xaxt="n",main=river)
  days<-seq(from=data_ins$first_DOY-50,by=1,length.out=length(log_means))
polygon(c(days,rev(days)),c(exp(log_means+1.96*log_mean_sds),rev(exp(log_means-1.96*log_mean_sds))),border=F,col="grey")

days_yrlng<-seq(from=data_ins$first_DOY-50,by=1,length.out=length(log_means_yrlngs))
polygon(c(days_yrlng,rev(days_yrlng))+365,c(exp(log_means_yrlngs+1.96*log_mean_sds_yrlngs),rev(exp(log_means_yrlngs-1.96*log_mean_sds_yrlngs))),border=F,col="grey")

points(days,exp(log_means),type="l")
points(days_yrlng+365,exp(log_means_yrlngs),type="l")

if(!is.null(breaks)){
abline(v=c(breaks[1],breaks[2])-50,col="red")
#abline(v=c(365)-data_ins$first_DOY,col="red")
text(c(91,182,300,520)-data_ins$first_DOY,y=ymax,pos=1,labels=c("Fry","Summer\nparr","Fall\nParr","Smolts"),cex=.9)
}
days<-c(91,182,274,366,(365+91))-50
labs<-c("Apr","Jul","Oct","Jan","Apr")
axis(1,at=days,labels=labs)
mtext("Emigrants/ day",2,3)
mtext("Average brood year",1,3)
}
if(save_plot) dev.off()
return(list(log_means=log_means,log_means_yrlngs=log_means_yrlngs,log_mean_sds=log_mean_sds,log_mean_sds_yrlngs,log_mean_sds_yrlngs,first_DOY=data_ins$first_DOY,N_day=data_ins$N_day))
}



#function to plot geometric mean daily emigrant time series for all streams
plot_all_geomean_daily_emigrants<-function(all_emigrants_estimates,breaks=NULL){
 chiwawa<- plot_average_ts_func(all_emigrants_estimates[["chiw_subs"]],all_emigrants_estimates[["chiw_yrlngs"]],all_data_lists[[1]],river="Chiwawa",breaks=breaks)
 nason<-plot_average_ts_func(all_emigrants_estimates[["nason_subs"]],all_emigrants_estimates[["nason_yrlngs"]],all_data_lists[[3]],river="Nason",breaks=breaks)
 white<- plot_average_ts_func(all_emigrants_estimates[["white_subs"]],all_emigrants_estimates[["white_yrlngs"]],all_data_lists[[5]],river="White",breaks=breaks)

#create matrix with average daily emigrants for each stream and age
all_log_means<-matrix(NA,nrow=365,ncol=6)
for ( i in 1:3){
  data_in<-get(c("chiwawa","nason","white")[i]) 
all_log_means[seq(from=data_in$first_DOY,by=1,length.out=length(data_in$log_means)),i]<-data_in$log_means
all_log_means[seq(from=data_in$first_DOY,by=1,length.out=length(data_in$log_means_yrlngs)),(3+i)]<-data_in$log_means_yrlngs
}

return(list(chiwawa=chiwawa,nason=nason,white=white,all_log_means=all_log_means))
 }


#function to plot sums and standard errors


#end of script
