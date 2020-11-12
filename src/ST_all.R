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
load(file=here("results","Rdata",paste0("emigrant_estimates Sep 08 2020",".Rdata")))
ts<-Sys.time() 
all_emigrants_estimates_4<-fit_all(all_data_lists=all_data_lists)
Sys.time()-ts
lapply(all_emigrants_estimates_4,function(x)(x$fit3$AIC))

#save as ".Rdata" objects
save(all_emigrants_estimates,file=here("results",paste0("emigrant_estimates",substr(date(),4,10),substr(date(),20,25),".Rdata")))

#plot geomeans of emigrants for each day of year across years
png(file=here("results","plots","average-timing.png"),units="in",width = 6,height=6,res=300)
par(mfrow=c(2,2),mar=c(3,3,2,1),oma=c(2.5,2.5,1,0))
geo_means<-plot_all_geomean_daily_emigrants(all_emigrants_estimates)
mtext("Emigrants/ day",2,1,outer=T)
#mtext("Average brood year",1,1,outer=T)
dev.off()
#geometric mean subyearlings across streams
across_stream_geomean<-data.frame(x=1:365,y=exp(rowMeans(geo_means$all_log_means[,1:3],na.rm=T)))
across_stream_geomean<-data.frame(x=1:600,y=exp(rowMeans(log(geo_means$geomean_all[1:600,]),na.rm=T)))
plot(across_stream_geomean,type="l")
#expand to daily observation for each (rounded) fish
observation_for_each_average_fish<-rep(across_stream_geomean$x[which(!is.na(across_stream_geomean$y))],times=round(na.exclude(across_stream_geomean$y)))

#fit three-normal mixture distribution
mix_geomean<-mixtools::normalmixEM(observation_for_each_average_fish,lambda=c(.2,.4,.5,.5),mu=c(100,200,300,375),)

#calculate density
mix_geomean_dens<-mix_geomean$lambda[1]* dnorm(seq(50,565,by=1),mix_geomean$mu[1],mix_geomean$sigma[1])+
  mix_geomean$lambda[2]*dnorm(seq(50,565,by=1),mix_geomean$mu[2],mix_geomean$sigma[2])+
  mix_geomean$lambda[3]*dnorm(seq(50,565,by=1),mix_geomean$mu[3],mix_geomean$sigma[3])+
  mix_geomean$lambda[4]*dnorm(seq(50,565,by=1),mix_geomean$mu[4],mix_geomean$sigma[4])
#plot
png(file=here("Methods","Splits.png"),res=300,units="in",height=4,width=5)
#plot(mix_geomean,2,breaks=100,main="")
hist(observation_for_each_average_fish,breaks=100,freq=T,main="",ylab="Density (Average emigrants)",xlab="Day of year")
points(seq(50,565,by=1),mix_geomean_dens*sum(across_stream_geomean$y,na.rm=T),type="l",lwd=2,col="firebrick3")
breaks<-numeric(3)
for ( i in 1:3) breaks[i]<-which.min(mix_geomean_dens[((round(sort(mix_geomean$mu)[i]):round(sort(mix_geomean$mu)[(i+1)]))-49)])+round(sort(mix_geomean$mu)[i])
abline(v=breaks,lwd=1,col="firebrick3",lty=3)
dev.off()


png(file=here("results","plots","average-timing_breaks.png"),units="in",width = 6,height=6,res=300)
par(mfrow=c(2,2),mar=c(3,3,2,1),oma=c(2.5,2.5,1,0))
x<-plot_all_geomean_daily_emigrants(all_emigrants_estimates,breaks=breaks,lab_LH = TRUE)
points(seq(50,550,by=1),mix_geomean_dens*sum(across_stream_geomean$y,na.rm=T),type="l",lwd=2,col="firebrick3")
mtext("Emigrants/ day",2,1,outer=T,adj=.7)

dev.off()



png(file=here("results","plots","temp_dis.png"),units="in",width = 4.5,height=3,res=300)
plot_dis_temp_func()
dev.off()

png(file=here("results","plots","temp_dis.png"),units="in",width = 6,height=4,res=300)
ggplot_timing_func()
dev.off()
#------------------------------------------------------------------------------------
# Function for analysis
#------------------------------------------------------------------------------------

#functions to read and process data
source(here("src","Chiwawa Data Proc 2.R")) # Chiwawa 
source(here("src","Nason White Data Proc 2.R")) # Nason and White



#function to make data lists for model
make_screw_trap_model_data<-function(data_in,stream="Chiwawa", lifestage="yrlng",Use_NB=1,subyearlings=0,breaks=c(139,262)){
  
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
             Use_NB=Use_NB,
             breaks=breaks,                                        # flag to use negative binomial observation model.
             stream_length = ifelse(stream=="Chiwawa",32.5,
                                    ifelse(stream=="Nason",15.4,
                                           16.1)))                 # tributary lengths (km) for standardizing juvenile emigrants  
  
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
             log_phi_NB=0)

return(params)
}



#function to fit model 
fit_model<-function(data_in,get_jp=FALSE ){
setwd(here("src","TMB"))
TMB::compile("screw_trap_LP_3.cpp") #compile TMB model
dyn.load("screw_trap_LP_3") #load TMB model
params<-make_screw_trap_model_inits(data_in) #initial parameters
str(params)
if(data_in$Use_NB) map<-list() else map=list(logit_p_NB=factor(rep(NA,1))) #if using Poisson, dont optimize NB "prob" param 
mod<-TMB::MakeADFun(data_in,params,random=c("epsilon","delta"),DLL="screw_trap_LP_3",silent=T,map=map)# construct model
fit3<-TMBhelper::fit_tmb(mod,mod$fn,mod$gr, getsd = TRUE,newtonsteps = 1,getJointPrecision = get_jp) #optimize model

#extract log sums of emigrants within life history migration windows, and standard errors
LH_sums<-LH_sums_sd<-NULL
try({
 LH_sums<-matrix(fit3$SD$value[names(fit3$SD$value)=="LH_sums"],ncol=ifelse(data_in$subyearlings,3,1))
 LH_sums_sd<-matrix(fit3$SD$sd[names(fit3$SD$value)=="LH_sums"],ncol=ifelse(data_in$subyearlings,3,1))
 rownames(LH_sums)<-rownames(LH_sums_sd)<-levels(as.factor(data_in$years))
})

#bootstrap log juvenile sums and sds
boot<-bootstrap_juves(data_in, mod$env$last.par.best, fit3$SD$jointPrecision, n_sim=10000,seed=1234)

 
return(list(LH_sums=LH_sums,LH_sums_sd=LH_sums_sd,boot=boot,fit3=fit3,mod=mod,M_hat=mod$report()$M_hat))
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

#function to boostrap juvenile abundances. IN PROGRESS
bootstrap_juves<-function(dat, last_best , joint_precis ,n_sim, seed){
test_sim<-rmvnorm_prec(last_best ,
                        joint_precis ,
                        n_sim,seed)
Nyears<-dat$Nyears
Ndays<-dat$N_day
first_rand<-min(which(names(last_best)=="delta"))
first_year_rand<-first_rand+Ndays
first_day<-dat$first_DOY
breaks<-dat$breaks
stream_length<-dat$stream_length

sim_array<-array(NA,dim=c(Ndays,n_sim,Nyears))
for ( year in 1:Nyears){
  sim_array[,,year]<- exp(matrix(rep(test_sim[year,1:n_sim],each=Ndays),ncol=n_sim)+test_sim[first_rand:(first_year_rand-1),1:n_sim]+test_sim[seq((first_year_rand+((year-1)*Ndays)),length=Ndays),1:n_sim])
}

LH_test<-apply(sim_array,2:3,function(x){
  c(sum(x[1:(breaks[1]-first_day)]),
    sum(x[((breaks[1]-first_day)+1):(breaks[2]-first_day)]),
    sum(x[(breaks[2]-first_day+1):Ndays]))
})/stream_length

return(list(
  boot_log_means=t(apply(LH_test,c(1,3),function(x)mean(log(x)))),
  boot_log_sds=t(apply(LH_test,c(1,3),function(x)sd(log(x)))),
  big_array=LH_test
))

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
plot_average_ts_func<-function(sub_fit,yrlng_fit,data_ins,save_plot=FALSE,river,breaks=NULL,labels=FALSE,lab_LH=FALSE){
log_means<-sub_fit$fit3$SD$value[names(sub_fit$fit3$SD$value)!="LH_sums"]
log_mean_sds<-sub_fit$fit3$SD$sd[names(sub_fit$fit3$SD$value)!="LH_sums"]

log_means_yrlngs<-yrlng_fit$fit3$SD$value[names(yrlng_fit$fit3$SD$value)!="LH_sums"]
log_mean_sds_yrlngs<-yrlng_fit$fit3$SD$sd[names(yrlng_fit$fit3$SD$value)!="LH_sums"]

if(save_plot) png(here("results","plots","chiw_daily_ts.png"),units="in",height=3,width=4,res=300)
par(cex=.8)
ymax<-max(c(exp(log_means+1.96*log_mean_sds),exp(log_means_yrlngs+1.96*log_mean_sds_yrlngs)))
plot(1,type="n",xlim=c(0,500),ylim=c(0,ymax),xlab="",ylab="",xaxt="n",main=river)
  days<-seq(from=data_ins$first_DOY-50,by=1,length.out=length(log_means))
polygon(c(days,rev(days)),c(exp(log_means+1.96*log_mean_sds),rev(exp(log_means-1.96*log_mean_sds))),border=F,col="grey")

days_yrlng<-seq(from=data_ins$first_DOY-50,by=1,length.out=length(log_means_yrlngs))
polygon(c(days_yrlng,rev(days_yrlng))+365,c(exp(log_means_yrlngs+1.96*log_mean_sds_yrlngs),rev(exp(log_means_yrlngs-1.96*log_mean_sds_yrlngs))),border=F,col="grey")

points(days,exp(log_means),type="l")
points(days_yrlng+365,exp(log_means_yrlngs),type="l")

if(!is.null(breaks)){
abline(v=c(breaks[1],breaks[2],breaks[3])-50,col="firebrick3",lty=3)
}
#abline(v=c(365)-data_ins$first_DOY,col="red")
if(lab_LH) text(c(91,205,320,520)-data_ins$first_DOY,y=ymax,pos=1,labels=c("Fry","Summ.","Fall","Smolts"),cex=.9)

days<-c(91,182,274,366,(365+91))-50
labs<-c("Apr","Jul","Oct","Jan","Apr")
axis(1,at=days,labels=labs)
if(labels){
  mtext("Emigrants/ day",2,3)
mtext("Average brood year",1,3)
}
if(save_plot) dev.off()
return(list(log_means=log_means,log_means_yrlngs=log_means_yrlngs,log_mean_sds=log_mean_sds,log_mean_sds_yrlngs,log_mean_sds_yrlngs,first_DOY=data_ins$first_DOY,N_day=data_ins$N_day))
}



#function to calculate geometric mean daily emigrants accross all streams and years
across_stream_year_geommean_func<-function(all_em_ests=all_emigrants_estimates,all_dat=all_data_lists,breaks=NULL,labels=FALSE){

  out<-matrix(NA,nrow=600,ncol=0) # zero column matrix to cbind to
  for ( i in c(1,3,5)){ #loop through streams
    out_2<-matrix(NA,nrow=600,ncol=all_dat[[i]]$Nyears) # matrix (365 * nYears[stream])
    out_2[seq(all_dat[[i]]$first_DOY,by=1,length.out=all_dat[[i]]$N_day),]<-
      all_em_ests[[i]]$M_hat # fill matrix with days where subyearling emigrants estimated
    out_2[seq(all_dat[[(i+1)]]$first_DOY,by=1,length.out=all_dat[[(i+1)]]$N_day)+365,]<-
      all_em_ests[[(i+1)]]$M_hat # fill matrix with days where yearling emigrants estimated 
    
    out<-cbind(out,out_2)                # cbind to other streams
  }
  
  plot(exp(rowMeans(log(out),na.rm=T)),type="l",lwd=2,xlim=c(50,550),xlab="",ylab="",xaxt="n",main="Average")
  days<-c(91,182,274,366,(365+91))
  labs<-c("Apr","Jul","Oct","Jan","Apr")
  axis(1,at=days,labels=labs)
  if(!is.null(breaks)){
    abline(v=c(breaks[1],breaks[2],breaks[3]),col="firebrick3",lty=3)
  }
  
  if(labels){
  mtext("Emigrants/ day",2,3)
  mtext("Average brood year",1,3)}
  return(out)
}




#ggplot of average daily emigrants
## THis is an ugly functionw that relys on a bung of stuff from the global environement. TODO fix
ggplot_timing_func<-function(){
#long tibble of geometric mean number of emigrants on each day of year
daily_ave_mat<-tibble(juveniles=unlist(lapply(all_emigrants_estimates, function(x)x$fit3$SD$value[names(x$fit3$SD$value)=="mean_day_log_M"])),
                      juvenile_se=unlist(lapply(all_emigrants_estimates, function(x)x$fit3$SD$sd[names(x$fit3$SD$value)=="mean_day_log_M"])),
                      doy=unlist(lapply(all_data_lists,function(x) seq(from = x$first_DOY,length.out =  x$N_day)-365*(x$subyearlings-1))),
                      age=unlist(lapply(all_data_lists,function(x) rep(x$subyearlings, x$N_day))),
                      stream=apply(matrix(1:6),2,function(x)rep(c("Chiwawa","Nason","White")[ceiling(x/2)],unlist(lapply(all_data_lists,function(x) x$N_day))[x]))) %>% rename(stream=5)

#geometric mean accross all days and years from the global environment
#add the fit of the mixture distributions (also from the gloabl environement)
ave<-across_stream_geomean %>% as_tibble() %>% rename(doy=x,juveniles=y) %>% filter(doy>=53 & doy<=565) %>% mutate(stream="Average",juveniles=log(juveniles),age=ifelse(doy<365,0,1),mix_dist=(mix_geomean_dens*sum(across_stream_geomean$y,na.rm=T))[doy-49]) 

#combine average and stream specific
daily_ave_mat<-bind_rows(daily_ave_mat,ave) %>% mutate(stream=fct_relevel(stream,"Chiwawa","Nason","White","Average"))

#text of lif ehistories to go on first fecet pannel
dat_text <- data.frame(
  label = c("Fry", "Summ.", "Fall","Smolt"),
  stream   = factor(rep("Chiwawa",4),levels = c("Chiwawa","Nason","White","Average")),
  x     = c(91,200,325,535),
  y     = rep(635,4)
)

#plot
ggplot(data=daily_ave_mat,aes(x=doy,y=exp(juveniles)))+geom_ribbon(aes(ymin=exp(juveniles-1.96*juvenile_se   ),ymax=exp(juveniles+1.96*juvenile_se),group=age),fill=rgb(.5,.5,.5,.8))+
  geom_path(aes(group=age),size=1.01)+
  facet_wrap(~stream,scale="free_y")+scale_x_continuous(breaks=c(1,91,182,274,366,456,547),labels=c("Jan","Apr","Jul","Oct","Jan","Apr","Jul"))+xlab("")+ylab("Emigrants/ day")+geom_path(aes(x=doy,y=mix_dist),color=rgb(.8,.1,.1,.7),size=1.1)+geom_vline(xintercept=breaks,linetype=2,color=rgb(.8,.1,.1,.7))+ geom_text(
    data    = dat_text,
    mapping = aes(x = x, y = y, label = label)
  )

}


#function to plot geometric mean daily emigrant time series for all streams
plot_all_geomean_daily_emigrants<-function(all_emigrants_estimates,breaks=NULL,lab_LH=FALSE){
 chiwawa<- plot_average_ts_func(all_emigrants_estimates[["chiw_subs"]],all_emigrants_estimates[["chiw_yrlngs"]],all_data_lists[[1]],river="Chiwawa",breaks=breaks,lab_LH=lab_LH)
 nason<-plot_average_ts_func(all_emigrants_estimates[["nason_subs"]],all_emigrants_estimates[["nason_yrlngs"]],all_data_lists[[3]],river="Nason",breaks=breaks)
 white<- plot_average_ts_func(all_emigrants_estimates[["white_subs"]],all_emigrants_estimates[["white_yrlngs"]],all_data_lists[[5]],river="White",breaks=breaks)

#create matrix with average daily emigrants for each stream and age
all_log_means<-matrix(NA,nrow=365,ncol=6)
for ( i in 1:3){
  data_in<-get(c("chiwawa","nason","white")[i]) 
all_log_means[seq(from=data_in$first_DOY,by=1,length.out=length(data_in$log_means)),i]<-data_in$log_means
all_log_means[seq(from=data_in$first_DOY,by=1,length.out=length(data_in$log_means_yrlngs)),(3+i)]<-data_in$log_means_yrlngs
}

#geomean of all 
geomean<-across_stream_year_geommean_func(breaks=breaks)

return(list(chiwawa=chiwawa,nason=nason,white=white,all_log_means=all_log_means,geomean_all=geomean))
 }


#function to plot average discharge
plot_dis_temp_func<-function(breaks){
  
  source(here("src","Discharge data funcs.R"))
  if(file.exists(here("chiwDis.csv"))){
    load(here("chiwDis.csv"))}else{
      chiwDis<-Chiw_discharge_func()
      save(chiwDis,file=here("chiwDis.csv"))  
    }

years<- range(all_data_lists[[1]]$years)
Chiw_disch<-chiwDis %>% 
  filter(Year>=years[1]&Year<=years[2]) %>% 
  group_by(doy) %>% summarise(mean(flow,na.rm=T)) %>% 
  mutate(stream="Chiwawa")




if(file.exists(here("nasWhiteDis.csv"))){
  load(here("nasWhiteDis.csv"))}else{
    disch_Dat<-Nason_White_Discharge_Func()
    save(disch_Dat,file=here("nasWhiteDis.csv"))  
  }

years_Nas<-range(all_data_lists[[3]]$years)
years_Whi<-range(all_data_lists[[5]]$years)
nas_dis<-disch_Dat$Nason_Dis %>%
  mutate(doy=as.numeric(format(date,form="%j"))) %>% 
  filter(Year>=years_Nas[1]&Year<=years_Nas[2]) %>% 
  group_by(doy) %>% summarise(mean(flow,na.rm=T)) %>% 
  mutate(stream="Nason")

whi_dis<-disch_Dat$White_Dis %>%
  mutate(doy=as.numeric(format(date,form="%j"))) %>% 
  filter(Year>=years_Whi[1]&Year<=years_Whi[2]) %>% 
  group_by(doy) %>% summarise(mean(flow,na.rm=T)) %>% 
  mutate(stream="White")


all_dis<-bind_rows(Chiw_disch,nas_dis,whi_dis) %>% 
  mutate(`mean(flow, na.rm = T)`=`mean(flow, na.rm = T)`*0.0283168) %>% rename(flow=`mean(flow, na.rm = T)`)



temp_dat<-read_csv(here("data","Temp","Wen_data_day.csv")) %>% subset(StreamName %in% c("Chiwawa River","Nason Creek","White River")) %>% group_by(StreamName) %>% filter(Elevation==min(Elevation)) %>%   group_by(StreamName,JulianDate) %>% summarize(mean(AvgDailyTemp,na.rm=T)) %>% `colnames<-`(c("stream","doy","temp_c")) %>% mutate(stream=gsub(" .*$","",stream))

#range in years of temp data for each stream
read_csv(here("data","Temp","Wen_data_day.csv")) %>% subset(StreamName %in% c("Chiwawa River","Nason Creek","White River")) %>% group_by(StreamName) %>% filter(Elevation==min(Elevation)) %>% summarise(minY=min(Year),maxY=max(Year))

dish_temp<-full_join(all_dis,temp_dat) %>% pivot_longer(c(flow,temp_c)) %>% mutate(name=case_when(name=="flow"~"Discharge~(m^3*s^2)",TRUE~"Temperature~( degree*C )")) %>% rename(Stream=stream)


temp_disch_plot<-ggplot(data=dish_temp,aes(x=doy,y=value))+facet_wrap(~name,scales = "free_y",nrow=2,strip.position = "left",labeller = label_parsed) + geom_line(size=1.,aes(color=Stream))+ ylab("")+xlab("") + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335,366),labels=c("Jan","","","Apr","","","Jul","","","Oct","","","Jan"))+ theme(strip.background = element_blank(),strip.placement = "outside")+ scale_color_viridis(option="D",discrete=T,end=.7)


temp_disch_plot
}


#end of script

