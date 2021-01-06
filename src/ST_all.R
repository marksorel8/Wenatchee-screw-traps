# Mark Sorel

# Function estimates the number of juvenile emigrants from natal streams in the Wenatchee River basin 

#takes data collected from screw traps provided by Washington Department of Fish and Wildlife (Chiwawa River) and Yakama Nation Fisheries (Nason and White Rivers) and returns estimates of daily emigrants and log sums of emigrants within migrations windows representing juvenile life histories (Fry, summer parr, fall parr, yearling smolt), and standard errors.

ST_all_func<-function(refit=FALSE){

#------------------------------------------------------------------------------------
# Load function for analysis
#------------------------------------------------------------------------------------

#functions to read and process data
source(here("src","Chiwawa Data Proc 2.R")) # Chiwawa 
source(here("src","Nason White Data Proc 2.R")) # Nason and White

# functions to fit screw trap model for all streams and ages, calculate geometric means of daily emigrants, fit a mixture distribution and come up with breakpoints between LHPs, and bootstrap total emigrants within each window.
source(here("src","ST_helper_funcs.R")) # Chiwawa 

#------------------------------------------------------------------------------------
# Calling functions to conduct analysis.
#------------------------------------------------------------------------------------

#process data
chiw_data<-Chiw_dat_Proc() #chiwawa
nas_whi_data<-Nason_White_data_Func() #nason and white

#make data lists for all models
all_data_lists<-make_data_list_func(chiw_data = chiw_data, nas_whi_data = nas_whi_data)


#load most recent evaluation of daily juvenile migrant abundance model if exists
if(!refit & length(list.files(here("results"))[substr(list.files(here("results")),1,18)=="emigrant_estimates"])>0){

  last_fit_file<-list.files(here("results"))[substr(list.files(here("results")),1,18)=="emigrant_estimates"][which.max(lubridate::mdy(substr(list.files(here("results"))[substr(list.files(here("results")),1,18)=="emigrant_estimates"],20,30)))]
  
  print(paste("Loading file:", last_fit_file ))
  
load(file=here("results",last_fit_file))
  
}else{
  
#fit all models (takes ~3 minutes)
ts<-Sys.time() 
all_emigrants_estimates<-fit_all(all_data_lists=all_data_lists,do_boot=FALSE)
Sys.time()-ts

#save as ".Rdata" object
save(all_emigrants_estimates,file=here("results",paste0("emigrant_estimates",substr(date(),4,10),substr(date(),20,25),".Rdata")))
}

#calculate geomeans of emigrants for each day of year across years
geo_means<-across_stream_year_geommean_func(all_em_ests=all_emigrants_estimates,all_data_lists=all_data_lists)
  
#transform the average number of emigrants on each day of year to a data set with an obervation for each average fish (rounded), for fitting a mixture distribution.
##extract geometric mean number of emigrants per day of year
across_stream_geomean<-data.frame(x=1:600,y=exp(rowMeans(log(geo_means),na.rm=T)))
##expand to daily observation for each (rounded) fish
observation_for_each_average_fish<-rep(across_stream_geomean$x[which(!is.na(across_stream_geomean$y))],times=round(na.exclude(across_stream_geomean$y)))

#fit three-normal mixture distribution
LL<- -Inf
for ( i in 1:5){ #fit three times to ensure convergence
mix_geomean_i<-mixtools::normalmixEM(observation_for_each_average_fish,lambda=c(.2,.4,.5,.5),mu=c(100,200,300,375),)
if(mix_geomean_i$loglik>LL){
  mix_geomean<-mix_geomean_i
  LL<-mix_geomean$loglik
}
}

#calculate density
mix_geomean_dens<-mix_geomean$lambda[1]* dnorm(seq(50,565,by=1),mix_geomean$mu[1],mix_geomean$sigma[1])+
  mix_geomean$lambda[2]*dnorm(seq(50,565,by=1),mix_geomean$mu[2],mix_geomean$sigma[2])+
  mix_geomean$lambda[3]*dnorm(seq(50,565,by=1),mix_geomean$mu[3],mix_geomean$sigma[3])+
  mix_geomean$lambda[4]*dnorm(seq(50,565,by=1),mix_geomean$mu[4],mix_geomean$sigma[4])

# find saddle points (local minima) between means of mixture componeent distributions
breaks<-numeric(3)
for ( i in 1:3) breaks[i]<-which.min(mix_geomean_dens[((round(sort(mix_geomean$mu)[i]):round(sort(mix_geomean$mu)[(i+1)]))-49)])+round(sort(mix_geomean$mu)[i])



#calculate bootstrapped distributions of total emigrants expressing each LHP in each stream (this takes ~ 10 minutes) if not already done
if(is.null(all_emigrants_estimates[[1]]$boot)){
  ts<-Sys.time() 
  all_emigrants_estimates<-lapply(all_emigrants_estimates,function(x){
  x$dat$breaks<-breaks[1:2]
  x$boot<-bootstrap_juves(x$dat,x$mod$env$last.par.best,x$fit$SD$jointPrecision,n_sim=50000,seed=1234)
  return(x)
})
  Sys.time()-ts
  #save as ".Rdata" objects
  save(all_emigrants_estimates,file=here("results",paste0("emigrant_estimates",substr(date(),4,10),substr(date(),20,25),".Rdata")))
}

return(list(all_emigrants_estimates=all_emigrants_estimates,
            all_data_lists=all_data_lists,
            across_stream_geomean=across_stream_geomean,
            mix_geomean_dens=mix_geomean_dens,
            breaks=breaks))

}

