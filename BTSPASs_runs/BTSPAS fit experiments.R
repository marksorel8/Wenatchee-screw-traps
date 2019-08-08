pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}#end of function

pkgTest("here")

#Load data
source(here("src","Load Screw Trap Data.R"))
screw_trap_dat<-load_dat()

#load function for making time series
source(here("src","ts_and_plotting_funcs.R"))



#make and plot chiwawa subyearling time series and discharge time series
Chiw_sub_ts<-make_ts(screw_trap_dat$chiw$chiw_catch$EndDate,screw_trap_dat$chiw$chiw_catch$subCatch,"Chiwawa River PARR",plot=TRUE,plotDis = TRUE,dis_dates =screw_trap_dat$chiw$chiwDis$Date ,dis_vals=screw_trap_dat$chiw$chiwDis$X_00060_00003 )

#make and plot chiwawa yearling time series and discharge time series
Chiw_yrlng_ts<-make_ts(screw_trap_dat$chiw$chiw_catch$EndDate,screw_trap_dat$chiw$chiw_catch$yrlngCatch,"Chiwawa River YEARLINGS",plot=TRUE,plotDis = TRUE,dis_dates =screw_trap_dat$chiw$chiwDis$Date ,dis_vals=screw_trap_dat$chiw$chiwDis$X_00060_00003 )




#make inputs for PTSPAS
time<-125:(365+175)

n1<-c(window(make_ts(screw_trap_dat$chiw$chiw_effic$Date2, screw_trap_dat$chiw$chiw_effic$rel)$catch ,start=c(2015,125),end=c(2015,365)),window(make_ts(screw_trap_dat$chiw$chiw_effic$Date2, screw_trap_dat$chiw$chiw_effic$rel)$catch ,start=c(2016,1),end=c(2016,175)))

n1[is.na(n1)]<-0

m2<-c(window(make_ts(screw_trap_dat$chiw$chiw_effic$Date2, screw_trap_dat$chiw$chiw_effic$recap)$catch ,start=c(2015,125),end=c(2015,365)),window(make_ts(screw_trap_dat$chiw$chiw_effic$Date2, screw_trap_dat$chiw$chiw_effic$recap)$catch ,start=c(2016,1),end=c(2016,175)))

m2[is.na(m2)]<-0

u2<-c(window(Chiw_sub_ts$catch,start=c(2015,125),end=c(2015,365)),window(Chiw_yrlng_ts$catch,start=c(2016,1),end=c(2016,175)))

u2[is.na(u2)]<-0

Chiw_2015_BTSPAS<-TimeStratPetersenDiagError_fit("Chiw_15-16",
                                             prefix="test",
                                             time=time,
                                             n1=n1,
                                             m2=m2,
                                             u2=u2,
                                             InitialSeed = 6231988,
                                             logitP.cov = rep(1, length(125:(365+175))),
                                             bad.n1 =which(is.na(n1))+124 ,
                                             bad.u2 = which(is.na(u2))+124,
                                             sampfrac = rep(1,length(time)))
