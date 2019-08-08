pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}#end of function

pkgTest("here")

par(mfrow=c(1,1))

load_dat<-function(){
  
  source(here("src","Chiwawa Data Proc.R"))
  
  chiw_dat<-Chiw_dat_Proc()
  
  source(here("src","Nason White Data Proc.R"))
  
  Nas_Whi_Dat<-Nason_White_data_Func()
  
  return(list(chiw=chiw_dat,nas=list(nas_effic=Nas_Whi_Dat$effic_dat$Nason_ET,nas_catch=Nas_Whi_Dat$catch_ops_dat$Nason_catch,nasDis=Nas_Whi_Dat$effic_dat$disch_Dat$Nason_Dis),whi=list(nhi_effic=Nas_Whi_Dat$effic_dat$White_ET,whi_catch=Nas_Whi_Dat$catch_ops_dat$White_catch,whiDis=Nas_Whi_Dat$effic_dat$disch_Dat$White_Dis)))
}
