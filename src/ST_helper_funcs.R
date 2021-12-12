#function to make data lists for model
make_screw_trap_model_data<-function(data_in,stream="Chiwawa", lifestage="yrlng",Use_NB=1,subyearlings=0){
  
  #subset data to days when there was catch data (including 0 catch)
  data_in<-filter(data_in,!is.na(data_in[,paste0("count.",lifestage)]))
  
  #design matrix for trap efficiency model
  
  pdat<-glmmTMB::glmmTMB(DOY~scale(Disch.cfs)+diag(1+scale(Disch.cfs)|year_factor) +(1|Week)+(1|Week:year_factor),data=data_in,dispformula = ~0,doFit=FALSE)
  #
  #list of data inputs for model
  data_list<-list(subyearlings=subyearlings,                                       #flag indicating if model is of subyearlings or yearling migrants
                  rel=c(na.exclude(data_in[,paste0(lifestage,"_rel")])),   # releases in efficiency trils
                  rec=c(na.exclude(data_in[,paste0(lifestage,"_recap")])), # recaps
                  efficIndex=which(!is.na(
                    data_in[,paste0(lifestage,"_rel")])),               # index of vector of trap days that efficiency trials correspond too
                  Catch=data_in[,paste0("count.",lifestage)],           # vector of catch by trap day
                  catch_DOY = data_in$DOY-min(data_in$DOY),             # vector of DOY of trap day
                  first_DOY = min(data_in$DOY),                         # first day of year, for reference
                  seriesFac = data_in$year_factor-1,                    # vector of year index of trap day
                  years=data_in$Year,                             # years, for reference
                  beta= pdat$parameters$beta,
                  b= pdat$parameters$b,
                  theta = pdat$parameters$theta,
                  terms_p = pdat$data.tmb$terms,
                  X = pdat$data.tmb$X,
                  Z = pdat$data.tmb$Z,
                  
                  # design matrix for model of efficiency on each trap day
                  N_trap=length(data_in[,paste0("count.",lifestage)]),  # number of trap days
                  N_day=diff(range(data_in$DOY))+1,                     # number of days per year to estimate emigrant abundance
                  Nyears=length(unique(data_in$year_factor)),           # number of years
                  Use_NB=Use_NB,                             # flag to use negative binomial observation model.
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



#function to make initial parameters for a model for a given data set
make_screw_trap_model_inits<-function(data_in){
  params<-list(#pCoefs=rep(-.5,ncol(data_in$pDat)),
               beta=data_in$beta,
               b=data_in$b,
               theta=data_in$theta,
               mu_M=log(c(tapply(data_in$Catch *3,data_in$seriesFac,mean))),
               logit_phi_d=qlogis(.9),
               logit_phi_e=qlogis(.9),
               ln_tau_d=0,
               ln_tau_e=0,
               delta=numeric((data_in$N_day)-1),
               epsilon=matrix(0,nrow=data_in$N_day,ncol=data_in$Nyears),
               log_phi_NB=0,
               ln_tau_edp=0,
               logit_phi_edp=qlogis(.9),
               eps_doy_p=rep(.1,data_in$N_day))
  
  return(params)
}



#function to fit model 
fit_model<-function(data_in,get_jp=TRUE ,do_boot=FALSE){
  setwd(here("src","TMB"))
  TMB::compile("screw_trap_LP_4.cpp") #compile TMB model
  dyn.load("screw_trap_LP_4") #load TMB model
  params<-make_screw_trap_model_inits(data_in) #initial parameters
  str(params)
  if(data_in$Use_NB) map<-list() else map=list(log_phi_NB=factor(rep(NA,1))) #if using Poisson, dont optimize NB "prob" param
  # map$logit_phi_d<-factor(NA)
  # map$logit_phi_e<-factor(NA)
  mod<-TMB::MakeADFun(data_in,params,random=c("epsilon","delta","b"),DLL="screw_trap_LP_4",silent=T,map=map)# construct model
  fit3<-TMBhelper::fit_tmb(mod,mod$fn,mod$gr, getsd = TRUE,newtonsteps = 1,getJointPrecision = get_jp) #optimize model
  
  #extract log sums of emigrants within life history migration windows, and standard errors
  LH_sums<-LH_sums_sd<-NULL
  try({
    LH_sums<-matrix(fit3$SD$value[names(fit3$SD$value)=="LH_sums"],ncol=ifelse(data_in$subyearlings,3,1))
    LH_sums_sd<-matrix(fit3$SD$sd[names(fit3$SD$value)=="LH_sums"],ncol=ifelse(data_in$subyearlings,3,1))
    rownames(LH_sums)<-rownames(LH_sums_sd)<-levels(as.factor(data_in$years))
  })
  
    
    
  if(do_boot){
  #bootstrap log juvenile sums and sds
  boot<-bootstrap_juves(data_in, mod$env$last.par.best, fit3$SD$jointPrecision, n_sim=10000,seed=1234)
  }else{boot=NULL}
  
  return(list(LH_sums=LH_sums,LH_sums_sd=LH_sums_sd,boot=boot,fit3=fit3,mod=mod,M_hat=mod$report()$M_hat, dat=data_in))
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

#function to boostrap juvenile abundances
bootstrap_juves<-function(dat, last_best , joint_precis ,n_sim, seed){
  # browser()
  test_sim<-rmvnorm_prec(last_best ,
                         joint_precis ,
                         n_sim,seed)
  
  
  Nyears<-dat$Nyears
  Ndays<-dat$N_day
  first_rand<-min(which(names(last_best)=="delta"))
  first_year_rand<-first_rand+(Ndays-1)
  first_day<-dat$first_DOY
  breaks<-dat$breaks
  stream_length<-dat$stream_length
  
  sim_array<-array(NA,dim=c(Ndays,n_sim,Nyears))
  for ( year in 1:Nyears){
    sim_array[,,year]<- (matrix(rep(test_sim[year,1:n_sim],each=Ndays),ncol=n_sim)+
                              rbind(matrix(0,nrow=1,ncol=n_sim),test_sim[first_rand:(first_year_rand-1),1:n_sim])+ 
                              test_sim[seq((first_year_rand+((year-1)*Ndays)),length=Ndays),1:n_sim])
  }
  
  LH_test<-apply(sim_array,2:3,function(x){
    c(sum(exp(x[1:(breaks[1]-first_day)])),
      sum(exp(x[((breaks[1]-first_day)+1):(breaks[2]-first_day)])),
      sum(exp(x[(breaks[2]-first_day+1):Ndays])))
  })/stream_length
  
  quants_geom_mean_DAY<-apply(exp(apply(sim_array,1:2,mean)),1,quantile,probs=c(.025,.5,.975))
  
  return(list(
    boot_log_means=t(apply(LH_test,c(1,3),function(x)mean(log(x)))),
    boot_log_sds=t(apply(LH_test,c(1,3),function(x)sd(log(x)))),
    quants_geom_mean_DAY=quants_geom_mean_DAY#,
    # big_array=LH_test
  ))
  
}


#function to fit models for all streams and ages
fit_all<-function(all_data_lists,do_boot){
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



#function to calculate geometric mean daily emigrants accross all streams and years
across_stream_year_geommean_func<-function(all_em_ests=all_emigrants_estimates,all_data_lists,do_plot=FALSE,breaks=NULL,labels=FALSE){
  
  out<-matrix(NA,nrow=600,ncol=0) # zero column matrix to cbind to
  for ( i in c(1,3,5)){ #loop through streams
    out_2<-matrix(NA,nrow=600,ncol=all_data_lists[[i]]$Nyears) # matrix (365 * nYears[stream])
    out_2[seq(all_data_lists[[i]]$first_DOY,by=1,length.out=all_data_lists[[i]]$N_day),]<-
      all_em_ests[[i]]$M_hat # fill matrix with days where subyearling emigrants estimated
    out_2[seq(all_data_lists[[(i+1)]]$first_DOY,by=1,length.out=all_data_lists[[(i+1)]]$N_day)+365,]<-
      all_em_ests[[(i+1)]]$M_hat # fill matrix with days where yearling emigrants estimated 
    
    out<-cbind(out,out_2)                # cbind to other streams
  }
  
  if(do_plot){
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
  }
  return(out)
}


