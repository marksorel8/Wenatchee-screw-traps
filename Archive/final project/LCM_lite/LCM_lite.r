library(here)
load(here("emigrants.Rdata"))
Chiw_redds<-read.csv(here("ChiwEscapement.csv"))
Chiw_spawners<-rowSums(Chiw_redds[7:29,2:3])


dim(all_boots_trans)

em_ins<-array(NA,dim=c(4,22,2000))


load(here("lognorm_pars.Rdata"))

for ( i in 1:4){
  em_ins[i,,]<-all_boots_trans[,,i]
}
dim(em_ins)
the_data <- list(wild_spawners=Chiw_redds[7:29,2],
                 hatchery_spawners=Chiw_redds[7:29,3],
                 em_obs_log_means=lognorm_pars[[1]],
                 em_obs_log_sds=lognorm_pars[[2]])
                 
                 #emigrants=em_ins)


load(here("par_outs.Rdata"))

init_cov_mat<-matrix(0.1,nrow=4,ncol=4)

diag(init_cov_mat)<-unlist(lapply(par_outs,function(x)x$ proc_sigma))


inits<-list(log_S_hat=log(the_data$wild_spawners+the_data$hatchery_spawner),
                 S_hat_CV=.05,
                 log_em_hat=log(matrix(unlist(lapply(par_outs,function(x)x$R_hat)),ncol=4)),
                 em_hat_CV=matrix(unlist(lapply(par_outs,function(x)x$R_obs_sd)),ncol=4),
                 em_proc_cov=init_cov_mat,
                 alpha=unlist(lapply(par_outs,function(x)x$alpha)),
                 R_max=unlist(lapply(par_outs,function(x)x$R_max)),
                 d=unlist(lapply(par_outs,function(x)x$d))[-4]
)


library(rstan),init=list(inits)
fit <- stan(file = 'LCM_lite.stan', data = the_data, 
            iter = 10000, chains = 1,verbose=F,control=list(adapt_delta=.9))


print(fit,pars=c("alpha","R_max","d"))

inits$log_S_hat

la <- extract(fit,permuted = TRUE)


apply(la$R_max,2,mean)

em<-apply(la$em_pred,2:3,mean)
sp<-exp(apply(la$log_S_hat,2,mean))

plot(sp[-1],em[,4])

par(mfrow=c(2,2))
hist(la$em_proc_cov[,3,3])
hist(la$mean_beta)