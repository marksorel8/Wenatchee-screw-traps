



#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#Function to  plot latent spawners vs juveniles with functional relationship and process error envelope
# returns a ggplot
#----------------------------------------------------------------------------------------

ggplot_spawner_juveniles<-function(mod_fit , mod_dat, mod_rep){

  #tibble of estimated spawner and emigrant abundances
  sum_out<- tibble(
    juveniles=mod_fit$SD$value[names(mod_fit$SD$value)=="log(J_pred)"],
    juveniles_sd=mod_fit$SD$sd[names(mod_fit$SD$value)=="log(J_pred)"],
    stream=  mod_dat$s_i,
    LH= mod_dat$l_i,
    t=mod_dat$t_i,
    S_fac=mod_dat$st_i) %>% left_join(  
      tibble(
        spawners=mod_fit$SD$value[names(mod_fit$SD$value)=="log_S_hat"],
        spawner_sd=mod_fit$SD$sd[names(mod_fit$SD$value)=="log_S_hat"],
        S_fac=mod_dat$st_i %>% levels())) %>% 
    mutate(stream=c("Chiwawa","Nason","White")[stream+1],
           stream=fct_relevel(stream,c("Chiwawa","Nason","White")),
           LH=c("Spr.0","Sum.0","Fal.0","Spr.1")[LH+1],
           LH=fct_relevel(LH,c("Spr.0","Sum.0","Fal.0","Spr.1")))
  
  
  #tibble of model prediction of emigrant abundances over a range of spawner abundances
  preds<-tibble(
    alpha=mod_rep$alpha,
    gamma=mod_rep$gamma,
    Jmax=mod_rep$Jmax,
    stream=rep(0:2,each=4),
    LH=rep(0:3,times=3),
    loadings=abs(c(log(mod_rep$Loadings_pf[1]),mod_rep$Loadings_pf[-1])),
    idio_var=mod_rep$sigma_eta^2) %>% 
    mutate(eps=sqrt(loadings+idio_var)) %>% 
    mutate(stream=c("Chiwawa","Nason","White")[stream+1],
           stream=fct_relevel(stream,c("Chiwawa","Nason","White")),
           LH=c("Spr.0","Sum.0","Fal.0","Spr.1")[LH+1],
           LH=fct_relevel(LH,c("Spr.0","Sum.0","Fal.0","Spr.1")))%>% 
    left_join(
      sum_out %>% group_by(stream) %>% summarize (S_max=max(exp(spawners+2*spawner_sd)))) %>% 
    crossing(
      spawners=seq(from=0,to=max(exp(sum_out$spawners+2*sum_out$spawner_sd)),by=.1)) %>% 
    filter(spawners<=S_max) %>% 
    mutate(juveniles=((alpha*(spawners)^(gamma)/
                         (1+alpha*(spawners)^(gamma)/Jmax))))
  
  
  
  
  #Begin plot of spawners vs juveniles 
  scale_juv<-100
  #facet_wrap plot
  SR_plot<-ggplot(data=preds,aes(x= spawners,y=juveniles/scale_juv))+facet_wrap(~LH+stream,scales = "free",nrow=4) + geom_ribbon(aes(ymin=juveniles*exp(-1.96*eps)/scale_juv,
                                                                                                                                     ymax=juveniles*exp(1.96*eps)/scale_juv,fill=rgb(.7,.1,.1,.2)),show.legend = FALSE)+
    geom_line(color="firebrick3",size=1.25)+
    geom_point(data=sum_out,aes(x=exp(spawners),y=exp(juveniles )/scale_juv))+ 
    geom_linerange(data=sum_out,aes(x=exp(spawners),y=exp(juveniles )/scale_juv,ymin=exp(qnorm(.025,juveniles ,juveniles_sd ))/scale_juv,ymax=exp(qnorm(.975,juveniles ,juveniles_sd ))/scale_juv))+
    geom_linerange(data=sum_out,aes(x=exp(spawners),y=exp(juveniles )/scale_juv,xmin=exp(qnorm(.025,spawners,spawner_sd)),xmax=exp(qnorm(.975,spawners,spawner_sd)))) + 
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      plot.margin = margin(25.5, 50, 5.5, 5.5, "pt"),
      panel.spacing.y=unit(-.25,"cm"),
      legend.box.spacing=unit(1.5,"cm"))+
    xlab("Spawners")+ylab("Emigrants x100")+
    scale_fill_manual(values=c(rgb(.7,.1,.1,.2),rgb(.7,.1,.1,.3)),labels=c("Idiosyncratic"," Correlated"),name="Prediction interval")+guides(fill = guide_legend(override.aes= list(alpha = c(0.2,.3))))
  
  
  #grid_wrap plot to get facet labels grobs from
  p_grid<-ggplot(data=preds,aes(x= spawners,y=juveniles))+facet_grid(LH~stream,scales = "free")+geom_ribbon(aes(ymin=juveniles*exp(-1.96*eps),
                                                                                                                ymax=juveniles*exp(1.96*eps)),fill=rgb(.7,.1,.1,.2)) + 
    theme(
      plot.margin = margin(2, 2, 2, 2, "cm")
    ) 
  
  
  # function to selectively remove certain grobs
  gtable_filter_remove <- function (x, name, trim = TRUE){
    matches <- !(x$layout$name %in% name)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    if (trim) 
      x <- gtable_trim(x)
    x
  }
  
  #convert facet_wrap plot to Grob list
  p_grid_tab<-ggplotGrob(p_grid)
  #convert facet_grid plot to Grob list
  SR_plot_tab<-ggplotGrob(SR_plot)
  
  # remove bottom axes from all but smolt facets
  SR_plot_filtered<-gtable_filter_remove(SR_plot_tab,
                                         name = c(paste0("axis-b-",rep(1:3,each=3) ,"-", 1:3),
                                                  paste0("strip-t-",rep(1:3,each=3) ,"-", 2:4)),
                                         trim = FALSE)
  
  #add facet labels for columns (Stream)
  SR_plot_filtered<-gtable::gtable_add_grob(SR_plot_filtered, p_grid_tab$grobs[grep('strip-t', p_grid_tab$layout$name)],t=rep(1,3), l= c(5,9,13))
  
  #add facet labels for rows (LHs)
  SR_plot_filtered<-gtable::gtable_add_grob(SR_plot_filtered, p_grid_tab$grobs[grep('strip-r', p_grid_tab$layout$name)],t=c(8,13,18,23),l=rep(17,4) )
  
  #adjust placement of facet labels
  SR_plot_filtered$heights[7]<-unit(-2,"pt")
  SR_plot_filtered$widths[15]<-unit(-.45,"cm")
  
  #render plots
  grid::grid.newpage()
  grid::grid.draw(SR_plot_filtered)

return(list(sum_out=sum_out,preds=preds))
}


#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Function to plot process error correlations
## returna a plot
#----------------------------------------------------------------------------------------

plot_cor<-function(mod_rep,mod_fit, mod_dat,nsim=100000){

lh_s_names<-paste(rep(c("Spr-0","Sum-0","Fall-0","Spr-1"),times=3),rep(c("Chiw.","Nason","White"),each=4),sep=" ")
loadings2<-cbind(mod_rep$Loadings_pf,diag(mod_rep$sigma_eta))
# loadings2<-mod_rep$Loadings_pf
#marg_var<-rowSums(loadings2^2)
cor_mat<-cov2cor( loadings2%*%t(loadings2)) 
dimnames(cor_mat)<-list(lh_s_names,lh_s_names)


bootstrap_corr<-function(mle=mod_fit$par,cov_mat=mod_fit$SD$cov.fixed,n_f=mod_dat$n_f,n_sl=mod_dat$n_sl){
  
  #empty array to hold bootstrapped correlation matrices  
  corr_array<-array(NA,dim=c(n_sl,n_sl,nsim))  
  
  #empty matrix to hold fact loadings and idiosyncratic error loadings
  loadings_mat<-matrix(0,n_sl,n_f+n_sl)
  #indices of lower triangular components of vector of factor loadings matrix to fill
  loadings_mat_ind<-which(lower.tri(loadings_mat[,1:n_f],diag=TRUE))
  #indices of factor loading in vector of parameters
  load_vec_ind<-which(names(mle)=="Loadings_vec")
  #indices idosyncratic process errors loadings in vector of paramaters
  log_sig_eta<-which(names(mle)=="log_sigma_eta")
  
  #indices of very small process errors with very large standard errors to fix at zero
  # fix_zero_ind<-which(mle[log_proc_sig_ind]<(-10))
  
  # draw parameter values from multivariate normal defined by mle and inverted hessian (covariance matrix)
  sim<-mvtnorm::rmvnorm(n=nsim,mean=mle,sigma=cov_mat,checkSymmetry=FALSE)
  
  #calculate and store error correlation for each parameter vector drawn above
  for ( i in 1:nsim){
    loadings_mat[loadings_mat_ind]<-sim[i,load_vec_ind] # fill in factor loadings
    
    diag(loadings_mat[,(n_f+1):(n_f+n_sl)])<-exp(sim[i,log_sig_eta]) # fill in idosyncratic loadings
    diag(loadings_mat[,(n_f+1):(n_f+n_sl)])[diag(loadings_mat[,(n_f+1):(n_f+n_sl)])==Inf]<-10000
    #fix certain idisyncratic loadings at zero
    # loadings_mat[(fix_zero_ind),(fix_zero_ind+n_f)]<-0
    
    corr_array[,,i]<-cov2cor( loadings_mat%*%t(loadings_mat)) 
  }

  corr_array
}

# corr_boot<-bootstrap_corr()

# p_corr<-apply(corr_boot,1:2,function(x){sum(sign(x)<=0)/nsim})*sign(cor_mat)+as.numeric(sign(cor_mat)<=0)

corrplot::corrplot(cor_mat,type="lower",method = "circle",
                   # p.mat =p_corr*2,  insig = "p-value", sig.level = -.05,
                   
                   mar=c(0,0,0,0),oma=c(0,0,0,0),tl.cex=1.5,tl.col="black",cl.cex=1.5,tl.srt=45,diag=FALSE)

}

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#Functon to plot spawners and juveniles by year
#----------------------------------------------------------------------------------------
plot_spawn_em_year<-function(mod_dat,sum_out){
# spawners/Km
spawner_dat<-tibble(Spawners= exp(sum_out$spawners) ,stream=c("Chiwawa","Nason","White")[(mod_dat$s_i+1)],BY=mod_dat$BY) %>% left_join(pHOS %>% rename(BY=1,stream=Stream)) %>% distinct()

sum_spawn<-spawner_dat %>% group_by(BY) %>% summarise(Spawners=sum(Spawners)/n(),n=n(),pHOS_weighted=sum(pHOS_weighted)/n()) %>% ungroup() %>% mutate(stream="Average") %>% filter(n==3)

spawner_dat<-bind_rows(spawner_dat,sum_spawn)%>% 
  #drop birst and last brood years which have incomplete juvenile data
  group_by(stream) %>% mutate(max_BY=max(BY),min_BY=min(BY)) %>% ungroup() %>% filter(BY<max_BY,BY>min_BY) %>% mutate(Wild=Spawners*(1-pHOS_weighted),Hatchery=Spawners*pHOS_weighted) %>% pivot_longer(c(Wild,Hatchery),names_to = "LH") %>% select(-Spawners) %>% rename(Spawners=value)


#juveniles/km
juvenile_dat<-bind_cols(Juveniles  = exp(sum_out$juveniles)/100,stream=c("Chiwawa","Nason","White")[mod_dat$s_i+1],BY=mod_dat$BY,LH=c("Spr-0","Sum-0","Fall-0","Spr-1")[mod_dat$l_i+1]) %>% mutate(LH=fct_relevel(LH,"Spr-0","Sum-0","Fall-0","Spr-1")) %>%   #drop first and last brood years which have incomplete juvenile data
  group_by(stream) %>% mutate(max_BY=max(BY),min_BY=min(BY)) %>% ungroup() %>% filter(BY<max_BY,BY>min_BY)

sum_juv<-juvenile_dat %>% group_by(BY,LH) %>% summarise(Juveniles = sum(Juveniles)/n(),count=n()) %>% ungroup() %>% mutate(stream="Average") %>% filter(count==3)

juvenile_dat<-bind_rows(juvenile_dat,sum_juv) %>% group_by(BY,stream) %>% mutate(prop=Juveniles/sum(Juveniles))

all_dat<-bind_rows(juvenile_dat,spawner_dat) %>%  pivot_longer(c(Juveniles,prop,Spawners))%>% mutate(name=fct_relevel(name,"Spawners","Juveniles","prop")) %>%  mutate(stream=fct_relevel(stream,"Chiwawa","Nason","White"))%>%  mutate(stream=fct_relevel(stream,"Chiwawa","Nason","White")) %>% mutate(LH=fct_relevel(LH,"Spawner","Spr-0","Sum-0","Fall-0","Spr-1"))

#stacked counts
juv_2<-ggplot(data=all_dat,
              aes(fill=LH, x=BY, y=value))+ geom_bar(stat="identity", width=.8)+facet_grid(name~ stream, scales="free_y", switch="y", labeller = labeller(name=function(x){c("Spawners/ km)", "Emigrants x100/ km","Proportion")}))+  scale_fill_discrete(name="Life Stage", type="viridis")+ theme_grey()+  ylab(NULL) + labs(x="Brood Year")+
  theme(strip.background = element_blank(),
        strip.placement = "outside")+ scale_x_continuous(breaks =seq(2000,2020,10))+ scale_fill_manual(values=c(viridis(option="B",end=.9,begin=.225,n=4),"grey40","black"))

juv_2
}


#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#Functon to plot expected emigrants vs. spawners
#----------------------------------------------------------------------------------------

expec_spawn_em_ploft_func<-function(preds){
##expected emigrants vs. spawners

pred_mat<-preds %>%  mutate(juveniles =juveniles /100,LH=fct_relevel(LH,"Spr-0","Sum-0","Fall-0","Spr-1")) %>% group_by(stream,spawners) %>%  mutate(prop=juveniles/sum(juveniles)) %>% pivot_longer(c(juveniles,prop)) %>% mutate(name=as_factor(name),name= fct_relevel(name,c("prop","juveniles")))
  

juvenile_plot<-ggplot(data=pred_mat,
                      aes(x=spawners/10,y=value,fill=LH))+geom_bar(stat="identity",width=.02)+facet_grid(name~stream, scales="free",switch="y", space="free_x", labeller = labeller(name=function(x){c( "Proportion","Emigrants x100/ km")}))+ scale_fill_discrete( name="Life Stage")+ theme_grey()+  ylab(NULL) + labs(x="Spawners (x10/ km)")+
  theme(strip.background = element_blank(),
        strip.placement = "outside")+ scale_fill_viridis(option="B",discrete=TRUE,end=.9,begin=.225)

juvenile_plot
}


#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#function to bootstrap and plot environmental covariates coefficient values
## returns a ggplot
#----------------------------------------------------------------------------------------

bootstrap_env_cov<-function(dat, last_best=fit$par , precis=fit$SD$cov.fixed ,n_sim=50000){
  
  test_sim<-mvtnorm::rmvnorm(n_sim,last_best ,
                             precis )
  
  row_beta_e<-which(names(last_best)=="beta_e")[1:8]
  
  out<-cbind(LH=c("Spr-0","Sum-0","Fall-0","Spr-1","Sum-0","Fall-0","Spr-1"),
             season=c(rep("Winter 1 max discharge",4),rep("Summer 1 mean discharge",3),"Winter 2\n max disch."),t(test_sim[,row_beta_e])) %>% as_tibble() %>%  pivot_longer(!c(LH,season),names_to=NULL,values_to="value") %>% mutate(value=as.numeric(value)) %>% mutate(LH=fct_relevel(LH,"Spr-0","Sum-0","Fall-0","Spr-1"),season=fct_relevel(season,"Winter 1 max discharge","Summer 1 mean discharge" ,"Winter 2\n max disch."))
  
  out<-ggplot(data=out,aes(x=LH,y=value)) + facet_grid(~season,scale="free_x", space = "free_x")+ geom_hline(yintercept=0,linetype=2)+geom_violin(fill="black")+ xlab("Life History") +ylab("Coefficient value")
  
  return(out)
}

