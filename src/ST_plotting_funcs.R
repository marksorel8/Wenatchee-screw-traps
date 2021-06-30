
#function to create a ggplot of average daily emigrants with mixture distribution and breaks for LHP emigration windows
ggplot_timing_func<-function(all_emigrants_estimates,all_data_lists,across_stream_geomean,mix_geomean_dens, breaks){
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
    label = c("Spr-0","Sum-0","Fall-0","Spr-1"),
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

#------------------------------------------------------------------------------

#function to plot average discharg and temperature
plot_dis_temp_func<-function(all_data_lists){
  
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
  
  dish_temp<-full_join(all_dis,temp_dat) %>% pivot_longer(c(flow,temp_c)) %>% mutate(name=case_when(name=="flow"~"Discharge~(m^3*s^{-1})",TRUE~"Temperature~( degree*C )")) %>% rename(Stream=stream)
  
  
  temp_disch_plot<-ggplot(data=dish_temp,aes(x=doy,y=value))+facet_wrap(~name,scales = "free_y",nrow=2,strip.position = "left",labeller = label_parsed) + geom_line(size=1.,aes(color=Stream))+ ylab("")+xlab("") + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335,366),labels=c("Jan","","","Apr","","","Jul","","","Oct","","","Jan"))+ theme(strip.background = element_blank(),strip.placement = "outside")+ scale_color_viridis(option="D",discrete=T,end=.7)
  
  
  temp_disch_plot
}
