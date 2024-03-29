pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}#end of function

Chiw_dat_Proc<-function(){
  #load  packages
  pkgTest("here")
  pkgTest( "tidyverse" )
  
  #download discharge data from internet
  source(here("src","Discharge data funcs.R"))
  if(file.exists(here("chiwDis.csv"))){
    load(here("chiwDis.csv"))}else{
    chiwDis<-Chiw_discharge_func()
    save(chiwDis,file=here("chiwDis.csv"))  
  }
  
 ################################################################## 
# process efficiency trial data
  
  ChiwEfficTrials<-read.csv(here("data","Chiwawa","Chiwawa_efficiency_trials_2.csv"))   %>% 
    select(c(-1,-7)) %>% #drop model year and efficiency

  # change column names
  rename("Lifestage"="Species" ,"Date"="X1st.day.of.collections","rel"="X..Released","recap"="X..Recaptured","Disch.cfs"="Discharge..cfs.") %>% 
    

#drop some rows that have comments suggesting that the data is not valid
  filter( !Notes %in% c("Cant find any trial similar to flow, sample size in 2007?","Cant find anything close to these groups","crossed out on datasheet","Nothing close in records","YCW from Lake Trap trial","Only used YCW. From Lake Trap trial, not sure why it was used in Chiwawa Model - JW")) %>% 
  
    
   #format recapture dates 
mutate(Date= as.Date(as.character(Date),format = ifelse(grepl("-",as.character(Date)),"%d-%b-%y","%m/%d/%Y"))) %>% 
  
    #add day, week, and year
    mutate(DOY=as.numeric(format(Date,form="%j"))) %>% 
    mutate(Week=as.numeric(format(Date,form="%W"))) %>% 
    mutate(Year=as.numeric(format(Date,form="%Y"))) %>% 
  
  
#rename lifestages YCW=yearlings, SBC=subyearling chinook
    mutate(LifeStage=recode(Lifestage,"YCW only"="YCW","Only used YCW"="YCW","SBC only"="SBC","Most likely SBC based on date"="SBC"),
           LifeStage=case_when(Notes%in%c("Only used YCW","YCW only")~"YCW",
                                      Notes%in%c("Most likely SBC based on date","SBC only")~"SBC",
                               DOY<150~"YCW",
                                      TRUE~LifeStage
                  )) %>% 
    filter(LifeStage%in%c("YCW","SBC")) %>% 
    
    arrange(LifeStage) %>% #sort so that "YCW" and "SBC" comes before "YCW & SBC" so when duplicates removed, "YCW & SBC" are removed first
  
  #drop duplicate rows
 distinct(Date,.keep_all = TRUE) %>% 
    
    #make columns with releases and recaps for trials that included subyearling or yearlings, with trials that had both included in both columns. 
    mutate(sub_rel=ifelse(LifeStage%in%c("SBC"),rel,NA),
           sub_recap=ifelse(LifeStage%in%c("SBC"),recap,NA),
           yrlng_rel=ifelse(LifeStage%in%c("YCW"),rel,NA),
           yrlng_recap=ifelse(LifeStage%in%c("YCW"),recap,NA)) %>% 
    
  
  #make a new "position" column where "low flow" is changed to "upper", because only a few data points for "low flow"
  mutate(Position2 = as.character(recode(Position,"Low Flow"="Upper")))
  
  

###########################################################
###########################################################
# process dialy catch and operations data
  
  #read data on trap operations
  trapOps<-read.csv(here("data","Chiwawa","Chiw.trap.ops.csv")) %>% 
  #format "endDate" to date 
  mutate(EndDate=as.Date(EndDate,format="%m/%d/%Y")) %>% 
    rename(Date=EndDate) %>% 
  
  #get rid of days when trap was stopped or partially stopped or not running at all.
  filter(Status!="S" &Status!="P" & Position!="Screw Stopper" & Position!="Partial Trapping"& Position!="Out"& Position!="IN") %>% 
  
  #make capitalization consistant
  mutate(Position=recode_factor(Position,"LOWER"="Lower","UPPER"="Upper")) %>% 
  #make empty position cells = "unknown"
    mutate(Position=recode_factor(Position,"Lower"="Lower","Upper"="Upper","Low Flow"="Low Flow",.default = "Unknown")) %>% 
    
     #make new "posotion" column where "low flow" is  changed to "upper", because only a few data points for "low flow".
  mutate(Position2 = as.character(recode(Position,"Low Flow"="Upper"))) %>% 
    
  #add day, week, and year
  mutate(DOY=as.numeric(format(Date,form="%j"))) %>% 
  mutate(Week=as.numeric(format(Date,form="%W"))) %>% 
  mutate(Year=as.numeric(format(Date,form="%Y")))  %>% 
                  
  #drop levels
  droplevels()
  

  #*********   Guess at missing trap positions based on discharge   **********
  #*********   should be estimated in model in future   ***
  
  #  tapply(ChiwEfficTrials$disch.cfs,ChiwEfficTrials$Position2,summary)
  # 
  # trapRunDays$Position2[trapRunDays$Position2==""& trapRunDays$dis<=475]<-"Upper"
  # trapRunDays$Position2[trapRunDays$Position2==""]<-"Lower"
  # trapRunDays$Position2<-factor(trapRunDays$Position2,levels=c("Lower","Upper"))
  # 
  # trapRunDays<-droplevels(trapRunDays)
  
  
  
  ################################################  
  # load catch data
  
  #function to get rid of leap day
  drop_leap_day<-function(x){
    if(length(which(x[8,-1]>0))>0){
    x[7,which(x[8,-1]>0)+1]<-sum( x[7:8,which(x[8,-1]>0)+1])}
    x<-x[-8,]
  }
  
  #function to reformat chiwawa catch data
  munge_chiw_catch<-function(x){
    #drop leap day
    drop_leap_day(x) %>% 
      ##change from wide to long format
      #fry
      gather(year,count,-1) %>% 
      #adjust column-name and year values
      rename("day"=X,Year=year) %>% 
      mutate(Year=as.numeric(substr(Year,2,5))) %>% 
      mutate(Date= as.Date(paste(day,Year,sep="-"),format = "%d-%b-%Y") ) %>% 
      select(-1)
    
  }
  
  
  ##load daily catch data by "lifestage"
  #fry
  Chiw_Cnt<-read.csv(here("data","Chiwawa","ChiwFryCnt.csv")) %>% 
    munge_chiw_catch() %>% 
  
  #load subyearlings and join with fry
    full_join((read.csv(here("data","Chiwawa","ChiwSubCatch.csv")) %>% munge_chiw_catch() ),by=c("Year", "Date"),suffix=c(".fry",".parr")) %>% 
      
  #add a column of fry+subyearlings
      mutate(count.sub=count.fry+count.parr) %>% 

    #load yearlings and join with subyearlings  
    full_join((read.csv(here("data","Chiwawa","ChiwYrlngCnt.csv")) %>% munge_chiw_catch() ),by=c("Year", "Date")) %>% 
    rename(count.yrlng=count) %>% 
    
  #Merge catch data with operations data
    full_join(trapOps,by=c("Date","Year")) %>% 
    
    mutate(DOY=as.numeric(format(Date,form="%j"))) %>% 
    mutate(Week=as.numeric(format(Date,form="%W"))) %>%
    
    #fill in the few missing discharges with data obtained with dataRetrieval package
    mutate(Disch.cfs=coalesce(Mean.discharge..CFS., as.integer(round(chiwDis$flow[match(Date,chiwDis$date-1)])))) %>% 
    
    
  dplyr::full_join(select(ChiwEfficTrials,sub_rel:yrlng_recap,DOY,Week,Year)) %>%  #join with efficiency trial data
    
  mutate(year_factor=as.numeric(as.factor(Year))) %>% #year as factor
    select(sub_rel:yrlng_recap,count.sub,count.yrlng,DOY,Week,Year,year_factor,Position2,Disch.cfs) #select certain columns that will be used in modeling
  
  

  return(list(dat=Chiw_Cnt))

  }#end of function


