
Chiw_dat_Proc<-function(){
  
  
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}#end of function

#start function to process efficiency trial data
Chiw_effic_dat_proc<-function(drop_unk_lifestage=FALSE){
  
  #read data
  pkgTest("here")
  library(here)
  ChiwEfficTrials<-read.csv(here("data","Chiwawa","Chiw_Efficiency_Trials.csv"))
  
  # change column names
  colnames(ChiwEfficTrials)<-c("modYear","lifeStage","Date","position","rel","recap","effic","disch.cfs")
  dim(ChiwEfficTrials)
  
  #reorider the rows so that the ones with actual dates come first, that way when we remove duplicates we dont remove the ones with full dates
  ChiwEfficTrials<-ChiwEfficTrials[order(nchar(as.character(ChiwEfficTrials$Date)),decreasing = T),]
  
  #drop duplicate rows
  ChiwEffic<-ChiwEfficTrials[!duplicated(paste(ChiwEfficTrials$rel,ChiwEfficTrials$recap,ChiwEfficTrials$disch.cfs,ChiwEfficTrials$position)),]
  
  
  
  if(drop_unk_lifestage==TRUE){
  #drop efficiency trials for both lifestages combined
  
  ChiwEffic<-droplevels(subset(ChiwEffic,lifeStage!="YCW & SBC"))
  }

  #make a new "position" column where "low flow" is changed to "upper", because only a few data points for "low flow"
  ChiwEffic$position2<-ChiwEffic$position

  ChiwEffic$position2[ChiwEffic$position=="Low Flow"]<-"Upper"

  ChiwEffic<-droplevels(ChiwEffic)
  
  #meand and standard deviation of discharge
  disMean<-attributes(scale(ChiwEffic$disch.cfs))$`scaled:center`
  
  disScale<-attributes(scale(ChiwEffic$disch.cfs))$`scaled:scale`
  
  return(list(ChiwEffic=ChiwEffic,scaleDis=c(disMean,disScale)))
  
}#end function



#start function to process dialy catch and operations data
Chiw_catch_dat_Proc<-function(disch_Scale){
  pkgTest("here")
  library(here)
  #read data on trap operations
  trapOps<-read.csv(here("data","Chiwawa","Chiw.trap.ops.csv"))
  #format "endDate" to date 
  trapOps$EndDate<-as.Date(trapOps$EndDate,format="%m/%d/%Y")
  
  #get rid of days when trap was stopped or partially stopped or not running at all.
  trapRunDays<-droplevels(subset(trapOps,Status!="S" &Status!="P" & Position!="Screw Stopper" & Position!="Partial Trapping"& Position!="Out"& Position!="IN"))
  
  #make capitalization consistant
  trapRunDays$Position[trapRunDays$Position=="LOWER"]<-"Lower"
  trapRunDays$Position[trapRunDays$Position=="UPPER"]<-"Upper"
  
  #make new "posotion" column where "low flow" is  changed to "upper", because only a few data points for "low flow".
  trapRunDays$Position2<-trapRunDays$Position
  trapRunDays$Position2[trapRunDays$Position=="Low Flow"]<-"Upper"
  
  #change factor to character
  trapRunDays$Position<-as.character(trapRunDays$Position)
  trapRunDays$Position2<-as.character(trapRunDays$Position2)
  #change "NA" in position columns to "unknown"
  trapRunDays$Position[which(is.na(trapRunDays$Position))]<-"Unknown"
  trapRunDays$Position2[which(is.na(trapRunDays$Position2))]<-"Unknown"
  
  trapRunDays<-droplevels(trapRunDays)
  
  
  #fill in the few missing discharges with data obtained with dataRetrieval package
  
  #download discharge data from internet
  source(here("src","Discharge data funcs.R"))
  chiwDis<-Chiw_discharge_func()
  
  #fill in missing discharges with downloaded data
  trapRunDays$dis<-trapRunDays$Mean.discharge..CFS.
  trapRunDays$dis[is.na(trapRunDays$Mean.discharge..CFS.)]<-
    chiwDis$X_00060_00003[match(trapRunDays$EndDate[is.na(trapRunDays$Mean.discharge..CFS.)],chiwDis$date)]
  
################################################  
  ##load daily catch data by "lifestage
  #fry
  ChiwFryCnt<-read.csv(here("data","Chiwawa","ChiwFryCnt.csv"))
  
  #get rid of leap day
  drop_leap_day<-function(x){
    if(length(which(x[8,-1]>0))>0){
    x[7,which(x[8,-1]>0)+1]<-sum( x[7:8,which(x[8,-1]>0)+1])}
    x<-x[-8,]
  }
  ChiwFryCnt<-drop_leap_day(ChiwFryCnt)
  
  #subyearlings
  ChiwSubCnt<-read.csv(here("data","Chiwawa","ChiwSubCatch.csv"))
  
  ChiwSubCnt<-drop_leap_day(ChiwSubCnt)
  
  #yearling
  ChiwYrlngCnt<-read.csv(here("data","Chiwawa","ChiwYrlngCnt.csv"))
  
  ChiwYrlngCnt<-drop_leap_day( ChiwYrlngCnt)
  
  ##change from wide to long format
  pkgTest("reshape")
  library("reshape")
  #fry
  ChiwFryCntLn<-melt(ChiwFryCnt,id="X")
  #adjust column-name and year values
  colnames(ChiwFryCntLn)<-c("day","year","count")
  ChiwFryCntLn$year<-as.numeric(substr(ChiwFryCntLn[,2],2,5))
  
  #subyearling
  ChiwSubCntLn<-melt(ChiwSubCnt,id="X")
  colnames(ChiwSubCntLn)<-c("day","year","count")
  ChiwSubCntLn$year<-as.numeric(substr(ChiwSubCntLn[,2],2,5))
  #yearlings
  ChiwYrlngCntLn<-melt(ChiwYrlngCnt,id="X")
  colnames(ChiwYrlngCntLn)<-c("day","year","count")
  ChiwYrlngCntLn$year<-as.numeric(substr(ChiwYrlngCntLn[,2],2,5))
 
  #######################################################
  #Merge catch- with operations data
  
  ChiwFryCntLn$date<-as.Date(paste(ChiwFryCntLn$day,ChiwFryCntLn$year,sep="-"),format = "%d-%b-%Y")
  
  ChiwSubCntLn$date<-as.Date(paste(ChiwSubCntLn$day,ChiwSubCntLn$year,sep="-"),format = "%d-%b-%Y")
  
  ChiwYrlngCntLn$date<-as.Date(paste(ChiwYrlngCntLn$day,ChiwYrlngCntLn$year,sep="-"),format = "%d-%b-%Y")
  
  
  
  #add fry to trop ops
  trapRunDays$fryCatch<-ChiwFryCntLn$count[match(trapRunDays$EndDate,ChiwFryCntLn$date)]
  
  #add subs
  trapRunDays$subCatch<-ChiwSubCntLn$count[match(trapRunDays$EndDate,ChiwSubCntLn$date)]

  #subs+fry
  trapRunDays$allSubs<-trapRunDays$fryCatch+trapRunDays$subCatch

  #add yrlngs
  trapRunDays$yrlngCatch<-ChiwYrlngCntLn$count[match(trapRunDays$EndDate,ChiwYrlngCntLn$date)]

  
  #add a column of scaled discharge (scaling based on the mean and sd in the efficiency trial data)
  trapRunDays$scaleDis<-(trapRunDays$dis- disch_Scale[1])/disch_Scale[2]
  
  #add year columns
  trapRunDays$year<-format(trapRunDays$EndDate,form="%Y")
  
  #add a column of day of year
  trapRunDays$DOY<-format(format(trapRunDays$EndDate,form="%j"))
  
  
  return(trapRunDays)
  
}#end of function




  Chiw_Effic<-Chiw_effic_dat_proc()
  Chiw_Catch_Ops<-Chiw_catch_dat_Proc(disch_Scale=Chiw_Effic$scaleDis)
  
  return(list(chiw_effic=Chiw_Effic$ChiwEffic,
              chiw_catch=Chiw_Catch_Ops))
  
}#end of function


