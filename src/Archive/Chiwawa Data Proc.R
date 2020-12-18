pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}#end of function

Chiw_dat_Proc<-function(){
# process efficiency trial data

  #read data
  pkgTest("here")
  library(here)
  ChiwEfficTrials<-read.csv(here("data","Chiwawa","Chiwawa_efficiency_trials_2.csv"))
  
  # change column names
  colnames(ChiwEfficTrials)<-c("modYear","lifeStage","Date","position","rel","recap","effic","disch.cfs","notes")
dim(ChiwEfficTrials)
  
#drop some rows that have comments suggesting that the data is not valid
  ChiwEfficTrials<-ChiwEfficTrials[!(ChiwEfficTrials$notes %in% c("Cant find any trial similar to flow, sample size in 2007?","Cant find anything close to these groups","crossed out on datasheet","Nothing close in records","YCW from Lake Trap trial","Only used YCW. From Lake Trap trial, not sure why it was used in Chiwawa Model - JW")),]
  
  
  #format dates
  ChiwEfficTrials$Date<-as.character(ChiwEfficTrials$Date)
  
  ChiwEfficTrials$Date<-as.Date(ChiwEfficTrials$Date,format = ifelse(grepl("-",ChiwEfficTrials$Date),"%d-%b-%y","%m/%d/%Y"))
  
    #reorder the rows so that the ones with actual dates and knowcome first, that way when we remove duplicates we dont remove the ones with full dates
  ChiwEfficTrials<-ChiwEfficTrials[order(nchar(as.character(ChiwEfficTrials$lifeStage)),decreasing = F),]
  
  
  #drop duplicate rows
  ChiwEffic<-ChiwEfficTrials[!duplicated(ChiwEfficTrials$Date),]
  
  ChiwEffic[ChiwEffic$notes %in% c("YCW only","Only used YCW"),"lifeStage"]<-"YCW"
  
  ChiwEffic[ChiwEffic$notes %in% c("SBC only","Most likely SBC based on date"),"lifeStage"]<-"SBC"
  
  
  #make a new "position" column where "low flow" is changed to "upper", because only a few data points for "low flow"
  ChiwEffic$position2<-ChiwEffic$position

  ChiwEffic$position2[ChiwEffic$position=="Low Flow"]<-"Upper"
  
  ChiwEffic<-droplevels(ChiwEffic)
  
  #add day, week, and year
  ChiwEffic$day<-as.numeric(format(ChiwEffic$Date,form="%j"))
  ChiwEffic$week<-ceiling(ChiwEffic$day/7)
  ChiwEffic$year<-as.numeric(format(ChiwEffic$Date,form="%Y"))
  
  #meand and standard deviation of discharge
  disMean<-attributes(scale(ChiwEffic$disch.cfs))$`scaled:center`
  
  disScale<-attributes(scale(ChiwEffic$disch.cfs))$`scaled:scale`
  


# process dialy catch and operations data
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

  
  #fill in the few missing discharges with data obtained with dataRetrieval package
  
  #download discharge data from internet
  source(here("src","Discharge data funcs.R"))
  chiwDis<-Chiw_discharge_func()
  
  #fill in missing discharges with downloaded data
  trapRunDays$dis<-trapRunDays$Mean.discharge..CFS.
  trapRunDays$dis[is.na(trapRunDays$Mean.discharge..CFS.)]<-
    chiwDis$X_00060_00003[match(trapRunDays$EndDate[is.na(trapRunDays$Mean.discharge..CFS.)],chiwDis$date-1)]
  
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

  
  
  #*********   Guess at missing trap positions   **********
  #*********   should be estimated in model in future   ***
  tapply(ChiwEffic$disch.cfs,ChiwEffic$position2,summary)
  trapRunDays$Position2[trapRunDays$Position2==""& trapRunDays$dis<=475]<-"Upper"
  trapRunDays$Position2[trapRunDays$Position2==""]<-"Lower"
  trapRunDays$Position2<-factor(trapRunDays$Position2,levels=c("Lower","Upper"))
  
  trapRunDays<-droplevels(trapRunDays)
  
  
  
  #add year columns
  trapRunDays$year<-format(trapRunDays$EndDate,form="%Y")
  
  #add a column of day of year
  trapRunDays$DOY<-format(format(trapRunDays$EndDate,form="%j"))
  
  
  #data for design matrix for efficiency model
  N_catch_obs<-nrow(trapRunDays)
  
  catch_ops_effic<-merge(trapRunDays,ChiwEffic ,by.x=c("EndDate"),by.y = c("Date"),all=TRUE,sort=FALSE,suffixes = c("",".y"))
  
  
  #add rows for efficiency trial dates with no catch data ( not sure how this is possible because it seems like you would need catch for , but it is) or for efficiency trials from a different lifestage.
  
  catch_ops_effic[(N_catch_obs+1):nrow(catch_ops_effic),"Position2" ]<-catch_ops_effic[(N_catch_obs+1):nrow(catch_ops_effic),"position2" ]#position2
  catch_ops_effic[(N_catch_obs+1):nrow(catch_ops_effic),"dis" ]<-catch_ops_effic[(N_catch_obs+1):nrow(catch_ops_effic),"disch.cfs" ]#dis
  catch_ops_effic[(N_catch_obs+1):nrow(catch_ops_effic),"year" ]<-catch_ops_effic[(N_catch_obs+1):nrow(catch_ops_effic),"year.y" ]#year
  catch_ops_effic[(N_catch_obs+1):nrow(catch_ops_effic),"DOY" ]<-catch_ops_effic[(N_catch_obs+1):nrow(catch_ops_effic),"day" ]#DOY
  
  #scale discharge
  catch_ops_effic$log_scale_Dis<-scale(log(catch_ops_effic$dis))
  
  return(catch_ops_effic)

  }#end of function


