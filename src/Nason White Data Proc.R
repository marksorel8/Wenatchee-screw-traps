Nason_White_data_Func<-function(){


pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}



White_Nason_Effic_Func<-function(){
  pkgTest("here")
  
  dat<-read.csv(here("data","Nason and White","Sorel_compiled_ET.csv"))
  
  # all the different names that refure to sping Chinook
  spChk_names<-c(
    "11W",
    "Chinook",
    "Spring Chinook (unknown r/t)" ,
    "Wild Spring Chinook" ,
    "Wild Spring Chinook (0)" ,
    "WILD SPRING CHINOOK (0)",
    "Wild Spring Chinook (Y)",
    "WILD SPRING CHINOOK (Y)" ,
    "HATCHERY SPRING CHINOOK (Y)",
    "Hatchery Spring Chinook",
    "Hat. Spring Chinook")
  
  #subset to spring chinook
  spChk_ET<-droplevels(dat[!is.na(match(dat$Species,spChk_names)),])
  
  
  #Names refering to nason creek
  nasNames<-c("Nason","NASON","Nason Creek",
              "Nason ","Nason 1.5m")
  
  #simplify trap names
  spChk_ET$TrapNames<-ifelse(!is.na(match(spChk_ET$Trap,nasNames)),"Nason","White")
  
  #Names referinging subyearlings
  subNames<-c("Parr","Subyearling")
  
  #add a column of simplified "lifestage" names
  spChk_ET$Lifestage2<-ifelse(!is.na(match(spChk_ET$Lifestage,subNames)),"Subyearling",ifelse(spChk_ET$Lifestage=="","Unk", "Yearling"))
  
  #making column of simplified "lifestage names continued
  subNames<-c("Wild Spring Chinook (0)","WILD SPRING CHINOOK (0)")
  yrlngNames<-c("Wild Spring Chinook (Y)","WILD SPRING CHINOOK (Y)","HATCHERY SPRING CHINOOK (Y)")
  
  spChk_ET$Lifestage2[!is.na(match(spChk_ET$Species,subNames))]<-"Subyearling"
  
  spChk_ET$Lifestage2[!is.na(match(spChk_ET$Species,yrlngNames))]<-"Yearling"
  
  
  #Names refuring to hatchery-origin fish
  hatNames<-c("HATCHERY SPRING CHINOOK (Y)","Hatchery Spring Chinook","Hat. Spring Chinook")
  wilNames<-c("Wild Spring Chinook (Y)","Wild Spring Chinook (0)","Wild Spring Chinook"  ,"11W","WILD SPRING CHINOOK (0)","WILD SPRING CHINOOK (Y)")
  
  #Make a column of simplified fish origin information
  spChk_ET$origin<-ifelse(!is.na(match(spChk_ET$Species,hatNames)),"Hatchery",ifelse(!is.na(match(spChk_ET$Species,wilNames)),"Wild","Unk"))
  
  
  #change discharge info from factor to numeric
  spChk_ET$release.day.cfs<-as.numeric(as.character(spChk_ET$CFS.mean.release.day.))
  
  
  #add a new column of discharge data downloaded from the Wa Dept of Ecology guages
  
  #import data
  source(here("src","Discharge data funcs.R"))
  disch_Dat<-Nason_White_Discharge_Func()
  
  
  #reformat release date to "Date" format
  spChk_ET$Rel.Date<-as.Date(spChk_ET$Release.Date,format="%m/%d/%Y")
  
  #add Nason discharge data 
  spChk_ET$discharge2[spChk_ET$TrapNames=="Nason"]<-disch_Dat$Nason_Dis$discharge[match(spChk_ET$Rel.Date[spChk_ET$TrapNames=="Nason"],as.Date(disch_Dat$Nason_Dis$date,format="%m/%d/%Y"))]
  #add white discharge data
  spChk_ET$discharge2[spChk_ET$TrapNames=="White"]<-disch_Dat$White_Dis$discharge[match(spChk_ET$Rel.Date[spChk_ET$TrapNames=="White"],as.Date(disch_Dat$White_Dis$date,format="%m/%d/%Y"))]
  
  #change form character to numeric
  spChk_ET$discharge2<-as.numeric(  spChk_ET$discharge2)
  
  #fix weird outliner
  spChk_ET$discharge2[spChk_ET$discharge2==2370.0]<-spChk_ET$release.day.cfs[spChk_ET$discharge2==2370.0]
  

  #add a column of proportion recaputed
  spChk_ET$p<-spChk_ET$number.recaptured/spChk_ET$number.released
  
  #add a column of release year and DOY
  spChk_ET$year<-format(spChk_ET$Rel.Date,format="%Y")
  spChk_ET$DOY<-as.numeric(format(spChk_ET$Rel.Date,format="%j"))
  
  #Drop rown where the "interrupted" column is year
  interpted<-c("yes","YES","Yes")
  spChk_ET<-droplevels(subset(spChk_ET,is.na(match(spChk_ET$Interruption,interpted))))
  
#Drop a row where it says 0 fish were released
    spChk_ET<-subset(spChk_ET,p<Inf)

#Subset to just Nason 
 Nas_spChk_ET<-droplevels(subset(spChk_ET,TrapNames=="Nason"))
  
#Make year a numeric
  Nas_spChk_ET$year<-as.numeric(Nas_spChk_ET$year)
#Make a column for a "year factor" starting at 0
Nas_spChk_ET$yearFac<-Nas_spChk_ET$year-min(Nas_spChk_ET$year)
  
#make a column of day number from a continuous sequence starting at Jan 1 on the first year
Nas_spChk_ET$daySeq<-Nas_spChk_ET$yearFac*365+Nas_spChk_ET$DOY
  
#add a factor column for whether the trial was on or after the date when the Nason trap was moved.   
Nas_spChk_ET$moved<-as.factor(ifelse(Nas_spChk_ET$Rel.Date>=as.Date("2014-06-30"),1,0))
  
#Subset White
White_spChk_ET<-droplevels(subset(spChk_ET,TrapNames=="White"))
  

return(list(Nason_ET=Nas_spChk_ET,
            White_ET=White_spChk_ET))

}



########################################################################################################################################################





Nason_White_Catch_Ops_Dat_func<-function(){
  
  #read in data
  bioD<-read.csv(here("data","Nason and White","Compiled Biodata.csv"))
  
  #change a column name that reads in weird sometimes
  colnames(bioD)[1]<-"Trap"
  
  #add a column of stream name
  bioD$stream<-ifelse(substr(bioD$Trap,1,1)=="N","Nason","White")
  
  #supstet to wild spring Chinook
  wSpC<-droplevels(subset(bioD,Species=="Wild Spring Chinook"))
  
  #drop precocial and adults
  wSpC<-droplevels(subset(wSpC,Stage!="A"&Stage!="PR"))
  
  #add a reformated "Date" column
  wSpC$Date2<-as.Date(wSpC$Date,format="%m/%d/%Y")
  
  #add DOY
  wSpC$DOY<-format(wSpC$Date2,format="%j")
  wSpC$DOY<-as.numeric(wSpC$DOY)
  
  #make "length" numeric
  wSpC$Length<-as.numeric(as.character(wSpC$Length))

  
  #Have to assign lifestages based on length and DOY at capture

  #White River
  plot(wSpC$DOY,as.numeric(as.character(wSpC$Length)),pch=19,cex=.4,type="n",ylab="Length",xlab='DOY',main="White River")
  
  points(wSpC$DOY[wSpC$stream!="Nason"],wSpC$Length[wSpC$stream!="Nason"], type="p",col=rgb(.1,.1,.1,.3),pch=19,cex=.4)
 
   #this line is what I came up with, where everythng below is a subyearling and everything above is a yearling
  segments(c(0,110,200),c(55,55,100),c(110,200,350),
           c(55,100,137.5),col="red")
  abline(h=50,col="red")
  
  #Nason Creek
  plot(wSpC$DOY,as.numeric(as.character(wSpC$Length)),pch=19,cex=.4,type="n",ylab="Length",xlab='DOY',main="Nason Creek")
  
  points(wSpC$DOY[wSpC$stream=="Nason"],wSpC$Length[wSpC$stream=="Nason"], type="p",col=rgb(.1,.1,.1,.3),pch=19,cex=.4) 
  
  
  #this line is what I came up with, where everythng below is a subyearling and everything above is a yearling
  segments(c(0,110,200),c(55,55,100),c(110,200,350),
           c(55,100,137.5),col="red")
  
  abline(h=50,col="red")
  
  
  
  #add the stage based on my delineation
  wSpC$Stage_2<-ifelse(wSpC$Length<=50,"fry",ifelse((wSpC$DOY<110 & wSpC$Length<=55) |
                         ((wSpC$DOY>=110 & wSpC$DOY<200)&((wSpC$Length/wSpC$DOY)<.5))|
                         ((wSpC$DOY>=200) &((wSpC$Length-(.25*wSpC$DOY))<50)),
                       "sub","yrlng"))
  
  
  #get rid of fish captured in 2.4 m trap in years >=2017, when that trap was used only to catch fish for efficiency trials, not to estimate abundance.
  wSpC$year<-as.numeric(format(wSpC$Date2,form="%Y"))
  wSpC<-subset(wSpC,Trap!="White 2.4 m"|year<2017)
  
  
  ##"guessing" lifestage of fish without lengths
  
  #Fish I couldn't assign to lifestage
  NAs<-subset(wSpC,is.na(Stage_2))
  
  
  #For fish that didn’t have length information, assigned age/lifestage based on data in “Stage” or “ConditionalComment” column. 
  wSpC$Stage_2[is.na(wSpC$Stage_2)&(wSpC$Stage=="F"|wSpC$Stage=="      F")]<-"fry"
  wSpC$Stage_2[is.na(wSpC$Stage_2)&(wSpC$Stage=="P"|wSpC$Stage=="   P")]<-"sub"
  
  wSpC$Stage_2[is.na(wSpC$Stage_2)&wSpC$Stage=="T"]<-"yrlng"
  wSpC$Stage_2[is.na(wSpC$Stage_2)&(wSpC$Stage=="S"|wSpC$Stage=="      S")]<-"yrlng"
  
  wSpC$Stage_2[is.na(wSpC$Stage_2)&
                 wSpC$ConditionalComment=="Y"]<-"yrlng"
  
  wSpC$Stage_2[is.na(wSpC$Stage_2)&
                 wSpC$Weight>100]<-"yrlng"
  
  wSpC$Stage_2[is.na(wSpC$Stage_2)&
                 wSpC$ConditionalComment=="0"]<-"sub"
  
 #drop a couple fish that were recaps or "morts in sculpin" Haha
   wSpC<-wSpC[!is.na(wSpC$Stage_2)|wSpC$TextualComment!="RECAP", ]
  
  wSpC<-wSpC[wSpC$TextualComment!="MORTALITY, RECAP INSIDE OF SCULPIN STOMACH", ]
  
  #Those with no information in those columns were assigned as subyearlings if they were tagged on DOY> 200. If they were in a row with multiple fish they were assumed to be fry,assuming that smaller fish (fry?) were not measured and given individual rows). 
  wSpC$Stage_2[is.na(wSpC$Stage_2)&
                 wSpC$DOY>=200]<-"sub"
  
  wSpC$Stage_2[is.na(wSpC$Stage_2)&
                 wSpC$Count>=10]<-"fry"
  
  
    
  #summing daily counts
  
  
  ##some counts are 0 even though their is a fork length and wight, so I will convert those to 1. 
  wSpC$Count[wSpC$Count==0]<-1
  
  
  
  ##subsetting lifestages for processing
  
  fry<-subset(wSpC,Stage_2=="fry")
  subs<-subset(wSpC,Stage_2=="sub")
  yrlngs<-subset(wSpC,Stage_2=="yrlng")
  
  
  
  #Fish I couldn't assign to lifestage
  NAs<-subset(wSpC,is.na(Stage_2))


  
  DFbio<-t(tapply(wSpC$Count, wSpC[,c("stream","Date2")], sum))
  
  dailyFry<-t(tapply(fry$Count, fry[,c("stream","Date2")], sum))
  
  dailySubs<-t(tapply(subs$Count, subs[,c("stream","Date2")], sum))
  
  dailyYrlngs<-t(tapply(yrlngs$Count, yrlngs[,c("stream","Date2")],sum))
  
  #changing format to long
  pkgTest("reshape2")
  DFbio2<-melt(DFbio)
  dailyFry2<-melt(dailyFry)
  dailySubs2<-melt(dailySubs)
  dailyYrlngs2<-melt(dailyYrlngs)
  
  #reformat date columns
  DFbio2$Date2<-as.Date(DFbio2$Date2)
  dailyFry2$Date2<-as.Date(dailyFry2$Date2)
  dailySubs2$Date2<-as.Date(dailySubs2$Date2)
  dailyYrlngs2$Date2<-as.Date(dailyYrlngs2$Date2)
  
  
  #load daily operations data
  OpsData<-read.csv(here("data","Nason and White","DailyOps (copy) 4132019.csv"))
  
  #fix names
  OpsData$stream<-ifelse(OpsData$Trap.Location=="Nason","Nason","White")
  
  #reformat Date column
  OpsData$Date2<-as.Date( OpsData[,1],#date
                          format="%m/%d/%Y")
  
  
  #add catch data to daily ops data
  Ops_catch<- merge(OpsData,DFbio2,by.x=c(17,16),by.y=1:2, all.x = TRUE,all.y = FALSE)
  
  Ops_catch<- merge(Ops_catch,dailyFry2,by.x=1:2,by.y=1:2, all.x = TRUE,all.y = FALSE)
  
  Ops_catch<- merge(Ops_catch,dailySubs2,by.x=1:2,by.y=1:2, all.x = TRUE,all.y = FALSE)
  
  Ops_catch<- merge(Ops_catch,dailyYrlngs2,by.x=1:2,by.y=1:2, all.x = TRUE,all.y = FALSE)
  
  colnames(Ops_catch)[18:21]<-c("tot_catch","Fry_catch","Sub_catch","Yrlng_catch")
  
  # Change NAs to 0 on days when the trap was operating but no fish were captured of a given lifestage
  
  Ops_catch[,18:21][is.na(Ops_catch[,18:21])]<-0
  
  
  
  # some fish from bio data are excluded because there is no operations data from those days (i.e. all of 2003 in Nason). Look at the excluded fish
  test<-DFbio2[is.na(match(paste0(DFbio2$Date2,DFbio2$stream),
                           paste0(OpsData$Date2,
                                  OpsData$stream))),]
  
  #View(test)
  
  
  #Subset out days when trap was "pulled" or "stopped".
  
  # sort(table(Ops_catch$Trap.Stop.Time))
  # names(sort(table(Ops_catch$Trap.Stop.Time)))
  # sort(table(as.character(Ops_catch$Trap.Start.Time)))
  # names(sort(table(as.character(Ops_catch$Trap.Start.Time))))
   
  # indicators of incomplete data
  notTrap<-c("Pulled ","Stopped ")
  
  #drop incomplete data points (days)
  trapDatSub<-subset(Ops_catch,
                     is.na(match(Ops_catch$Trap.Stop.Time,notTrap ))&
                       is.na(match(Ops_catch$Trap.Start.Time,notTrap )) )
  
  
  #check for missing discharge data
  
  trapDatSub$Discharge2<-as.numeric(as.character(trapDatSub$Discharge)) #make discharge numeric
  trapDatSub$date2<-as.Date(trapDatSub$Date2,format="%m/%d/%Y")  #format date column
  trapDatSub<-trapDatSub[order(trapDatSub$date2),] # sort by date
  trapDatSub<-trapDatSub[order(trapDatSub$Trap.Location),] # sort by stream
  which(is.na(trapDatSub$Discharge2)) # location of NAs
  
  trapDatSub$Discharge2[which(is.na(trapDatSub$Discharge2))]<-
    trapDatSub$Discharge2[(which(is.na(trapDatSub$Discharge2)))-1] #replace NA with value form previous day
  
  #add year and DOY
  trapDatSub$year<-format(trapDatSub$date2,form="%Y")
  trapDatSub$DOY <-format(trapDatSub$date2,form="%j")
  
  
  # for ( i in unique(trapDatSub$Trap.Location)){
  #   streamSub<-subset(trapDatSub,Trap.Location==i)
  #   for ( j in unique(streamSub$year)){
  #     year_<-subset(streamSub,year==j)
  #     
  #     plot(year_$DOY,year_$Sub_catch,
  #          main=paste("sub",j,i))
  #     
  #     plot(year_$DOY,year_$Yrlng_catch,
  #          main=paste("yrlng",j,i))
  #     
  #   }
  # }
  
  #Subset Nason
  trapDatSub_Nas<-droplevels(subset(trapDatSub,
                                    stream=="Nason"))
  #subset White
  trapDatSub_White<-droplevels(subset(trapDatSub,
                                      stream=="White"))

  return(list(Nason_catch=trapDatSub_Nas,White_catch=trapDatSub_White))
  
}

effic_dat<-White_Nason_Effic_Func()
catch_ops_dat<-Nason_White_Catch_Ops_Dat_func()

return(list(effic_dat=effic_dat,catch_ops_dat=catch_ops_dat))
}
