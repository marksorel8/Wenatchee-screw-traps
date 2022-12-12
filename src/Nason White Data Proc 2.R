pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}


Nason_White_data_Func<-function(plot_length_day=FALSE){
  ###########################################################################
  #Packages
  pkgTest("here")
  pkgTest( "tidyverse" )
  
 here::i_am("src/Discharge data funcs.R")
  ###########################################################################  
  #Discharge data 
  
  #functiont to process Washington department of ecology discharge data
  source(here("src","Discharge data funcs.R"))
  
  
  #download discharge data from internet
  source(here("src","Discharge data funcs.R"))
  if(file.exists(here("nasWhiteDis.csv"))){
    load(here("nasWhiteDis.csv"))}else{
      disch_Dat<-Nason_White_Discharge_Func()
      save(disch_Dat,file=here("nasWhiteDis.csv"))  
    }
  
  #combine and munge discharge data
  ecol_disch<-rbind(((disch_Dat$Nason_Dis) %>% mutate(TrapNames ="Nason")),
        ((disch_Dat$White_Dis) %>% mutate(TrapNames ="White"))) %>% slice(-1) %>%  mutate(Rel.Date=as.Date(date,format="%m/%d/%Y")) %>% mutate(flow=as.numeric(flow)) %>% select(Rel.Date,TrapNames,flow)
         
 

########################################################################################################################################################
#Biological catch data
  
load(here("data","processed","cutoffs_and_props.Rdata")) #load delineation line

  #read in data
  bioD<-read.csv(here("data","Nason and White","Compiled Biodata.csv")) %>% 
  
  #change a column name that reads in weird sometimes
  rename(Trap=1) %>% 
  
  #add a column of stream name
  mutate(TrapNames=case_when(
    Trap%in%c("White 1.5 m","White 2.4 m")~"White",
    TRUE~"Nason")) %>% 
  
  #supstet to wild spring Chinook
  filter(Species=="Wild Spring Chinook") %>% 
  
  #drop precocial and adults
  filter(Stage!="A"&Stage!="PR") %>% 
  
  #add a reformated "Date" column
  mutate(Date2=as.Date(Date,format="%m/%d/%Y")) %>% 
  
  #add DOY
  mutate(DOY=as.numeric(format(Date2,form="%j"))) %>% 
  
  #make "length" numeric
  mutate(Length=as.numeric(as.character(Length))) %>% 
  
  #fix weird stage names
  mutate(Stage=gsub(" ","",Stage)) %>% 
  
  #add the stage based on delineation
  mutate(Stage_2=ifelse(DOY>179,"sub", #if DOY > 179 then subyearling
                        ifelse(!is.na(Length), #if Length available
                               ifelse(Length>=cutoffs_and_props[[1]]$y[(DOY-49)],"YCW","sub"), #assign age based on cutoff rule
                               ifelse(Stage!="", #if length not available but field call available
                                      ifelse(Stage%in%c("S","T"),"YCW","sub"), #call "smolts" and "transitionals" yearlings
                                      ifelse(Count>1,"sub", #if count >1 assume subyearling
                                             ifelse(cutoffs_and_props[[2]]$prop_lo_pred[(DOY-49)]>=.5,"YCW","sub")))))) %>% #otherwise assign based on the lifestage that was more common caught on that day  
  
 
  #get rid of fish captured in 2.4 m trap in years >=2017, when that trap was used only to catch fish for efficiency trials, not to estimate abundance.
 mutate(year=as.numeric(format(Date2,form="%Y"))) %>% 
  filter(Trap!="White 2.4 m"|year<2017) %>% 
  


  #drop a couple fish that were recaps or "morts in sculpin""
  filter(!TextualComment%in%c("RECAP","MORTALITY, RECAP INSIDE OF SCULPIN STOMACH")) %>% 
  
  ##some counts are 0 even though their is a fork length and wight, so I will convert those to 1. 
  mutate(Count=case_when(Count== 0 ~ as.integer(1),
                         TRUE ~ Count)) %>% 
  
         #sum daily catch by stream
    dplyr::group_by(Stage_2,TrapNames,Date2) %>%
    dplyr::summarize (catch=sum(Count)) %>%
    tidyr::spread(Stage_2,catch)
    


#########################################################################
  #load daily operations data
  OpsData<-read.csv(here("data","Nason and White","DailyOps (copy) 4132019.csv")) %>% 
  
  #fix names
  mutate(stream=case_when(Trap.Location=="Nason"~"Nason",TRUE~"White")) %>% 
  
  #reformat Date column
  mutate(Date2=as.Date(.[[1]],#date
                          format="%m/%d/%Y")) %>% 
  
  #Subset out days when trap was "pulled" or "stopped".
 filter(!(Trap.Stop.Time%in%c("Pulled","Stopped ")|Trap.Start.Time%in%c("Pulled ","Stopped "))) %>% 
    
  #make discharge data numeric
  mutate(Discharge=as.numeric(as.character(Discharge))) %>%  #make discharge numeric
  
  #fill in missing discharge value
  mutate(Discharge=case_when(is.na(Discharge)~488,
            TRUE~Discharge))  %>% 
  #
  
  #add year and DOY
  mutate(year=as.numeric(format(Date2,form="%Y")),
          DOY=as.numeric(format(Date2,form="%j")))


  
################################################################################
#efficiency trial data
  
  spChk_ET<-read.csv(here("data","Nason and White","Sorel_compiled_ET.csv")) %>% 
    
    #subset to spring chinook
    filter(Species%in%c(
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
      "Hat. Spring Chinook")) %>% droplevels() %>% 
    
    #simplify trap names
    mutate(Trap=as.character(Trap)) %>% 
    mutate(TrapNames=case_when(
      Trap%in%c("Nason","NASON","Nason Creek",
                "Nason ","Nason 1.5m")~"Nason",
      TRUE~"White")) %>% 
   
     #reformat release date to "Date" format
    mutate(Rel.Date=as.Date(Release.Date,format="%m/%d/%Y")) %>%  
    
    #add a column of recapture year and DOY. Day after release date would correspond to the day when (most) recaptures presumaby happen and the day of catch which each trial should presumably refer to.
    
    mutate(year=as.numeric(format(Rel.Date,format="%Y")),Rel.Date.plus.one=as.Date(Rel.Date+1), rel_plus_one_DOY=(as.numeric(format(Rel.Date,format="%j"))+1)) %>% 
    
    #add a column of simplified "lifestage" names
    mutate(Lifestage2=case_when(
      Lifestage%in%c("Parr","Subyearling")~"sub",
      Lifestage%in%c("Yearling","Smolt")~"YCW",
      Species%in%c("Wild Spring Chinook (0)","WILD SPRING CHINOOK (0)")~"sub",
      Species%in%c("Wild Spring Chinook (Y)","WILD SPRING CHINOOK (Y)","HATCHERY SPRING CHINOOK (Y)")~"YCW",
      TRUE~ifelse(rel_plus_one_DOY>175,"sub","YCW"))) %>% 
    
    #Make a column of simplified fish origin information
    mutate(origin = case_when(#combine age classes of suckers (Marie agreed to this)
      Species%in%c("HATCHERY SPRING CHINOOK (Y)","Hatchery Spring Chinook","Hat. Spring Chinook")~"Hatchery",
      Species%in%c("Wild Spring Chinook (Y)","Wild Spring Chinook (0)","Wild Spring Chinook"  ,"11W","WILD SPRING CHINOOK (0)","WILD SPRING CHINOOK (Y)")~"Wild",
      TRUE~"Unk")) %>% 
    
    
    
    #Drop rows where the "interrupted" column is yes
    filter(!Interruption %in% c("yes","YES","Yes")) %>% 
    
    #Drop a row where it says 0 fish were released
    filter(number.released > 0 ) %>%    
    
    
    #change discharge info from factor to numeric
    mutate(release.day.cfs=as.numeric((as.character(CFS.mean.release.day.)))) %>% 
    
    #add a new column of discharge data downloaded from the Wa Dept of Ecology guages
    left_join(ecol_disch,by=c("TrapNames","Rel.Date")) %>% 
    
    #fix outlier point 
    mutate(discharge2=case_when(flow==2370.0~release.day.cfs,
                                TRUE~flow)) %>% 
    
    # collapse by date and trap (a couple days had multiple releases, usually on different sides of the river)
    # sum daily catch by stream
    dplyr::group_by(Rel.Date,Rel.Date.plus.one,rel_plus_one_DOY,year,TrapNames,Lifestage2,discharge2) %>%
    dplyr::summarize (number.released=sum(number.released),number.recaptured=sum(number.recaptured)) %>% 
    
    
     droplevels() %>% 
    
    
    
    mutate(non_recaps=number.released-number.recaptured) 
  
  
  test<-glm(cbind(number.recaptured,non_recaps)~TrapNames*Lifestage2,data=spChk_ET,family="binomial" )
  
  summary(test)
  
  #########################################################################
  #Merge catch data and daily operations data  
  
  Ops_catch<-dplyr::left_join(OpsData,bioD,by=c("Date2","stream"="TrapNames" )) %>%
    rename("TrapNames" ="stream") %>% 
    #add a column of subs +fry called SBC
    dplyr::mutate(SBC=sub) %>%
    # Change NAs to 0 on days when the trap was operating but no fish were captured of a given lifestage
    dplyr::mutate_at(.vars=c("YCW","SBC"),~replace_na(.,0)) %>% 
    
    dplyr::left_join(spChk_ET,by=c("Date2"="Rel.Date.plus.one","TrapNames" )) %>% 
    
    #make columns with releases and recaps for trials that included subyearling or yearlings
    mutate(sub_rel=ifelse(Lifestage2=="sub",number.released,NA),
           sub_recap=ifelse(Lifestage2=="sub",number.recaptured,NA),
           yrlng_rel=ifelse(Lifestage2=="YCW",number.released,NA),
           yrlng_recap=ifelse(Lifestage2=="YCW",number.recaptured,NA)) %>%   
    #rename
    rename(count.sub=SBC,count.yrlng=YCW,Year=year.x,Disch.cfs=Discharge) %>% 
    
    select(TrapNames,sub_rel:yrlng_recap,count.sub,count.yrlng,DOY,Year,Disch.cfs) #select certain columns that will be used in modeling
  
  
  
  
  #Subset Nason
  trapDat_Nas<-droplevels(subset(Ops_catch,
                                    TrapNames=="Nason")) %>% 
  #add year factor
  droplevels() %>% mutate(year_factor=as.numeric(as.factor(Year))) %>% #year as factor
  # add a factor column for whether the trial was on or after the date when the Nason trap was moved.  
  mutate(moved=Year>= as.numeric(format(as.Date("2014-06-30"),format="%Y"))&
           DOY>= as.numeric(format(as.Date("2014-06-30"),format="%j")))
           
  
  #subset White
  trapDat_White <- droplevels(subset(Ops_catch,
                                      TrapNames=="White")) %>%
  #add year factor
  droplevels() %>% 
  mutate(year_factor = as.numeric(as.factor(Year))) #year as factor
  
  
  
  return(list(Nason=trapDat_Nas,White=trapDat_White))

}

