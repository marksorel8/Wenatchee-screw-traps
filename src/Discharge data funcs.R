Chiw_discharge_func<-function(){
  pkgTest("dataRetrieval")
  require("dataRetrieval")
  pkgTest("here")
  require("here")
  Chiw.flowDV<- readNWISdv(siteNumber="12456500",
                           "00060", "1990-01-01", "2018-12-31",statCd="00003")
  #munging
  Chiw.flowDV$date2<-as.Date(Chiw.flowDV$Date,format="%Y-%m-%d")
  
  Chiw.flowDV$Year<-format(Chiw.flowDV$date2,"%Y")
  Chiw.flowDV$week<-format(Chiw.flowDV$date2,"%W")
  Chiw.flowDV$month<-format(Chiw.flowDV$date2,"%m")
  
  Chiw.flowDV$day<-format(Chiw.flowDV$date2,"%d")
  
  Chiw.flowDV$doy<-format(Chiw.flowDV$date2,"%j")
  return(Chiw.flowDV)
}


#####################################################


Nason_White_Discharge_Func<-function(){
  pkgTest("here")
  ##Nason Flow
  #downloaded from https://fortress.wa.gov/ecy/eap/flows/station.asp?sta=45J070
  
  #function to load discharge data: "loc" is the folder location with the files (one for each yar), and "breakYr" is the file where the format of the data changes (related to the number of rows of header). 
  loadDis<-function(loc,breakYr){
    
    nasFiles<-list.files(loc)
    
    Nas<-data.frame()
    for ( i in breakYr:length(nasFiles)){
      test<-read.delim(paste(loc,nasFiles[i],sep="/"),skip=10,header=T,sep="")[,-4]
      Nas<-rbind(Nas,test)
    }
    
    colnames(Nas)<-c("date","discharge","data.code")
    #View(Nas)
    
    Nas$discharge<-as.character(Nas$discharge)
    Nas$data.code<-as.character(Nas$data.code)
    for( i in 1:(breakYr-1)){
      test<-read.fwf(paste(loc,nasFiles[i],sep="/"),
                     c(13,11,14,3),
                     skip=4,header=F)
      test<-test[,c(1,3,4)]
      colnames(test)<-colnames(Nas)
      Nas<-rbind(Nas,test)
    }
    
               
    Nas<-unique(Nas)
    return(Nas)
  }
  
  
  Nas2<-loadDis(loc=here("data","Nason and White","Nason Discharge"),breakYr = 4)
  
  
  ## White Flow
  #downloaded from https://fortress.wa.gov/ecy/eap/flows/station.asp?sta=45K090#block2
  
  White<-loadDis(loc=here("data","Nason and White","White Discharge"),breakYr = 13)
  
return(list(Nason_Dis=Nas2,White_Dis=White))
  
}
