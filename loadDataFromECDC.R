#these libraries are necessary

library(readxl)
library(httr)
library(tidyr)
library(ggplot2)
library (deSolve)
library(dplyr)
library(scales)
library (readr)

#Load the newest data
  #create the URL where the dataset is stored with automatic updates every day

  #Load data from Hopkins 
  urlfile="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv"
  mydata<-read_csv(url(urlfile))  


  url <- paste("https://www.ecdc.europa.eu/sites/default/files/documents/COVID-19-geographic-disbtribution-worldwide-",format(Sys.time(), "%Y-%m-%d"), ".xlsx", sep = "")
  #url <- paste("https://www.ecdc.europa.eu/sites/default/files/documents/COVID-19-geographic-disbtribution-worldwide-","2020-03-20.xlsx", sep = "")
  #download the dataset from the website to a local temporary file
  GET(url, authenticate(":", ":", type="ntlm"), write_disk(tf <- tempfile(fileext = ".xlsx")))
  #read the Dataset sheet into “R”
  OriginalData <- read_excel(tf)
  
  OriginalData$DateRep <- as.Date(OriginalData$DateRep)
  OriginalData$GeoId <- factor(OriginalData$GeoId)

#make total counts out of new cases
  OriginalData <- OriginalData %>% group_by(GeoId) %>% arrange(DateRep) %>% mutate(sumCases=cumsum(Cases)) %>% mutate(sumDeaths=cumsum(Deaths))

  #Add the number of people with Corona over 60 from manual data collection
  CoronaPositiveOver60 <- data.frame(
    DateRep = as.Date(
      c("2020-03-21",
        "2020-03-20",
        "2020-03-19",
        "2020-03-18",
        "2020-03-17",
        "2020-03-16",
        "2020-03-15",
        "2020-03-14",
        "2020-03-13",
        "2020-03-12",
        "2020-03-11"))
    ,
    GeoId="DE",
    nCoronaPositiveOver60 = c(2854,2342,1800,1337,1181,900,703,525,416,293,142),
    stringsAsFactors=FALSE
  )
  
  #Merge the two datafiles
  
  CombinedData <- left_join(OriginalData,CoronaPositiveOver60)  
  
#Plot of south korea, italy, spain, US and german data
  SelectedCountries <- filter(CombinedData, GeoId %in% c("DE","CH","IT","ES","US","KR"))
  SelectedCountries <- filter(SelectedCountries, DateRep >= as.Date("2020-03-01"))
  myplot_step1 <- ggplot(data=SelectedCountries,aes(x=DateRep, y=sumCases))+
    geom_line(aes(colour=factor(GeoId)),size=1)+
    geom_point(aes(colour=factor(GeoId)),size=1)+
    scale_y_log10(breaks = log_breaks())
  myplot_step1
  myplot_step2 <- ggplot(data=SelectedCountries,aes(x=sumCases, y=sumDeaths))+
    geom_line(aes(colour=factor(GeoId)),size=1)+
    geom_point(aes(colour=factor(GeoId)),size=1)+
    scale_x_log10(breaks = log_breaks())+
    theme_bw()
    #facet_wrap(~ factor(GeoId))
  myplot_step3 <- ggplot(data=SelectedCountries,aes(x=DateRep, y=sumCases))+
    geom_line(aes(colour=factor(GeoId)),size=1)+
    geom_point(size=1)+
    geom_line(aes(x=DateRep, y=nCoronaPositiveOver60,colour=factor(GeoId)),size=1,shape=1)+
    geom_point(aes(x=DateRep, y=nCoronaPositiveOver60),size=1,shape=1)+
    scale_y_log10(breaks = log_breaks())+
    facet_wrap(~ factor(GeoId))+
    theme_bw()+
    theme(legend.position ="none")
    
  
#next step:find out numbers of infected >60 in italy

#https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30627-9/fulltext
#Specifically, the percentage of patients admitted to intensive care units reported daily in Italy, from March 1, up until March 11, was consistently between 9% and 11% of patients who were actively infected.

odes_TurnoverModel_3_stimulation_of_buildup <- function (t, A, p) {  # ODE system 
  
  
  ##assign the current values
  Dose_Central = A[1]
  Drug_Central = A[2]
  Drug_Peri = A[3]
  Response = A[4]
  
  ##assign the parameters
  Central        <- p$Central
  Peripheral     <- p$Peripheral
  
  ka_Central <- p$ka
  ke_Central <- p$keC
  
  kcp <-p$kcp
  
  kin <- p$Rss*p$kout
  #kin <- p$kin
  kout <- p$kout
  #kpc <-p$kpc
  
  ##reactions
  
  ReactionFlux1 = ka_Central*Dose_Central
  ReactionFlux2 = (ke_Central*Drug_Central)*Central
  ReactionFlux4 = (kcp*Drug_Central)*Central-(kcp*Drug_Peri)*Peripheral
  
  
  
  ## right hand side
  
  ddtDose_Central = -ReactionFlux1
  ddtDrug_Central =  1/Central*(ReactionFlux1 - ReactionFlux2 - ReactionFlux4)
  ddtDrug_Peripheral = 1/Peripheral*(ReactionFlux4)
  ddtResponse = kin*(1+(p$Emax_turnover*Drug_Central/(p$EC50 +Drug_Central))) - kout*Response
  
  deriva <- list ( c (  ddtDose_Central, ddtDrug_Central, ddtDrug_Peripheral, ddtResponse) )
  
  #print(deriva)
  
  return ( deriva )
}
