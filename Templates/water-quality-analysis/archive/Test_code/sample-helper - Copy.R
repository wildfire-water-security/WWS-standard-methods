library(dplyr)
library(stringr)
library(lubridate)
library(writexl)



###for ISCO unknown ##

if (samp_type=="IS"){

sn <- sum(nsamp) #total sample amount from the vector that was filled out

#functions for data labels 

#first function - diffrent for known
data_label <- function(){
  sampleinfo <- data.frame(site=site,nsamp=nsamp,interval=interval,study=study,samp_type=samp_type,blanktime=blanktime, sn=sn, rn=rn)
  data_info <- function(sampleinfo){ #this function creates 
    return(data.frame(
      site =rep(sampleinfo$site, times = sampleinfo$nsamp), #repeats site for number of samples per site
      blank_time=sampleinfo$blanktime, #inserts blank time format 
      botnum = sprintf("%02d",rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp]), #bottle numbers 1-24 repeated for number of samples
      study=sampleinfo$study,samp_type=sampleinfo$samp_type
    ))}
  labels<-lapply(1:nrow(sampleinfo),function(x){data_info(sampleinfo[x,])})%>% #using data_info function to make data.frame
    bind_rows()%>%
    mutate(
      randomn=sample(rn:sn,sn, replace=F), #add random numbers
      sample_name = paste0(study, "_", "R", sprintf("%04d",randomn)), 
      field_sample_ID_blank= paste0(study, "_", site, "_", samp_type,"_",blanktime))
  return(labels)

}



#second function- diffrent for known
if(pt2==TRUE){
sampletimeinfo <- data.frame(nsamp=nsamp,interval=interval,updated_start=updated_start,DOC=DOC,DOM=DOM,CCAL=CCAL,Levoglucosan=Levoglucosan,Excess=Excess)   

#first part 
update_time <- function(){
  site_label_blank <- function(sampletimeinfo){
    timestart <- as.POSIXct(sampletimeinfo$updated_start)
    timelist <- seq(timestart, by=paste0(sampletimeinfo$interval, " hours"), length = sampletimeinfo$nsamp)
    return(
      data.frame(time_PST=as.character(timelist),datetime_PST=format(timelist,"%Y%m%d%_%H%M")))} #might need to split this for the new template
  
  data_labels <- lapply(1:nrow(sampletimeinfo),function(x){site_label_blank(sampletimeinfo[x,])}) %>% #makes the frame
    bind_rows()%>%
    mutate(blank_labels)%>%
    mutate(field_sample_ID= paste0(study, "_", site, "_",samp_type,"_",datetime_PST))%>%
    arrange(sample_name)
  return(data_labels)
}

lab_label <- function(){
  analysis_opt <- c("CN", "AQ", "NU", "LV", "EX")
  analysis_choice <- c(DOC, DOM, CCAL, Levoglucosan, Excess)
  analysiscodes <- analysis_opt[analysis_choice] #analytes used
  analysisn <- length(analysiscodes) #how many analytes there are
  datacodes<- data.frame(code = rep(analysiscodes,times=sn)) %>% #add labcodes
    mutate(analysis = case_when(
      code == "CN" ~ "Shimadzu",
      code == "AQ" ~ "Aqualog",
      code == "NU" ~ "CCAL",
      code == "LV" ~ "Levoglucosan",
      code == "EX" ~ "Excess",
      TRUE ~ NA_character_))
  
  
  data_labels <- update_time()%>% #putting this here so code doesn't have to get confused
    mutate(blank_labels)%>% #add on your data frame that had blank labeles
    slice(rep(1:n(), each = analysisn))%>% #multiply for each analyte
    mutate(datacodes, #add in sample codes
         lab_sample_ID =
           paste0(field_sample_ID, "_",
                  datacodes$code )) #add in lab labels 
  #add the data codes 
  #NEED TO ADD ENDING CODE TO SAMPLE ID

  return(data_labels) #add print here
} #lab_label
} #if pt 2
} #if ISCO



















####For grab samples ##

if (samp_type=="GB"){
  
sn <- nrow(site) #how many sites there are  


sampleinfo<- data.frame(site=site,study=study,samp_type=samp_type,blanktime=blanktime,sn=sn,rn=rn)

data_label <- function(){ #this will add blank labels in 
  datablank<- data.frame(site)%>%
    mutate(study=study,
           samp_type=samp_type,
           blanktime=blanktime,
           randomn= sample(rn:sn,sn, replace=F),
           date_text_PST = "",
           time_text_PST ="") 
  blank_labels <- datablank %>%
    mutate(sample_name=paste0( study, "_", "R", sprintf("%04d", randomn)), #adds sample name
           field_sample_ID_blank = paste0(study, "_",
                                          site, "_",
                                          samp_type,"_",
                                          blanktime))#adds blank sample ID
  data_label <- blank_labels %>% select(
    sample_name,field_sample_ID_blank,randomn,study,samp_type,site,date_text_PST,time_text_PST) #reorder so it's easier for data entry
  return(data_label)
}  
  


if(pt2==TRUE){

update_time<- function(){
  updated_labels<-blank_labels %>%
    mutate(datetime_PST=
             paste0(sprintf("%06d",blank_labels$date_text_PST),"_",   #manual find 
                    sprintf("%04d",blank_labels$time_text_PST)),
           field_sample_ID =
             paste0(study, "_",
                    site, "_",
                    samp_type,"_",
                    datetime_PST))%>%
    arrange(sample_name)
  

  return(updated_labels)
}
  
  

  

  

} #if pt 2 

}#if GB 



#this might be the same for ISCO 
lab_label<- function(){  
  analysis_opt <- c("CN", "AQ", "NU", "EX", "LV")
  analysis_choice <- c(DOC, DOM, CCAL, Excess, Levoglucosan)
  
  analysiscodes <- analysis_opt[analysis_choice]
  analysiscodes
  
  analysisn <- length(analysiscodes) #how many analytes there are
  samples <- analysisn*sn #how many total filtered samples
  
  datacodes<- data.frame(code = rep(analysiscodes,times=sn)) %>% #add labcodes
    mutate(analysis = case_when(
      code == "CN" ~ "Shimadzu",
      code == "AQ" ~ "Aqualog",
      code == "NU" ~ "CCAL",
      code == "LV" ~ "Levoglucosan",
      code == "EX" ~ "Excess",
      TRUE ~ NA_character_
    ))
  datacodes
  #repeat samples for each analysis
  
  data_labels <- update_time()%>% #putting this here so code doesn't have to get confused <- updated_labels %>% #put in so code isn't confused 
    slice(rep(1:n(), each = analysisn)) %>%
    mutate(datacodes,
           lab_sample_ID = #add in lab labels 
             paste0(field_sample_ID, "_",
                    datacodes$code ))
  
  return(data_labels)
}