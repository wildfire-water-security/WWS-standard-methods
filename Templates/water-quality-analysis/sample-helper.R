library(dplyr)
library(stringr)
library(lubridate)
library(writexl)
library(openxlsx2)




###for ISCO unknown ##





if (samp_type=="IS"){

sn <- sum(nsamp) #total sample amount from the vector that was filled out

#functions for data labels 
if (knowntime==TRUE){ #blank data frame
  field_label <- function(){ #makes data frame
    sampleinfo <- data.frame(site=site,nsamp=nsamp,interval=interval,study=study,samp_type=samp_type,start=start, sn=sn, rn=rn)
    data_info <- function(sampleinfo){ #this function creates
      sitelist <- rep(sampleinfo$site, times = sampleinfo$nsamp)
      timestart <- as.POSIXct(sampleinfo$start)
      timelist <- seq(timestart, by=paste0(sampleinfo$interval, " hours"), length = sampleinfo$nsamp)
      botnum <- rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp] 
      data <-data.frame(site =sitelist, timelist=timelist, 
                        botnum = sprintf("%02d",botnum), 
                        study=sampleinfo$study,
                        samp_type=sampleinfo$samp_type,
                        datetime_PST=format(timelist,"%Y%m%d%_%H%M"))
      return(data)
    }
    field_data<-lapply(1:nrow(sampleinfo),function(x){data_info(sampleinfo[x,])})%>%
      bind_rows()%>%
      mutate(
        randomn=sample(rn:sn,sn, replace=F),
        sample_name = paste0(study, "_", "R", sprintf("%04d",randomn)),
        field_sample_ID = paste0(study, "_", site, "_", samp_type,"_",datetime_PST)
      )
    return(field_data)
  }

} #field_label #knowntime=TRUE

else if(knowntime==FALSE){ #blank data frame
  #first function - diffrent for known
  field_label <- function(){
  sampleinfo <- data.frame(site=site,nsamp=nsamp,interval=interval,study=study,samp_type=samp_type,blanktime=blanktime, sn=sn, rn=rn)
  data_info <- function(sampleinfo){ #this function creates 
    data <- data.frame(
      site =rep(sampleinfo$site, times = sampleinfo$nsamp), #repeats site for number of samples per site
      blank_time=sampleinfo$blanktime, #inserts blank time format 
      botnum = sprintf("%02d",rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp]), #bottle numbers 1-24 repeated for number of samples
      study=sampleinfo$study,samp_type=sampleinfo$samp_type)
    return(data)
    }
  
  field_data <-lapply(1:nrow(sampleinfo),function(x){data_info(sampleinfo[x,])})%>% #using data_info function to make data.frame
    bind_rows()%>%
    mutate(
      randomn=sample(rn:sn,sn, replace=F), #add random numbers
      sample_name = paste0(study, "_", "R", sprintf("%04d",randomn)), 
      field_sample_ID_blank= paste0(study, "_", site, "_", samp_type,"_",blanktime))
  return(field_data )
}#field_label #knowntime=FALSE
}
} #if ISCO





#part two 
if(pt2==TRUE & samp_type =="IS"){
  if (knowntime==TRUE){
    data_sheet<- function(){
      data_arranged <- field_label() %>%
        arrange(sample_name)
      return(data_arranged)
    }
}else if (knowntime==FALSE){
  sampleinfo <- data.frame(nsamp=nsamp,interval=interval,start=start,DOC=DOC,DOM=DOM,CCAL=CCAL,Levoglucosan=Levoglucosan,Excess=Excess) 
  #first part 
  data_sheet <- function(){
    site_label_blank <- function(sampleinfo){
      timestart <- as.POSIXct(sampleinfo$start)
      timelist <- seq(timestart, by=paste0(sampleinfo$interval, " hours"), length = sampleinfo$nsamp)
      return(
        data.frame(time_PST=as.character(timelist),datetime_PST=format(timelist,"%Y%m%d%_%H%M")))} #might need to split this for the new template
    
    updated_labels <- lapply(1:nrow(sampleinfo),function(x){site_label_blank(sampleinfo[x,])}) %>% #makes the frame
      bind_rows()%>%
      mutate(blank_labels)%>%
      mutate(field_sample_ID= paste0(study, "_", site, "_",samp_type,"_",datetime_PST))%>%
      arrange(sample_name)
    return(updated_labels)
  } #unknowntime
  
}

} #if pt 2 and IS 

#lab_label






####For grab samples ##

if (samp_type=="GB"){
site<- read_xlsx(file,"Sample-DateTime")%>%select(site)
sn <- nrow(site) #how many sites there are  
sampleinfo<- data.frame(site=site,study=study,samp_type=samp_type,blanktime=blanktime,sn=sn,rn=rn)
field_label <- function(){ #this will add blank labels in 
  datablank<- data.frame(site)%>%
    mutate(study=study,
           samp_type=samp_type,
           blanktime=blanktime,
           randomn= sample(rn:sn,sn, replace=F),
           date_text_PST = "",
           time_text_PST ="") 
  field_data <- datablank %>%
    mutate(sample_name=paste0( study, "_", "R", sprintf("%04d", randomn)), #adds sample name
           field_sample_ID_blank = paste0(study, "_",
                                          site, "_",
                                          samp_type,"_",
                                          blanktime))#adds blank sample ID
  field_data
  return(  field_data%>%select("field_sample_ID_blank",
                               "sample_name"))
}  



}#if GB 

if(pt2==TRUE & samp_type=="GB"){
  sampleinfo<- data.frame(study=study,samp_type=samp_type,blanktime=blanktime,sn=sn)
  data<- read_xlsx(file,"Sample-DateTime")
  sample_ID <- function(){
    updated_labels <-data %>%
      mutate(
        sample_ID =
          paste0(study, "_",
                 data$site, "_",
                 samp_type,"_",
                 data$date_text_PST, "_",
                 data$time_text_PST))
    return(updated_labels)
    }
    
  }
  
  

if (pt2==TRUE){ 
data_sheet <- function(){  
  data <-sample_ID()%>%
    arrange(sample_name)}
return(data_sheet)
}
  #if pt 2 


if (pt2==TRUE){ 
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
  
  lab_labels <- data_sheet() %>% #putting this here so code doesn't have to get confused <- updated_labels %>% #put in so code isn't confused 
    slice(rep(1:n(), each = analysisn)) %>%
    mutate(datacodes,
           lab_sample_ID = #add in lab labels 
             paste0(sample_ID, "_",
                    datacodes$code ))
  
  return(lab_labels%>%
           select(lab_sample_ID,sample_name,analysis)
         )
}
} #pt 3










#katies function 
create_datasheets <- function(df, saveloc, samptype, type, open=FALSE, overwrite=FALSE){
  stopifnot(type %in% c("datetime", "fieldlabel", "lablabel", "datasheet"),
            samptype %in% c("isco", "grab"))
  
  #locate template/workbook
  if(file.exists(saveloc) & !overwrite){
    #load existing workbook
    path <- saveloc
    
  }else{
    #get template
    path <- list.files(getwd(), pattern = "^wqsample-label-datasheet-templates[.]xlsx", full.names = TRUE, recursive = TRUE)
    if(length(path) == 0){
      #download from github as a backup
      path <- "https://github.com/wildfire-water-security/WWS-standard-methods/raw/refs/heads/main/Templates/water-quality-analysis/wqsample-label-datasheet-templates.xlsx"
    }
  }
  
  #read in template
  tab <- wb_load(path)
  
  #add datetime
  if(type == "datetime"){
   
    #don't overwrite dates/times
    exist_data <- tab %>% wb_to_df(sheet="Sample-DateTime")
    if(any(!is.na(exist_data$`Date (YYYYMMDD)`)) | any(!is.na(exist_data$`Time (HHMM)`))){
    df <- df %>% select(-c(date_text_PST, time_text_PST))
    }
    tab <- tab %>% wb_add_data(sheet= "Sample-DateTime",
                               x=df, dims="A2",col_names = FALSE)
  }
  
  #add labels
  if(type == "lablabel"){
    tab <- tab %>% wb_add_data(sheet= "Labels-Lab", x= df, dims="A2", col_names=FALSE)
  }
  
  if(type == "fieldlabel"){
    tab <- tab %>% wb_add_data(sheet= "Labels-Field", x= df, dims="A2", col_names=FALSE)
  }
  
  #add datasheets
  if(type == "datasheet"){
    ids <- tab %>% wb_to_df(sheet="Sample-DateTime") %>% select('sample_name','sample_ID') %>% arrange(`sample_name`)
    tab <- tab %>% wb_add_data(sheet= "DS-Bottle Numbers", dims="A3", x=ids, col_names=FALSE)
    tab <- tab %>% wb_add_data(sheet= "DS-Filter Weight", dims="A3", x=ids, col_names=FALSE)
    
    if (samp_type == "IS"){
      tab <- tab %>% wb_add_data(sheet= "DS-ISCO Bottle Weight", dims="A3", x=ids, col_names=FALSE)}
  }
  
  #save
  wb_save(tab, file=saveloc)
  
  
  #open if requested
  if(open){
    shell.exec(file.path(getwd(), saveloc))
  }else{
    cat(paste0("information populated and saved to:\n", saveloc))
  }
  
  
  
  
}


























