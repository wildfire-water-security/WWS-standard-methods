# Functions to safely write slsxs and csvs without overwriting anything
safe.write_xlsx <- function(x, path, ...) {
  if(file.exists(path)) {
    stop("The file you are trying to write already exists. DO NOT OVERWRITE UNLESS YOU MEAN TO!\n",
         "We don't want to overwrite randomized sample names\n",
         "Check ", path, " to see if the already written file is ok.")
  }
  write_xlsx(x, path, ...)
}


#' Title
#' 
#' Description here
#'
#' @param x 
#' @param path 
#' @param ... 
#'
#' @returns
#'
safe.write_csv <- function(x, path, ...) {
  if(file.exists(path)) {
    stop("The file you are trying to write already exists. DO NOT OVERWRITE UNLESS YOU MEAN TO!\n",
         "We don't want to overwrite randomized sample names\n",
         "Check ", path, " to see if the already written file is ok.")
  }
  write.csv(x, path, ...)
} #had to add this because couldn't get blank files to work wth xlsx



#functions for data labels 

#function 1  
data_label <- function(){
  sampleinfo <- data.frame(site=site,nsamp=nsamp,interval=interval,study=study,samp_type=samp_type,blanktime=blanktime, sn=sn, rn=rn)
  data_info <- function(sampleinfo){ #this function creates 
    return(data.frame(
      site =rep(sampleinfo$site, times = sampleinfo$nsamp),
      blank_time=sampleinfo$blanktime,
      botnum = sprintf("%02d",rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp]), study=sampleinfo$study,samp_type=sampleinfo$samp_type
    ))}
  x<-lapply(1:nrow(sampleinfo),function(x){data_info(sampleinfo[x,])})%>%
    bind_rows()%>%
    mutate(
      randomn=sample(rn:sn,sn, replace=F),
      sample_name = paste0(study, "_", "R", sprintf("%04d",randomn)),
      field_sample_ID_blank= paste0(study, "_", site, "_", samp_type,"_",blanktime))
}

#second version

lab_label <- function(){
  sampletimeinfo <- data.frame(nsamp=nsamp,interval=interval,start=start,DOC=DOC,DOM=DOM,CCAL=CCAL,Levoglucosan=Levoglucosan,Excess=Excess)   
  blank_labels <- data.frame(read.csv(blank_fileformat))
  site_label_blank <- function(sampletimeinfo){
    timestart <- as.POSIXct(sampletimeinfo$start)
    timelist <- seq(timestart, by=paste0(sampletimeinfo$interval, "hours"), length = sampletimeinfo$nsamp)
    return(
      data.frame(time_PST=as.character(timelist),datetime_PST=format(timelist,"%Y%m%d%_%H%M")))}

  analysis_opt <- c("CN", "AQ", "NU", "LV", "EX")
  analysis_choice <- c(DOC, DOM, CCAL, Levoglucosan, Excess)
  analysiscodes <- analysis_opt[analysis_choice]
  analysisn <- length(analysiscodes) #how many analytes there are
  datacodes<- data.frame(code = rep(analysiscodes,times=sn)) %>% #add labcodes
    mutate(analysis = case_when(
      code == "CN" ~ "Shimadzu",
      code == "AQ" ~ "Aqualog",
      code == "NU" ~ "CCAL",
      code == "LV" ~ "Levoglucosan",
      code == "EX" ~ "Excess",
      TRUE ~ NA_character_))
  
  data_labels <- lapply(1:nrow(sampletimeinfo),function(x){site_label_blank(sampletimeinfo[x,])}) %>% #makes the frame
    bind_rows()%>%
    mutate(blank_labels)%>%
    mutate(field_sample_ID= paste0(study, "_", site, "_",samp_type,"_",datetime_PST))%>%
    arrange(randomn)%>%
    slice(rep(1:n(), each = analysisn))%>%
    mutate(datacodes)
  
  return(data_labels)
  
}


