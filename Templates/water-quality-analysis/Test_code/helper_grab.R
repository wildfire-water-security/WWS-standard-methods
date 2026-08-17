library(dplyr)
library(stringr)
library(lubridate)
library(writexl)



#safewrite function 
safe.write_xlsx <- function(x, path, ...) {
  if(file.exists(path)) {
    stop("The file you are trying to write already exists. DO NOT OVERWRITE UNLESS YOU MEAN TO!\n",
         "We don't want to overwrite randomized sample names\n",
         "Check ", path, " to see if the already written file is ok.")
  }
  write_xlsx(x, path, ...)
}
safe.write_csv <- function(x, path, ...) {
  if(file.exists(path)) {
    stop("The file you are trying to write already exists. DO NOT OVERWRITE UNLESS YOU MEAN TO!\n",
         "We don't want to overwrite randomized sample names\n",
         "Check ", path, " to see if the already written file is ok.")
  }
  write.csv(x, path, ...)
} 



sn <- nrow(site) #how many sites there are
sampleinfo<- data.frame(site=site,study=study,samp_type=samp_type,blanktime=blanktime,sn=sn,rn=rn)


labels <- function(){
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
  print <- blank_labels %>% select(
    sample_name,field_sample_ID_blank,randomn,study,samp_type,site,date_text_PST,time_text_PST) #reorder so it's easier for data entry
return(safe.write_csv(print, "example-grab-blank.csv", row.names=FALSE))
}


labdata<- function(){
  data<-datalabels %>%
  mutate(datetime_PST=
           paste0(sprintf("%06d",datalabels$date_text_PST),"_",
                  sprintf("%04d",datalabels$time_text_PST)),
         field_sample_ID =
           paste0(study, "_",
                  site, "_",
                  samp_type,"_",
                  datetime_PST))
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
datasortlab<-data[order(data$randomn),]

labdata <- datasortlab %>%
  slice(rep(1:n(), each = analysisn)) %>%
  mutate(datacodes,
         lab_sample_ID =
           paste0(study, "_",
                  site, "_",
                  samp_type,"_",
                  datetime_PST,"_", 
                  datacodes$code ))
write_data.frames <- list(
  lab_labels   = labdata %>% select(sample_name,lab_sample_ID,analysis),
  lab_filter_sheets = datasortlab %>% select(sample_name,field_sample_ID)
)
return(safe.write_xlsx(write_data.frames,
                path = "example-grab-finished.xlsx"))

}







