#load libraries
library(dplyr)
library(stringr)
library(lubridate)
library(writexl)

############ Fill out##############

#sample information:
project <- "WWS" #the large project code
study <- "SSEMG" #the study code
samp_type <- "IS" #(GB,IS, SS, etc.)
rn <- 1 #the stat of a random number sequence
site <- c("STA","COA")
nsamp <- c(24,24) #amount of samples per site. Matches the site vector above.
interval <- c(2,2) #interval of sampling ISCO is programmed for, usually 2hrs
sc <- "20260220" #what sample campaign this is, the date it starts "YYYYMMDD". This is for file name purposes. If it is an unknown then do something like "study_campaign_#"
TZone <- "Etc/GMT+8"


#sample times:
knowntime <- TRUE #fill TRUE When you have a time start and FALSE for if you don't know the time start. Knowing the time is the most efficient method.Important for If statements.  
#fill if you have known time:
start <- c("2026-02-23 15:00","2026-02-23 17:00")
#fill if you don't know time:
blanktime <- "____________" #either have it blank  "____________" or with sampling date "YYYYMMDD____" or "YYYYMM______", etc

#analysis, fill TRUE if using them, FALSE if not:
DOM <- TRUE
DOC <- FALSE
CCAL <- TRUE
Levoglucosan <- TRUE
Excess <- TRUE



#optional change:
#what format the files will end with. ex this returns: study_sapmpletype_date_label-table 
fileformat <- paste0(study,"_",samp_type,"_",sc,"_label-table.xlsx") #the file name
fileformat
blank_fileformat <- paste0("blank_",study,"_",samp_type,"_",sc,"_label-table.csv") 
blank_fileformat 


########## Code pt.1- don't change #######################

#If doing blank samples, rerun this section in part two
#start#
sn <- sum(nsamp) #total sample amount from the vector that was filled out



#adding sample information
sampleinfo <- data.frame(site=site,nsamp=nsamp,start=start,interval=interval,study=study,samp_type=samp_type,blanktime=blanktime) #all sample info for the function
#would need to alter if it includes grabs 




# Functions to safely write slsxs and csvs without overwriting anything
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
} #had to add this because couldn't get blank files to work wth xlsx
#end#


#function to fill out site, study, sample type, time stamps, and bottle numbers
#If statement, defines the site_label function depending on if there is a known time or not
if (knowntime == TRUE){
  #Known time 
  site_label <- function(sampleinfo){
    sitelist <- rep(sampleinfo$site, times = sampleinfo$nsamp) #if you turn into a df to return you wouldn't need this, it would automatically recycle
    timestart <- as.POSIXct(sampleinfo$start)
    timelist <- seq(timestart, by=paste0(sampleinfo$interval, " hours"), length = sampleinfo$nsamp)
    botnum <- rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp] #added interval so if we ever want to change it, it's flexible
    return(data.frame(site =sitelist, timelist=timelist, 
                      botnum = botnum, study=sampleinfo$study,
                      samp_type=sampleinfo$samp_type,
                      datetime_PST=format(timelist,"%Y%m%d%H%M"))) #could include other label things in here too
  }
} else if (knowntime == FALSE){
  site_label <- function(sampleinfo){
    sitelist <- rep(sampleinfo$site, times = sampleinfo$nsamp) #if you turn into a df to return you wouldn't need this, it would automatically recycle
    botnum <- rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp] #added interval so if we ever want to change it, it's flexible
    return(data.frame(site =sitelist,blank_time=sampleinfo$blanktime, botnum = botnum, study=sampleinfo$study,samp_type=sampleinfo$samp_type)) #could include other label things in here too
  }
}
#make the sample table from the function and add random numbers to it. Re-running will change the random numbers.
data <- lapply(1:nrow(sampleinfo),function(x){site_label(sampleinfo[x,])})%>%
  bind_rows()%>%
  mutate (randomn= sample(rn:sn,sn, replace=F))
data


####
#add sample names and IDs
if (knowntime == TRUE){
  datalabels<-data %>%
    mutate(sample_name=paste0( study, "_", "R", sprintf("%04d", randomn)),
           field_sample_ID = paste0(study, "_",
                                    site, "_",
                                    samp_type,"_",
                                    datetime_PST))
}else if (knowntime == FALSE ){
  datalabels<-data %>%
    mutate(sample_name=paste0( study, "_", "R", sprintf("%04d", randomn)),
           field_sample_ID_blank = paste0(study, "_",
                                          site, "_",
                                          samp_type,"_",
                                          blanktime))}
datalabels

#this prints your data labels out if you have blank samples, since you will need to rerun this code again. 
if(knowntime == FALSE){
  safe.write_csv(datalabels,
                 path = blank_fileformat,
                 row.names = FALSE)
  stop("DO NOT CONTINUE UNTIL YOU HAVE YOUR TIMES. SAVE YOUR CODE")}



#######Adding analytes- code pt.2########

#section for blank labels 
#if you have blank times: if the code isn't still loaded: rerun the 'fillout' section and the top section in 'code pt.1'. leave knownsamples = FALSE
#and import data set 

newtimes <- c("2026-02-23 15:00","2026-02-23 17:00")



if(knowntime == FALSE){
  datalabels <- read.csv(blank_fileformat)
}

datalabels



if(knowntime == FALSE){
  site_label_blank <- function(sampleinfo){
    # Time code - need to fix time zones
    timestart <- as.POSIXct(sampleinfo$start, tz = TZone)
    timelist <- seq(timestart, by=paste0(sampleinfo$interval, " hours"), length = sampleinfo$nsamp)
    return(data.frame(time_PST=as.character(timelist),datetime_PST=format(timelist,"%Y%m%d%H%M"))) #could include other label things in here too
  }} 
data_times <- lapply(1:nrow(sampleinfo),function(x){site_label_blank(sampleinfo[x,])}) %>%
  bind_rows %>%
  mutate(datalabels)

data_times <- mutate(data_times, field_sample_ID = paste0(study, "_",
                                                          site, "_",
                                                          samp_type,"_",
                                                          datetime_PST))


data_times
datalabels 
data
#print all the data just in case
# safe.write.csv(datalabels, paste0(alldata_file,fileformat), row.names=FALSE) 
#print csv for field sample labels
# safe.write.csv(datalabels %>% select(field_sample_ID, sample_name,botnum),
#           paste0(field_file,fileformat), row.names=FALSE) 



#######Adding analytes- code pt.2########

#adding in processing information:
#add analytes that are TRUE
analysis_opt <- c("CN", "AQ", "NU", "EX", "LV")
analysis_choice <- c(DOC, DOM, CCAL, Excess, Levoglucosan)
analysiscodes <- analysis_opt[analysis_choice]
analysiscodes
analysisn <- length(analysiscodes) #how many analytes there are
samples <- analysisn*sn #how many total filtered samples there will be



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

#sort the data by the randomn, since they should be filtered that way
datasortlab <- data_times[order(data_times$randomn),]
datasortlab #use this to print datasheets out

#repeat sample rows for each analysis
labdata <- datasortlab %>% slice(rep(1:n(), each = analysisn)) %>%
  mutate(datacodes,
         lab_sample_ID =
           paste0(field_sample_ID, "_", datacodes$code ))
print(labdata)




#
write_data.frames <- list(
  field_labels = data_times %>% select(field_sample_ID,sample_name,botnum),
  lab_labels   = labdata %>% select(lab_sample_ID,sample_name,analysis),
  label_metadata = data_times %>% select(site,time_PST,study,samp_type,randomn,botnum)
)
#THE TIME ISN'T PRINTING CORRECTLY ON HERE

safe.write_xlsx(write_data.frames,
                path = fileformat)
