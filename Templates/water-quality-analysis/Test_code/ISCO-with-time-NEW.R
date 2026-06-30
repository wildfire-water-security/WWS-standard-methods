#load libraries
library(dplyr)
library(stringr)
library(lubridate)
library(writexl)
library(readxl)
############ Fill out##############

#sample information:
study <- "SSEMG" #the study code
samp_type <- "IS" #(GB,IS, SS, etc.)
rn <- 1 #the stat of a random number sequence
site <- c("STA","COA")
nsamp <- c(24,24) #amount of samples per site. Matches the site vector above.
interval <- c(2,2) #interval of sampling ISCO is programmed for, usually 2hrs
sc <- "20260224" #what sample campaign this is, the date it starts "YYYYMMDD". This is for file name purposes
TZone <- "Etc/GMT+8"

#optional: add grab samples. Not sure how to code this yet.How to include in random numbers???? Make own dataframe? Just have blank IDs in the lab?

GRAB <- FALSE #fill TRUE for having grab samples, FALSE for not. 
grabsampn <- 4 #total number of grab samples expected
grabsamp <- c(2,2) #samples per site 
blk <- TRUE #add procedural blank?
blkname <-"mr_blanky"
#sample times:
#fill if you have known time:
start <- c("2026-02-23 15:00","2026-02-23 17:00")

#analysis, fill TRUE if using them, FALSE if not:
DOM <- TRUE
DOC <- FALSE
CCAL <- TRUE
Levoglucosan <- TRUE
Excess <- TRUE

#optional formatting- could delete this and just have it in the code
#what format the files will end with. ex this returns: study_sapmpletype_date_label-table 
fileformat <- paste0(study,"_",samp_type,"_",sc,"_label-table.xlsx") 
fileformat


########## Code pt.1- don't change #######################

#If doing blank samples, rerun this section in part two
#start#
sn <- sum(nsamp) #total sample amount from the vector that was filled out
#adding sample information
sampleinfo <- data.frame(site=site,nsamp=nsamp,start=start,interval=interval,study=study,samp_type=samp_type) #all sample info for the function
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


#function to fill out site, study, sample type, time stamps, and bottle numbers
#If statement, defines the site_label function depending on if there is a known time or not

site_label <- function(sampleinfo){
    sitelist <- rep(sampleinfo$site, times = sampleinfo$nsamp) #if you turn into a df to return you wouldn't need this, it would automatically recycle
    timestart <- as.POSIXct(sampleinfo$start, tz = TZone)
    timelist <- seq(timestart, by=paste0(sampleinfo$interval, " hours"), length = sampleinfo$nsamp)
    botnum <- rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp] #added interval so if we ever want to change it, it's flexible
    return(data.frame(site =sitelist, time_PST =timelist, botnum = botnum, study=sampleinfo$study,samp_type=sampleinfo$samp_type,datetime_PST=format(timelist,"%Y%m%d%H%M"))) #could include other label things in here too
  }

#add grabs
if(GRAB==TRUE){
  
}


#make the sample table from the function and add random numbers to it. Re-running will change the random numbers.
data <- lapply(1:nrow(sampleinfo),function(x){site_label(sampleinfo[x,])})%>%
  bind_rows()
  

  
  data
grabby <- data.frame(grabsamp=site, re)
grabby
  grabs <- append(data,)

#old
# data <- lapply(1:nrow(sampleinfo),function(x){site_label(sampleinfo[x,])})%>%
#   bind_rows() %>%
#   mutate (randomn= sample(rn:sn,sn, replace=F))
# data


#add sample names and IDs
  datalabels<-data %>%
    mutate(sample_name=paste0( study, "_", "R", sprintf("%04d", randomn)),
           field_sample_ID = paste0(study, "_",
                                    site, "_",
                                    samp_type,"_",
                                    datetime_PST))

datalabels


#######Adding analytes- code pt.2########


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
datasortlab <- datalabels[order(datalabels$randomn),]
datasortlab #use this to print datasheets out

#repeat sample rows for each analysis
labdata <- datasortlab %>% slice(rep(1:n(), each = analysisn)) %>%
  mutate(datacodes,
         lab_sample_ID =
           paste0(field_sample_ID, "_", datacodes$code ))
print(labdata)



#print csv for lab sample labels
# safe.write.csv(labdata %>% select(lab_sample_ID,sample_name,analysis), 
#           paste0(lab_file,fileformat), row.names=FALSE) 

# Named list of data.frames for writing to an excel file


#
write_data.frames <- list(
  field_labels = datalabels %>% select(sample_name,field_sample_ID,botnum), #no sample ID?
  filtering = datasortlab %>% select(sample_name,field_sample_ID),
  lab_labels   = labdata %>% select(sample_name,lab_sample_ID,analysis),
  label_metadata = datalabels %>% select(site,datetime_PST,study,samp_type,randomn,botnum) #jsut extra didn't know what to add
)

safe.write_xlsx(write_data.frames,
                path = fileformat)
