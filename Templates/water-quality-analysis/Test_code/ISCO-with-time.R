#load libraries
library(dplyr)
library(stringr)
library(kableExtra)
library(lubridate)
library(writexl)


######### Fill out##########
project <- "WWS" #the large project code
study <- "SSEMG" #the study code
samp_type <- "IS" #(GB,IS, SS, etc.)
rn <- 1 #the stat of a random number sequence
site <- c("STA","COA")
nsamp <- c(24,24) #amount of samples per site
start <- c("2026-02-23 15:00","2026-02-23 17:00")
interval <- c(2,2) #2 hrs
sc <- "20260224" #what sample campaign this is
#how to print off just the date???????
#timestart<- as.charcter(start,[1)
#timestart
####analysis- fill out#####
DOM <- TRUE
DOC <- FALSE
CCAL <- TRUE
Levoglucosan <- TRUE
Excess <- TRUE

####optional formatting- fill out#####
alldata_file <- "alldata" #what to call header for alldata file
field_file <-  "fieldsamples" #what to call header for field label file
lab_file <-  "labsamples" #what to call header for lab label file
#what format the files will end with. ex this returns:
#header_study_sapmpletype_date 
fileformat <- paste0(study,"_",samp_type,"_",sc,"_label-table.xlsx") 
fileformat



##########Code#######################
# Function to safely write csvs without overwriting anything
safe.write_xlsx <- function(x, path, ...) {
  if(file.exists(path)) {
    stop("The file you are trying to write already exists. DO NOT OVERWRITE UNLESS YOU MEAN TO!\n",
         "We don't want to overwrite randomized sample names\n",
         "Check ", path, " to see if the already written file is ok.")
  }
  write_xlsx(x, path, ...)
}

sampleinfo <- data.frame(site=site,nsamp=nsamp,start=start,interval=interval,study=study,samp_type=samp_type) #all sample info for function
sn <- sum(nsamp) #total sample amount

site_label <- function(sampleinfo){
  sitelist <- rep(sampleinfo$site, times = sampleinfo$nsamp) #if you turn into a df to return you wouldn't need this, it would automatically recycle
  timestart <- as.POSIXct(sampleinfo$start)
  timelist <- seq(timestart, by=paste0(sampleinfo$interval, " hours"), length = sampleinfo$nsamp)
  botnum <- rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp] #added interval so if we ever want to change it, it's flexible

  return(data.frame(site =sitelist, time=timelist, botnum = botnum, study=sampleinfo$study,samp_type=sampleinfo$samp_type)) #could include other label things in here too
}

#make data table and add random numbers to it. Re-running will change the random numbers
data <- lapply(1:nrow(sampleinfo),function(x){site_label(sampleinfo[x,])})%>%
  bind_rows() %>%
  mutate (randomn= sample(rn:sn,sn, replace=F))

data


#add sample names and IDs
datalabels<-data %>%
  mutate(sample_name=paste0( study, "_", "R", sprintf("%04d", randomn)),
         field_sample_ID = paste0(study, "_",
                           site, "_",
                           samp_type,"_",
                           format(time,"%Y%m%d%H%M")))
datalabels





#print all the data just in case
# safe.write.csv(datalabels, paste0(alldata_file,fileformat), row.names=FALSE) 
#print csv for field sample labels
# safe.write.csv(datalabels %>% select(field_sample_ID, sample_name,botnum),
#           paste0(field_file,fileformat), row.names=FALSE) 


#add analytes that are TRUE
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

#sort the data by the randomn, since they should be filtered that way
datasortlab <- datalabels[order(data$randomn),]
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

write_data.frames <- list(
  field_labels = datalabels %>% select(field_sample_ID,sample_name,botnum),
  lab_labels   = labdata %>% select(lab_sample_ID,sample_name,analysis),
  label_metadata = datalabels %>% select(site,time,study,samp_type,randomn,botnum)
)
#THE TIME ISN'T PRINTING CORRECTLY ON HERE

safe.write_xlsx(write_data.frames,
           path = fileformat)
