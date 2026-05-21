library(dplyr)
library(stringr)
library(kableExtra)
library(lubridate)
 
#this one is assuming you are going to the field with blank times and will
#fill the dates into the made spreadsheet and use that to print lab labels
#this is for unique site names!
# Function to safely write csvs without overwriting anything
safe.write.csv <- function(x, file, ...) {
  if(file.exists(file)) {
    stop("The file you are trying to write already exists. DO NOT OVERWRITE UNLESS YOU MEAN TO!\n",
         "We don't want to overwrite randomized sample names\n",
         "Check ", file, " to see if the already written file is ok.")
  }
  write.csv(x, file, ...)
}

######### Filling in codes##########
project <- "WWS" #the large project code
study <- "BDRK" #the study code
samp_type <- "GB" #(GB,IS, SS, etc.)
rn <- 1 #the stat of a random number sequence
site <- read.csv("example_sites.csv") #import your site list. Has the header "site"
#not sure how to do this without a specific header...
site
time <- "____________" #either have it blank or with sampling date "____________" or "YYYYMMDD____"
time
sc <- "20260823" #what sample campaign this is. Either date ex."20260223 or campaign ex. "1"


####optional formatting- fill out#####
alldata_blank_file <- "alldata_blank" #what to call header for blank times alldata file
alldata_file <- "alldata" #what to call header for alldata file you fill out
field_file <-  "fieldsamples" #what to call header for field label file
lab_file <-  "labsamples" #what to call header for lab label file

#what format the files will end with. ex this returns:
#header_study_sapmpletype_date 
fileformat <- paste0("_",study,"_",samp_type,"_",sc,".csv") 
fileformat

#save your entries (for option #4 in next step)
filenames <- data.frame(project=project,study=study,sc=sc,
                      alldata_file=alldata_file,field_file=field_file,lab_file=lab_file)
filenames
safe.write.csv(filenames,paste0("filenames","_",fileformat))


##########code########
sn <- nrow(site) #how many sites there are
sn #each label sheet prints 30 samples

#Adds random numbers to data frame
datablank<- data.frame(site) %>%
  mutate(study=study,
         samp_type=samp_type,
         blank_time=time,
         randomn= sample(rn:sn,sn, replace=F)
         ) 
datablank
#adds blank sample IDs, sample names, and adds text columns to fill out for date and time
blanklabels <- datablank %>%
  mutate(sample_name=paste0( study, "_", "R", sprintf("%04d", randomn)), #adds sample name
         field_sample_ID_blank = paste0(study, "_",
                                        site, "_",
                                        samp_type,"_",
                                        blank_time), #adds blank sample ID
         date_text = rep('', nrow(datablank)), time_text = rep('', nrow(datablank))) #adds column names to fill out. Puts NAs in them though
blanklabels
#print blank labels 
safe.write.csv(blanklabels, paste0("alldata_blank",fileformat), row.names=FALSE)
#Before filtering- fill out blank columns for date and time and print lab labels
#Make sure to do it in text only! 
  #For example: for 8:00am write 800 instead of 08:00, for 2pm write 1400

#print blank field labels
safe.write.csv(blanklabels %>% select(field_sample_ID_blank, sample_name, site),
          paste0("fieldsamples_blank",fileformat), row.names=FALSE)  #file format is fieldsamples_study_sapmpletype_date




############morning of filtering- should be it's own code########################

#fill out analysis
DOM <- TRUE
DOC <- TRUE
CCAL <- FALSE
Levoglucosan <- TRUE
Excess <- FALSE

DOM

##We don't want to accidentally re-run the top code because it will return new random numbers and overwrite files#

#choose option: (could make this an if statement but think one of these will get deleted)
#option 1: import blank csv, make copy and fill it out
#option 2:import the csv you've already filled out
#option 3: refill out all parts 
#option 4: import data for filenames and create them there

#Reminder:fill out blank columns for date and time and print lab labels
  #Make sure to do it in text only! 
  #For example: for 8:00am write 800 instead of 08:00, for 2pm write 1400

  ##option 1#####
#import your blank csv`
mycsv<- "alldata_BDRK_GB_20260823.csv"
header <- "alldata_blank_" #what the header was on the file (before site)
#make a copy
safe.write.csv(mycsvblank,sub(header,'alldata_filled_',mycsvblank))
#fill out the date_text and time_text sections of new csv and save 
mycsv <- sub(header,'alldata_filled_',mycsvblank)
mycsv

####optional formatting- fill out#####
###alldata_filled <- "alldata_filled" #what to call header for filled out alldata file
field_file <-  "fieldsamples_blank" #what to call header for field label file
lab_file <-  "labsamples" #what to call header for lab label file

  ##option 2#####
#Note: the date_text and time_text columns were already filled out.
#Note: rename the file to take the _blank part out
#import csv
mycsv<- "alldata_BDRK_GB_20260823.csv"
header <- "alldata_" #what the header was on the file (before site)

####optional formatting- fill out#####
###alldata_filled <- "alldata_filled" #what to call header for filled out alldata file
field_file <-  "fieldsamples_blank" #what to call header for field label file
lab_file <-  "labsamples" #what to call header for lab label file



######

#refill out this section: should match the top of the code. We don't want to re-run
#the top code because it will return new random numbers and overwrite files 

project <- "WWS" #the large project code
study <- "BDRK" #the study code
samp_type <- "GB" #(GB,IS, SS, etc.)
sc <- "20260823" #what sample campaign this is. Either date ex."20260223 or campaign ex. "1"

###optional formatting- fill out
fileformat <- paste0("_",study,"_",samp_type,"_",sc,".csv") 
fileformat




labsamplefile<-sub(header,'labsamples_',mycsv)
labsamplefile
newname








#########code#########
#what format the files will end with. ex this returns:#alldata_study_sapmpletype_date 
datafilled<-read.csv(mycsv) 
datafilled 
#make datetimes and the sampleID with times
sn <- nrow(datafilled ) #how many sites there are
sn 
alldata<-datafilled %>%
  mutate(datetime=
           paste0(sprintf("%06d",datafilled$date_text),
                  sprintf("%04d",datafilled$time_text)),
         field_sample_ID =
           paste0(study, "_",
                  site, "_",
                  samp_type,"_",
                  datetime))

alldata
#WHY WON'T YOU WORK*************
#datafilled2 <-mutate(datafilled,
                     #datetime=
                       #as.POSIXct(paste(date,time),format = "%Y%m%d%H%M"))


#texttime <- as.character(paste(datafilled$date,datafilled$time))
#texttime1 <- as.POSIXct(texttime)
##############pretend that worked




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
#repeat samples for each analysis
labdata <- alldata %>% slice(rep(1:n(), each = analysisn)) %>%
mutate(datacodes,
       lab_sample_ID =
        paste0(study, "_",
                      site, "_",
                      samp_type,"_",
                      datetime,"_", 
                      datacodes$code ))

labdata
#printing
safe.write.csv(labdata %>% select(lab_sample_ID,sample_name,analysis), 
          paste0("labsamples",fileformat), row.names=FALSE)

##print_data_lab <- labdata %>% select(lab_sample_ID,sample_name,analysis)
#print_data_lab #for lab
#paste0("fieldsamples",fileformat)
#safe.write.csv(print_data_lab, "lab.csv", row.names=FALSE)
