library(dplyr)
library(stringr)
library(kableExtra)
library(lubridate)
 
#this one is assuming you are going to the field with blank times and will
#fill the dates into the made spreadsheet and use that to print lab labels



######### Filling in codes##########
project <- "WWS" #the large project code
study <- "BDRK" #the study code
samp_type <- "GB" #(GB,IS, SS, etc.)
rn <- 1 #the stat of a random number sequence
site <- read.csv("example_sites.csv") #import site list. has the header "site"
site
sn <- 20 #how many sites there are
#each label sheet prints 30 samples

##########code########
time <- "____________" #either have it blank or with sampling date "____________" or "YYYYMMDD____"
time
#data 
datablank<- data.frame(site) %>%
  mutate(study=study,
         samp_type=samp_type,
         blank_time=time,
         randomn= sample(rn:sn,sn, replace=F)
         )

datablank

#add sample names and IDs and blank columns for entering date and time
datalabels<-datablank %>%
  mutate(sample_name=paste0( study, "_", "R", sprintf("%04d", randomn)),
         field_sample_ID_blank = paste0(study, "_",
                           site, "_",
                           samp_type,"_",
                           blank_time),
         date_text = rep(NA, nrow(datalabels)), time_text = rep(NA, nrow(datalabels)))
datalabels
write.csv(datalabels, "alldatablank.csv", row.names=FALSE)

#print blank labels 
fielddata <- datalabels #just re-naming it for clarity later
print_data_field_blank <- fielddata %>% select(field_sample_ID_blank, sample_name,site)
print_data_field_blank #for field
write.csv(print_data_field_blank, "blankfield.csv", row.names=FALSE)

#add times and dates as numbers only 

############next da########################
#not sure how to print with an empty column for date and time*******
datafilled<-read.csv("alldatablank.csv")
datafilled 

#make datetime and sampleID with times

datafilled<-datafilled %>%
  mutate(datetime=
           paste0(sprintf("%06d",datafilled$date_text),
                  sprintf("%04d",datafilled$time_text)),
         field_sample_ID =
           paste0(study, "_",
                  site, "_",
                  samp_type,"_",
                  datetime))

datafilled 
#WHY WON'T YOU WORK*************
#datafilled2 <-mutate(datafilled,
                     #datetime=
                       #as.POSIXct(paste(date,time),format = "%Y%m%d%H%M"))


#texttime <- as.character(paste(datafilled$date,datafilled$time))
#texttime1 <- as.POSIXct(texttime)
##############pretend that worked

####analysis#####
DOM <- TRUE
DOC <- FALSE
CCAL <- TRUE
Levoglucosan <- TRUE
Excess <- TRUE


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
labdata <- datafilled %>% slice(rep(1:n(), each = analysisn)) %>%
mutate(datacodes,
       lab_sample_ID =
        paste0(study, "_",
                      site, "_",
                      samp_type,"_",
                      datetime,"_", 
                      datacodes$code ))

labdata
#printing

print_data_lab <- labdata %>% select(lab_sample_ID,sample_name,analysis)
print_data_lab #for lab
write.csv(print_data_lab, "lab.csv", row.names=FALSE)
