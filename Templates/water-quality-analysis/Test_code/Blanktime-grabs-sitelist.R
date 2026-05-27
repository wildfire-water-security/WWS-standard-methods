library(dplyr)
library(stringr)
library(lubridate)
library(writexl)

#this one is assuming you are going to the field with blank times and will
#fill the dates into the made spreadsheet and use that to print lab labels
#this is for unique site names!
# Function to safely write csvs without overwriting anything

######### Filling in codes##########
study <- "BDRK" #the study code
samp_type <- "GB" #(GB,IS, SS, etc.)
rn <- 1 #the stat of a random number sequence
site <- read.csv("example_sites.csv") #import your site list. Has the header "site"
#not sure how to do this without a specific header...
site
blank_time <- "____________" #either have it blank or with sampling date "____________" or "YYYYMMDD____"
blank_time
sc <- "20260823" #what sample campaign this is. Either date ex."20260223 or campaign ex. "1"

#fill out analysis
DOM <- TRUE
DOC <- TRUE
CCAL <- FALSE
Levoglucosan <- TRUE
Excess <- FALSE

####optional formatting- fill out#####
alldata_blank_file <- "alldata_blank" #what to call header for blank times alldata file
alldata_file <- "alldata" #what to call header for alldata file you fill out
field_file <-  "fieldsamples" #what to call header for field label file
lab_file <-  "labsamples" #what to call header for lab label file

#what format the files will end with. ex this returns:
#header_study_sapmpletype_date 
fileformat <- paste0("_",study,"_",samp_type,"_",sc,".csv") 
fileformat
#filenames
safe.write.csv(filenames,paste0("filenames","_",fileformat))


fileformat <- paste0(study,"_",samp_type,"_",sc,"_label-table.xlsx") 
fileformat
blank_fileformat <- paste0("blank_",study,"_",samp_type,"_",sc,"_label-table.csv") 
blank_fileformat





##########code########



sn <- nrow(site) #how many sites there are
sn #each label sheet prints 30 samples

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

#Adds random numbers to data frame


datablank<- data.frame(site) %>%
  mutate(study=study,
         samp_type=samp_type,
         blank_time=blank_time,
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
safe.write_csv(blanklabels, paste0(blank_fileformat), row.names=FALSE)

stop("Before printing lab labels out- fill out blank columns for date and time and print lab labels (as text only! example: for 8:00am write 800 instead of 08:00, for 2pm write 1400), save file.")


################
#re-run top of code, library through safewrite functions  ~line 75
#import data 
datalabels <- read.csv(blank_fileformat) #bring in the saved data labels. 
datalabels
stop("check that date_text and time_text rows were filled out in the consol")


#########code#########
#what format the files will end with. ex this returns:#alldata_study_sapmpletype_date 

#make datetimes and the sampleID with times

data<-datalabels %>%
  mutate(datetime=
           paste0(sprintf("%06d",datalabels$date_text),
                  sprintf("%04d",datalabels$time_text)),
         field_sample_ID =
           paste0(study, "_",
                  site, "_",
                  samp_type,"_",
                  datetime))

data
#sort the data by the randomn, since they should be filtered that way
datasortlab <- data[order(data$randomn),]
datasortlab #use this to print datasheets out

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
labdata <- datasortlab %>% slice(rep(1:n(), each = analysisn)) %>%
mutate(datacodes,
       lab_sample_ID =
        paste0(study, "_",
                      site, "_",
                      samp_type,"_",
                      datetime,"_", 
                      datacodes$code ))

labdata
#printing
write_data.frames <- list(
  filter_labels = datasortlab %>% select(field_sample_ID,sample_name),
  lab_labels   = labdata %>% select(lab_sample_ID,sample_name,analysis),
  label_metadata = data %>% select(site,datetime,study,samp_type,randomn) #was just a dump for things
)
safe.write_xlsx(write_data.frames,
                path = fileformat)

