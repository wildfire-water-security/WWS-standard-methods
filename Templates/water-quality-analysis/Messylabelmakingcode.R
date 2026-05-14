library(dplyr)
library(stringr)
library(kableExtra)
library(lubridate)


######### Filling in codes##########
project <- "WWS" #the large project code
study <- "SSEMG" #the study code
samp_type <- "IS" #(GB,IS, SS, etc.)
rn <- 1 #the stat of a random number sequence
site <- c("STA","COA")
nsamp <- c(12,12) #amount of samples per site
start <- c("2026-02-23 15:00","2026-02-23 17:00")
interval <- c(2,2)

####analysis#####
DOM <- TRUE
DOC <- FALSE
CCAL <- TRUE
Levoglucosan <- TRUE
Excess <- TRUE

##########code########
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
fielddata <- datalabels #just re-naming it for clarity later

#add analytes that are TRUE
analysis_opt <- c("CN", "AQ", "NU", "EX", "LV")
analysis_choice <- c(DOC, DOM, CCAL, Excess, Levoglucosan)

analysiscodes <- analysis_opt[analysis_choice]
# if(DOC) {analysiscodes <- append(analysiscodes,"CN")}
# if(DOM) {analysiscodes <- append(analysiscodes,"AQ")}
# if(CCAL) {analysiscodes <- append(analysiscodes,"NU")}
# if(Excess) {analysiscodes <- append(analysiscodes,"EX")}
# if(Levoglucosan) {analysiscodes <- append(analysiscodes,"LV")}
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
labdata <- datalabels %>% slice(rep(1:n(), each = analysisn)) %>%
mutate(datacodes,
       lab_sample_ID =
        paste0(field_sample_ID, "_", datacodes$code ))

print(labdata)

#printing
print_data_field <- fielddata %>% select(field_sample_ID, sample_name,botnum )
print_data_field #for field
write.csv(print_data_field, "field.csv", row.names=FALSE)

print_data_lab <- labdata %>% select(lab_sample_ID,sample_name,analysis )
print_data_lab #for lab
write.csv(print_data_lab, "lab.csv", row.names=FALSE)
