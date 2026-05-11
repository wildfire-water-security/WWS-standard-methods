library(dplyr)
library(stringr)
library(kableExtra)


######### Filling in codes##########
#codes
project <- "WWS" #the large project code
study <- "SSEMG" #the study code 
samp_type <- "GB" #(GB,IS, SS, etc.)
#the stat of a random number sequence 
rn <- 1
####
#attempt of being a vector?
#is there a way to make this into a vector?
sps <- c(12,12) #amount of samples per site
sites <- c("STA","COA") #site 1, 2
#repeat sites names how many samples- seems to kinda work?
sitelist <- rep(sites, times = sps)
sitelist
#
#Not vector- can't get the code below to work with a vector. Don't like these
sn <- 24 #total sample amount. Is there a way to use the sps to print this???
samples1 <- 12 #samples for site 1
samples2 <-12 #samples for site 2
site1 <- "STA"
site2 <- "COA"
#time starts how to do vector????
timestart1 <- as.POSIXct("2026-02-23 15:00:00") #start time of site 1
timestart2 <- as.POSIXct("2026-02-23 17:00:00") #start time of site 2

#works with one site at a time I guess- how to do for loop?? 
#would like to use the "sps" and "sites" here
date1 <- seq(timestart1, by = "2 hours", length = samples1)
date1
date2 <- seq(timestart2, by = "2 hours", length = samples2)
date2
timecombined <- c(date1,date2)
timecombined
sitetime <- data.frame (time=timecombined,sites=sitelist)


#########analysis##################
DOC <- "TRUE"
DOM <- "FALSE"
CCAL <- "TRUE"
EXCESS <- "TRUE"
Levoglucosan <- "TRUE"
analyte <- c("DOC", "DOM", "CCAL", "LEVO", "EXCESS") #analysis

#how to make it less messy??? this is awful 
analysiscodes <- c()
if(DOC == "TRUE") {analysiscodes <- append(analysiscodes,"CN")}
if(DOM == "TRUE") {analysiscodes <- append(analysiscodes,"AQ")}
if(CCAL == "TRUE") {analysiscodes <- append(analysiscodes,"NU")}
if(EXCESS == "TRUE") {analysiscodes <- append(analysiscodes,"EX")}
if(Levoglucosan == "TRUE") {analysiscodes <- append(analysiscodes,"LV")}
analysiscodes
analysis <- c()
if(DOC == "TRUE") {analysis <- append(analysis,"Shimadzu")}
if(DOM == "TRUE") {analysis <- append(analysis,"Aqualog")}
if(CCAL == "TRUE") {analysis <- append(analysis,"CCAL")}
if(EXCESS == "TRUE") {analysis <- append(analysis,"Excess")}
if(Levoglucosan == "TRUE") {analysis <- append(analysis,"Levoglucosan")}
analysis
#putting it together- breaks if there is a false :/      
lab <- data.frame(analysis=analysis,analysiscodes = analysiscodes)
lab
#how to make it repeat all analysis for each sample?
lablabels <- rep(lab, times = sps)
#need a way to have ex. AQ added to end of ID and Aqualog added in another row




#####Random IDs ##########
#random ID number starting and sample numbers. W
#would like the sn to autofill itself somehow.
randomID <- sample(rn:sn,sn, replace=F)
print(randomID)
#add leading zeros using sprintf
randomn <- sprintf("%04d", randomID)
print(randomn)
#Random IDs
randomlabel <- paste0( study, "_", "R", randomn )
randomlabel


######### If statements for sample type and analysis??? ###########
if(samp_type == "GB") {}
if(samp_type == "IS") {}


#SampleID: project_site_type_date_time
#RandomID: project_r###
#ISCOID: just the number order it goes in the isco?
#paste function works well to join things together 

paste0("A", "B", "C")

paste("A", "B", "C")

paste("A", "B", "C", sep="_")
randomlabel <- paste("R", randomn, sep="_")



#Dates code again
start_time <- as.POSIXct("2023-01-01 02:00:00")
end_time <- as.POSIXct("2023-01-01 23:00:00")
timestamp_sequence <- seq.POSIXt(start_time, end_time, by = "hour")
print(timestamp_sequence)


#the seq function is good for generating a sequence 
seq(as.Date("2023-02-01"), as.Date("2024-01-03"), by="month")



#we can store all this in a data.frame similar to excel 
data <- data.frame(project = project, 
                   study=study, 
                   site=site1,
                   dates = seq(as.POSIXct("2023-02-01 12:00",), 
                               as.POSIXct("2024-02-02 12:00"), by="2 hours"))

data





#can't get the pipe function to work yet 
print_data <- data %>% select(project, study)
print_data
#write.csv(print_data, "filename.csv", row.names=FALSE)
