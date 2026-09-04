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





#FOR ISCO UNKNOWN

if (knowntime == TRUE){
  #Known time 
  site_label <- function(sampleinfo){
    sitelist <- rep(sampleinfo$site, times = sampleinfo$nsamp) #if you turn into a df to return you wouldn't need this, it would automatically recycle
    timestart <- as.POSIXct(sampleinfo$start)
    timelist <- seq(timestart, by=paste0(sampleinfo$interval, " hours"), length = sampleinfo$nsamp)
    botnum <- rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp] #added interval so if we ever want to change it, it's flexible
    return(data.frame(site =sitelist, timelist=timelist, 
                      botnum = sprintf("%02d",botnum), 
                      study=sampleinfo$study,
                      samp_type=sampleinfo$samp_type,
                      datetime_PST=format(timelist,"%Y%m%d%_%H%M"))) #could include other label things in here too
  }
} else if (knowntime == FALSE){
  site_label <- function(sampleinfo){
    sitelist <- rep(sampleinfo$site, times = sampleinfo$nsamp) #if you turn into a df to return you wouldn't need this, it would automatically recycle
    botnum <- rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp] #added interval so if we ever want to change it, it's flexible
    return(data.frame(site =sitelist,
                      blank_time=sampleinfo$blanktime,
                      botnum = sprintf("%02d",botnum),
                      study=sampleinfo$study,samp_type=sampleinfo$samp_type)) #could include other label things in here too
  }
}




#If statement, defines the site_label function depending on if there is a known time or not
if (knowntime == TRUE){
  #Known time 
  site_label <- function(sampleinfo){
    sitelist <- rep(sampleinfo$site, times = sampleinfo$nsamp) #if you turn into a df to return you wouldn't need this, it would automatically recycle
    timestart <- as.POSIXct(sampleinfo$start)
    timelist <- seq(timestart, by=paste0(sampleinfo$interval, " hours"), length = sampleinfo$nsamp)
    botnum <- rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp] #added interval so if we ever want to change it, it's flexible
    return(data.frame(site =sitelist, timelist=timelist, 
                      botnum = sprintf("%02d",botnum), 
                      study=sampleinfo$study,
                      samp_type=sampleinfo$samp_type,
                      datetime_PST=format(timelist,"%Y%m%d%_%H%M"))) #could include other label things in here too
  }
} else if (knowntime == FALSE){
  site_label <- function(sampleinfo){
    sitelist <- rep(sampleinfo$site, times = sampleinfo$nsamp) #if you turn into a df to return you wouldn't need this, it would automatically recycle
    botnum <- rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp] #added interval so if we ever want to change it, it's flexible
    return(data.frame(site =sitelist,
                      blank_time=sampleinfo$blanktime,
                      botnum = sprintf("%02d",botnum),
                      study=sampleinfo$study,samp_type=sampleinfo$samp_type)) #could include other label things in here too
  }
}


#add sample names and IDs
if (knowntime == TRUE){
  data_labels<-data %>%
    mutate(sample_name=paste0( study, "_", "R", sprintf("%04d", randomn)),
           field_sample_ID = paste0(study, "_",
                                    site, "_",
                                    samp_type,"_",
                                    datetime_PST))
}else if (knowntime == FALSE ){
  data_labels<-data %>%
    mutate(sample_name=paste0( study, "_", "R", sprintf("%04d", randomn)),
           field_sample_ID_blank = paste0(study, "_",
                                          site, "_",
                                          samp_type,"_",
                                          blanktime))}

if(knowntime == FALSE){
  safe.write_csv(data_labels,
                 path = blank_fileformat,
                 row.names = FALSE)
  stop("DO NOT CONTINUE UNTIL YOU HAVE YOUR TIMES. SAVE YOUR CODE")}










