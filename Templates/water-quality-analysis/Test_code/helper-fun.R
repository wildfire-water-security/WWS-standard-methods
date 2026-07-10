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

  #add sample names and IDs
  data_label <- function(sampleinfo){
    sitelist <- rep(sampleinfo$site, times = sampleinfo$nsamp) #if you turn into a df to return you wouldn't need this, it would automatically recycle
    botnum <- rep(1:24, times= ceiling(sampleinfo$nsamp/24))[1:sampleinfo$nsamp] #added interval so if we ever want to change it, it's flexible
    return(data.frame(site =sitelist,
                      blank_time=sampleinfo$blanktime,
                      botnum = sprintf("%02d",botnum),
                      study=sampleinfo$study,samp_type=sampleinfo$samp_type,
                      randomn= sample(rn:sn,sn, replace=F)
    )%>%
      mutate(sample_name = paste0( study, "_", "R", sprintf("%04d",randomn)),
             field_sample_ID_blank= paste0(study, "_",
                                                      site, "_",
                                                      samp_type,"_",
                                                       blanktime)))
  }
  

  
#section 2
  
if (pt2==TRUE){
  datasortlab

    site_label_blank <- function(sampletimeinfo){
      timestart <- as.POSIXct(sampletimeinfo$start)
      timelist <- seq(timestart, by=paste0(sampletimeinfo$interval, " hours"), length = sampletimeinfo$nsamp)
      return(
        data.frame(time_PST=as.character(timelist),datetime_PST=format(timelist,"%Y%m%d%_%H%M"))%>%
        cross_join(read.csv(blank_fileformat))%>%#imports file using previous file naming
        mutate(field_sample_ID = paste0(study, "_",
                                        site, "_",
                                        samp_type,"_",
                                        datetime_PST)) #add sample Ids with time
      )

      }


  

#if statements inside of functions


    samplelab <- data.frame(DOC=DOC,DOM=DOM,CCAL=CCAL,Levoglucosan=Levoglucosan, Excess=Excess, sn=sn)   

    adding in processing information:
    add analytes that are TRUE
    analysis_opt <- c("CN", "AQ", "NU", "LV", "EX")
    analysis_choice <- c(DOC, DOM, CCAL, Levoglucosan, Excess)

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
    datasortlab <- data_labels[order(data_labels$randomn),]
    datasortlab #use this to print datasheets out


    #repeat sample rows for each analysis
    lab_data <- datasortlab %>% slice(rep(1:n(), each = analysisn)) %>%
      mutate(datacodes,
             lab_sample_ID =
               paste0(field_sample_ID, "_", datacodes$code ))
    print(lab_data)


    
labfun <- function(samplelab){
  analysis_opt <- c("CN", "AQ", "NU", "LV", "EX")
  analysis_choice <- c(DOC, DOM, CCAL, Levoglucosan, Excess)
  analysiscodes <- analysis_opt[analysis_choice]
  analysisn <- length(analysiscodes) #how many analytes there are
  samples <- analysisn*sn #how many total filtered samples there will be
  lab_data <- datasortlab %>% slice(rep(1:n(), each = analysisn)) #duplicate sample so they can be repeated for each analyte
  return(
    data.frame(code = rep(analysiscodes,times=sn)) %>% #add labcodes
           mutate(analysis = case_when(
             code == "CN" ~ "Shimadzu",
             code == "AQ" ~ "Aqualog",
             code == "NU" ~ "CCAL",
             code == "LV" ~ "Levoglucosan",
             code == "EX" ~ "Excess",
             TRUE ~ NA_character_
           )) #%>%
           #cross_join(lab_data)
  )
}    
    
data_labels <- lapply(1:nrow(sampletimeinfo),function(x){site_label_blank(sampletimeinfo[x,])}) %>% #makes the frame
  bind_rows()
data_labels


lab_data <- datasortlab %>% slice(rep(1:n(), each = analysisn)) %>%
  mutate(datacodes,
         lab_sample_ID =
           paste0(field_sample_ID, "_", datacodes$code ))
print(lab_data)

  
  
  }
  