
#' Generate Standardized Sample IDs
#'
#' Uses provided information to create standardized sample ID's based on the the WWS Node 1
#' sample naming procedure.
#'
#' @param study 4-5 letter code associated with the study
#' @param sample_type type of sample being collected options include `GB`, `IS`, `SS`
#' @param sites either a file path to a `.csv` file or a vector of site names
#' @param start a vector of date/times or blank to use for sample collection
#' @param nLB number of lab blanks
#' @param nFB number of field blanks
#' @param rn number to start the random number sequence. This is used to assign a random analysis order.
#' @param nsamp number of samples that will be collected at each site, only used if `sample_type` is `IS`.
#' @param interval number of hours between each sample taken, only used if `sample_type` is `IS`.
#'
#' @returns a `data.frame` containing information required to generate the sampleIDs
#' @export
#'
#' @examples
generate_ids <- function(study, sample_type, sites, start, nLB, nFB, rn, nsamp=NULL, interval=NULL){
  stopifnot(is.character(study), sample_type %in% c("GB", "IS", "SS"),
            is.numeric(nLB), is.numeric(nFB), is.numeric(rn))

  #if sites are a .csv, extract info
    if(length(sites) == 1 && grepl(".csv", sites)){
      path <- list.files(here(), pattern=sites, recursive = TRUE)
      sites <- read.csv(file.path(here(), path))[,1]
    }

  #if isco and we know times, fill in
    if(samp_type == "IS"){
      sitedf <- data.frame(site=sites, nsamp = nsamp, start=start, interval=interval)
      start <- as.vector(apply(sitedf, MARGIN = 1, FUN = get_isco_times))
      if(any(grepl("^\\d{10}$", start))){start <- as.POSIXct(start)}

      #generate df of sample IDs
      ids <- data.frame(study=study,
                        site= rep(sites, times=nsamp),
                        samp_type = samp_type,
                        samp_time = start)
    }else{
      #generate df of sample IDs
      ids <- data.frame(study=study,
                        site= sites,
                        samp_type = samp_type,
                        samp_time = start)
    }

  #generate random IDs
    randomid <- seq(from=rn, by=1, length=nrow(ids))
    ids$randomn <- sample(randomid,nrow(ids), replace=F)

  #replace date and time if available
    knowntime <- is.POSIXct(ids$samp_time)
    if(knowntime){
      ids <- ids %>% mutate(date = strftime(samp_time, format="%Y%m%d"),
                            time = strftime(samp_time, format="%H%M"),
                            samp_time = paste(date,time, sep="_"))
    }else{
      ids <- ids %>% mutate(date = "",
                            time = "")
    }

  #create full ids
    ids <- ids %>% mutate(
      sample_name=paste0(study, "_", "R", sprintf("%04d", randomn)),
      field_sample_ID = paste0(study, "_",site, "_",samp_type,"_",samp_time))

  return(ids)

}

#get isco times for each site
#' Get ISCO times for a site
#'
#' @param df a data.frame with info required to generate the sequence of times.
#'
#' @returns a vector of times

get_isco_times <- function(df){
  df <- data.frame(as.list(df))

  #if time specified fill in
  if(grepl("^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}$", df$start)){
    df$start <- as.POSIXct(df$start)
    times <- seq(df$start, by=paste0(as.numeric(df$interval), " hours"), length = as.numeric(df$nsamp))
  }else{
    #if time not specified, just repeat
    times <- rep(df$start, times=df$nsamp)
  }

  return(times)
}

#packages needed to run code
library(dplyr)
library(stringr)
library(lubridate)
library(writexl)
library(openxlsx2)

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

#' Create sample IDs and names from sampling information
#'
#' @param sites a `data.frame`
#'
#' @returns

generateids <- function(sites){}

#' Populate excel datasheets for printing
#' @param df `data.frame` containing the info to populate the table with. See details for what df should contain for each type.
#' @param saveloc File path where file should be saved, including file name.
#' @param type Character specifying the info to enter into template. Either
#' "datetime", "lablabel", "fieldlabel", or "datasheet".
#' @param samptype Character specifying the sample type. Either "isco" or "grab".
#' @param open Logical, should excel file be opened on running?
#' @param overwrite Logical, should existing template be overwritten?
#' @md
#' @details
#' If `type` is
#' - "datetime" df should have five columns:
#'  - Sample ID
#'  - Sample Name
#'  - Site
#'  - Random ID
#'  - Bottle Number (will be blank for grab)
#'  - Date (may be blank)
#'  - Time (may be blank)
#' - "lablabel" df should have three columns:
#'  - Sample ID
#'  - Sample Name
#'  - Bottle Number
#' - "fieldlabel" df should have three columns:
#'  - Sample ID
#'  - Sample Name
#'  - Analysis
#' - "datasheet" no df is needed, but dates and times should be filled out
#' @returns Writes an `.xlsx` file to save location and opens up the file in excel.


create_datasheets <- function(df, saveloc, samptype, type, open=FALSE, overwrite=FALSE){
  stopifnot(type %in% c("datetime", "fieldlabel", "lablabel", "datasheet"),
            samptype %in% c("isco", "grab"))

  #locate template/workbook
    if(file.exists(saveloc) & !overwrite){
      #load existing workbook
      path <- saveloc

    }else{
      #get template
      path <- list.files(getwd(), pattern = "^wqsample-label-datasheet-templates[.]xlsx", full.names = TRUE, recursive = TRUE)
      if(length(path) == 0){
        #download from github as a backup
        path <- "https://github.com/wildfire-water-security/WWS-standard-methods/raw/refs/heads/main/Templates/water-quality-analysis/wqsample-label-datasheet-templates.xlsx"
      }
    }

  #read in template
    tab <- wb_load(path)

  #add datetime
    if(type == "datetime"){
      if(samptype == "isco"){
        df <- df %>% mutate(date_text_PST="", time_text_PST="",.after=botnum)}
      if(samptype == "grab"){
        df <- df %>% mutate(botnum="N/A", .before=randomn)}
      site_df <- df %>% select(c("sample_name", "field_sample_ID_blank",
                                 "site", "botnum","randomn", "date_text_PST", "time_text_PST"))
      #don't overwrite dates/times
      exist_data <- tab %>% wb_to_df(sheet="Sample-DateTime")
      if(any(!is.na(exist_data$`Date (YYYYMMDD)`)) | any(!is.na(exist_data$`Time (HHMM)`))){
           site_df <- site_df %>% select(-c(date_text_PST, time_text_PST))
         }
      tab <- tab %>% wb_add_data(sheet= "Sample-DateTime",
                                 x=site_df, dims="A2",col_names = FALSE)
    }

  #add labels
    if(type == "lablabel"){
      tab <- tab %>% wb_add_data(sheet= "Labels-Lab", x= df, dims="A2", col_names=FALSE)
    }

    if(type == "fieldlabel"){
      tab <- tab %>% wb_add_data(sheet= "Labels-Field", x= df, dims="A2", col_names=FALSE)
    }

  #add datasheets
    if(type == "datasheet"){
      ids <- tab %>% wb_to_df(sheet="Sample-DateTime") %>% select('Sample ID','Sample Name') %>% arrange(`Sample ID`)

      rows <- nrow(ids) #number of formatted lines we need #formatting in table goes to 500

      #prepare bottle numbers df
      n_page <- ceiling(rows / 37)
      clear <- (n_page*38 + 1)
      tab <- tab %>% wb_add_data(sheet= "DS-Bottle Numbers", dims="A3", x=ids, col_names=FALSE) %>%  #36/38 per page
        wb_clean_sheet(dims = paste0("A", clear, ":G500"))

      #prepare bottle weight df
      n_page <- ceiling(rows / 52)
      clear <- (n_page*53 + 1)
      tab <- tab %>% wb_add_data(sheet= "DS-ISCO Bottle Weight", dims="A3", x=ids, col_names=FALSE) %>% #51/53 per page
        wb_clean_sheet(dims = paste0("A", clear, ":G500"))

      #prepare filter weight df
      n_page <- ceiling(rows / 53)
      clear <- (n_page*54 + 1)
      tab <- tab %>% wb_add_data(sheet= "DS-Filter Weight", dims="A3", x=ids, col_names=FALSE) %>% #52/54 per page
        wb_clean_sheet(dims = paste0("A", clear, ":G500"))

    }

  #save
    wb_save(tab, file=saveloc)


  #open if requested
    if(open){
      shell.exec(file.path(getwd(), saveloc))
    }else{
      cat(paste0("information populated and saved to:\n", saveloc))
    }




  }
