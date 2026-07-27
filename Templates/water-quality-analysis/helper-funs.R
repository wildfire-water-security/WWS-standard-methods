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
#' @examples
#' df <- read.csv("Templates/water-quality-analysis/examples/Blank_IS.csv")
#' saveloc <- "Templates/water-quality-analysis/example-template-is.xlsx"
#' create_datasheets(df, saveloc, "isco", "datetime", open=FALSE, overwrite=TRUE)
#' df <- read_xlsx("Templates/water-quality-analysis/examples/Final_IS.xlsx", sheet="field_labels")
#' create_datasheets(df, saveloc, "isco", "fieldlabel", open=FALSE)
#' df <- read_xlsx("Templates/water-quality-analysis/examples/Final_IS.xlsx", sheet="lab_labels")
#' create_datasheets(df, saveloc, "isco", "lablabel", open=FALSE)
#' create_datasheets(df, saveloc, "isco", "datasheet", open=TRUE)

#' saveloc <- "Templates/water-quality-analysis/example-template-gb.xlsx"
#' df <- read.csv("Templates/water-quality-analysis/examples/Blank_GB.csv")
#' create_datasheets(df, saveloc, "grab", "datetime", open=FALSE, overwrite=FALSE)
#' #df <- read_xlsx("Templates/water-quality-analysis/examples/Final_GB.xlsx", sheet="field_labels")
#' #create_datasheets(df, saveloc, "grab", "fieldlabel", open=FALSE)
#' df <- read_xlsx("Templates/water-quality-analysis/examples/Final_GB.xlsx", sheet="lab_labels")
#' create_datasheets(df, saveloc, "grab", "lablabel", open=FALSE)
#' create_datasheets(df, saveloc, "grab", "datasheet", open=TRUE)

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
