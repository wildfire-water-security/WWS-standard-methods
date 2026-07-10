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
#' @param open Logical, should excel file be opened on running?
#' @md
#' @details
#' If `type` is
#' - "datetime" df should have five columns:
#'  - Sample ID
#'  - Sample Name
#'  - Site
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

create_datasheets <- function(df, saveloc, type, open=FALSE){
  stopifnot(type %in% c("datetime", "fieldlabel", "lablabel", "datasheet"))

  #locate template
    template <- list.files(getwd(), pattern = "^wqsample-label-datasheet-templates[.]xlsx", full.names = TRUE, recursive = TRUE)
    if(length(template) == 0){
      #download from github as a backup
      template <- "https://github.com/wildfire-water-security/WWS-standard-methods/raw/refs/heads/main/Templates/water-quality-analysis/wqsample-label-datasheet-templates.xlsx"
    }

  #read in template
    tab <- wb_load(template)

  #add datetime
    if(type == "datetime"){
      #check to see if it was empty, if yes, don't overwrite user added info
      if(all(is.na(site_df$datetext_PST)) & all(is.na(site_df$time_text_PST))){
        site_df <- site_df %>% select(-c(datetext_PST, time_text_PST))
      }

      tab <- tab %>% wb_add_data(sheet= "Sample-DateTime", x=site_df, dims="A2")
    }

  #add labels
    if(type == "lablabel"){
      tab <- tab %>% wb_add_data(sheet= "Labels-Lab", x= site_df, dims="A2")
    }

    if(type == "fieldlabel"){
      tab <- tab %>% wb_add_data(sheet= "Labels-Field", x= site_df, dims="A2")
    }

  #add datasheets
    if(type == "datasheet"){
      ids <- tab %>% wb_to_df(sheet="Sample-DateTime") %>% select('Sample Name', 'Sample ID')
      datasheets <- datasheets %>% wb_add_data(sheet= "DS-Bottle Numbers", dims="A3", x=ids)
      datasheets <- datasheets %>% wb_add_data(sheet= "DS-ISCO Bottle Weight", dims="A3", x=ids)
      datasheets <- datasheets %>% wb_add_data(sheet= "DS-Filter Weight", dims="A3", x=ids)
    }

  #save
    wb_save(tab, savloc)
    print(paste0("information populated and saved to:\n", save_loc))

  #open if requested
    if(open){
      wb_open(tab)
    }




  }
