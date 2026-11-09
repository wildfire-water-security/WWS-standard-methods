## uses a fieldform and photo metadata to rename field photos with site names
  #written by Katie Wampler

#step 0: load libraries and functions ------
  library(exifr)
  library(stringr)
  library(data.table)
  library(lubridate)
  library(dplyr)
  library(boxrdrive)
  library(tidyr)

  get_photo_site <- function(photo_time, field_notes, maxdiff = 60*15){
    #only consider sites visited on the same day as photo
      photo_date <- as.Date(photo_time, tz=tz(field_notes$datetime))

    #get diff in time and get nearest sites before and after
    time_diff <- field_notes %>% filter(as.Date(datetime, tz=tz(datetime)) == photo_date) %>%
      mutate(diff = datetime - photo_time,
             abs_diff = abs(diff)) %>%
      arrange(abs_diff) %>% filter(abs_diff <= maxdiff)

    if(nrow(time_diff)  == 0){return(NA)}
    site <- time_diff$site[which.min(abs(time_diff$abs_diff))]
    return(site)
  }

#step 1: load photo info and site visit info -------
  photo_dir <- file.path(box_drive(),"Wildfire_Water_Security/02_Nodes/01_Empirical/06_Projects/LSCS-streambed-soil-leaching/photos/soil-sampling-20261103/drop-folder")
  field_notes <- read.csv("../WWS-Node1-LSCS-streambed-soil-leaching/data/metadata/fieldnotes/LSCS-SS-fieldnotes-20261103.csv")

  #ensure field notes has a datetime column with site visit time and site column
    field_notes <- field_notes %>% mutate(datetime = as.POSIXct(paste(Date, Time, sep=" "), format="%Y%m%d %H%M", tz="Etc/GMT+8"),
                                          site = Site)
      #should be in standard time

  #load photo info
    photos <- list.files(photo_dir, pattern="jpg|jpeg|png", full.names = TRUE)
    photo_meta <- read_exif(photos, tags=c("FileName", "DateTimeOriginal"))
    photo_meta <- photo_meta %>% mutate(datetime = as.POSIXct(DateTimeOriginal, format= "%Y:%m:%d %H:%M:%S", tz="Etc/GMT+8"))
    #use "Etc/GMT+8" for photos taken in standard time, "Etc/GMT+7" for photos in daylight time

#step 2: link photo time to visit time -----
  photo_meta$site <- sapply(photo_meta$datetime, get_photo_site, field_notes)

  #check for errors (photos with sd of time greater than 10 minutes)
  site_check <- photo_meta %>% group_by(site) %>% drop_na() %>% summarise(time_sd = sd(datetime)) %>%
    filter(time_sd > 60*10)

#step 3: rename photos with sitename ------
  #only rename photos with a site
  rename_files <- photo_meta %>% drop_na(site)  %>% mutate(ext = tools::file_ext(FileName),
                                      date_char = strftime(datetime, format="%Y%m%d_%H%M%S"),
                                      new_name = paste0(site, "_", date_char, ".", ext))

  save_loc <- file.path(box_drive(),"Wildfire_Water_Security/02_Nodes/01_Empirical/06_Projects/LSCS-streambed-soil-leaching/photos/soil-sampling-20261103/SWI")
  file.copy(rename_files$SourceFile, file.path(save_loc, rename_files$new_name)) #copy over so we don't lose photos
