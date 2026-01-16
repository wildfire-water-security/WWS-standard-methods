## code to generate EEMs and absorbance long term standard and MDLs
library(eemanalyzeR)
library(pbapply)
library(dplyr)

#pull blank files from Aqualog folder and put in raw data -----
  dir <- "T:/Research/Aqualog_Data/2_PNNL_DOM"
  dates <- list.files(dir, "[0-9]{4}_[0-9]{2}_[0-9]{2}")
  #days with problems, obv blank cont. weird stuff
  date_rm <- c("2023_12_15", "2022_11_08", "2023_09_15b", "2023_03_24b", "2024_05_08", "2023_09_17", "2023_12_15b")

  dates <- dates[!(dates %in% date_rm)]
  output_dir <- "QAQC-files/Lab/EEMs/long-term-standards/blanks"
  get_blank_info <- function(file, output_dir){
    #get info about samples and metadata
    blk_files <- list.files(file.path(dir, file),"blk.*\\.dat|BLK.*\\.dat|Blk.*\\.dat", recursive = TRUE)

    meta_name <- list.files(file.path(dir, file), "metadata.rds", recursive = TRUE, full.names = TRUE)
    if(length(meta_name) > 0){
      meta_file <- readRDS(meta_name)

      #get metadata just for blanks
      meta_file <- meta_file[sapply(meta_file$data_identifier, function(x){any(grepl(x, blk_files))}),]
      meta_file <- meta_file[meta_file$data_identifier != "BLK",] #remove any samples with JUST BLK I don't want to deal
      meta_file <- meta_file[!(duplicated(meta_file$data_identifier) | duplicated(meta_file$data_identifier, fromLast = TRUE)), ]


      if(nrow(meta_file) > 0){
        #make unique name based on folder name (if run morn/even blanks are replicated)
        meta_file$long_term_name <- paste(file, meta_file$data_identifier, sep="_")

        blk_names <- blk_files
        for(x in 1:nrow(meta_file)){
          match <- grepl(paste0(meta_file$data_identifier[x], "_[0-9]_[0-9]s"), blk_names)
          blk_names[match] <- gsub(meta_file$data_identifier[x],
                                   meta_file$long_term_name[x], blk_names[match])
        }

        blk_name <- basename(blk_names)

        #copy files to new location
        abs <- grepl("1_Absorbance", blk_files)
        eem <- grepl("2_Blanks|3_Samples", blk_files)

        file.copy(file.path(dir, file, blk_files[abs]), paste(output_dir, "abs", blk_name[abs], sep="/"), overwrite = TRUE)
        file.copy(file.path(dir, file, blk_files[eem]), paste(output_dir, "eem", blk_name[eem], sep="/"), overwrite = TRUE)

        #write metadata
        meta_file <- meta_file %>% mutate(data_identifier = long_term_name) %>%
          dplyr::select(analysis_date, data_identifier, replicate_no, integration_time_s, dilution, RSU_area_1s) %>%unique()

        write.csv(meta_file, file.path(output_dir, paste0("metadata/", file, "_metadata.csv")), row.names = FALSE, quote = FALSE)

      }}
  }

  pblapply(dates, get_blank_info, output_dir=output_dir)

  #combine metadata and save
  meta <- list.files(file.path(output_dir, "metadata"))
  metadata <- lapply(meta, function(x){read.csv(file.path(output_dir, "metadata", x))}) %>% dplyr::bind_rows() %>%
    mutate(sample_type="sblank")
  write.csv(metadata, file.path(output_dir, "merged-blk-metadata.csv"), row.names=FALSE)

  #make mdl file for eemanalyzeR
    create_mdl("QAQC-files/Lab/EEMs/long-term-standards/blanks",
               recursive = TRUE, type="eem", iblank="_blank")

    create_mdl("QAQC-files/Lab/EEMs/long-term-standards/blanks",
               recursive = TRUE, type="abs", iblank="_blank")

#do the same for the tea samples -----
  #pull files from Aqualog folder and put in raw data
  dir <- "T:/Research/Aqualog_Data/2_PNNL_DOM"
  dates <- list.files(dir, "[0-9]{4}_[0-9]{2}_[0-9]{2}")
  output_dir <- "QAQC-files/Lab/EEMs/long-term-standards/tea-standards"
  #days with problems,
  date_rm <- c("2022_11_12", "2023_12_15", "2024_05_08", "2022_11_14", "2023_09_17", "2022_08_09", "2023_12_15b", "2022_08_06")
  dates <- dates[!(dates %in% date_rm)]

  get_tea_info <- function(file, output_dir){
    #get info about samples and metadata
    tea_files <- list.files(file.path(dir, file),"postTea.*\\.dat|preTea.*\\.dat", recursive = TRUE)

    meta_name <- list.files(file.path(dir, file), "metadata.rds", recursive = TRUE, full.names = TRUE)
    if(length(meta_name) > 0){
      meta_file <- readRDS(meta_name)

      #get metadata just for blanks
      meta_file <- meta_file[sapply(meta_file$data_identifier, function(x){any(grepl(x, tea_files))}),]

      if(nrow(meta_file) > 0){
        #make unique name based on folder name (if run morn/even blanks are replicated)
        meta_file$long_term_name <- paste(file, meta_file$data_identifier, sep="_")

        tea_names <- tea_files
        for(x in 1:nrow(meta_file)){
          tea_names <- gsub(meta_file$data_identifier[x], meta_file$long_term_name[x], tea_names)
        }

        tea_name <- basename(tea_names)

        #copy files to new location
        abs <- grepl("1_Absorbance", tea_files)
        eem <- grepl("2_Blanks|3_Samples", tea_files)

        file.copy(file.path(dir, file, tea_files[abs]), paste(output_dir, "abs", tea_name[abs], sep="/"), overwrite = TRUE)
        file.copy(file.path(dir, file, tea_files[eem]), paste(output_dir, "eem", tea_name[eem], sep="/"), overwrite = TRUE)

        #write metadata
        meta_file <- meta_file %>% mutate(data_identifier = long_term_name) %>%
          dplyr::select(analysis_date, data_identifier, replicate_no, integration_time_s, dilution, RSU_area_1s) %>%unique()

        write.csv(meta_file, file.path(output_dir, paste0("metadata/", file, "_metadata.csv")), row.names = FALSE, quote = FALSE)

      }}
  }

  pblapply(dates, get_tea_info, output_dir=output_dir)

  #combine metadata and save
  meta <- list.files(file.path(output_dir, "metadata"))
  metadata <- lapply(meta, function(x){read.csv(file.path(output_dir, "metadata", x))}) %>% dplyr::bind_rows() %>%
    mutate(sample_type="check")

  write.csv(metadata, file.path(output_dir, "merged-tea-metadata.csv"), row.names=FALSE)

  #make std file for eemanalyzeR
  create_std("QAQC-files/Lab/EEMs/long-term-standards/tea-standards", recursive = TRUE,
             type="eem", iblank="_blank", abs_pattern="Abs")

  create_std("QAQC-files/Lab/EEMs/long-term-standards/tea-standards", recursive = TRUE,
             type="abs", abs_pattern="Abs")
