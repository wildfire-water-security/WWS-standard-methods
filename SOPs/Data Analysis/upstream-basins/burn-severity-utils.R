#' Clean MTBS wildfire files
#'
#' MTBS downloads come in multiple zipped folders and the file names are long and meaningless. This
#' will pull all files in a single folder and rename files to be more user friendly.
#'
#'@importFrom sf st_read
#'@importFrom stringr str_detect str_replace_all str_split regex

#' @param wd the location with the raw MTBS files
#' @param save_loc the location to save all the files
#' @export
#' @examples
#'

clean_mtbs <- function(wd, save_loc){
  stopifnot(file.exists(wd), file.exists(save_loc))

  #get zipped files of fires
  files <- list.files(wd, recursive = T)
  files <- files[stringr::str_detect(files, ".zip")]

  pbapply::pblapply(1:length(files), function(f){
    x <- files[f]

    #get place to put zip files
    exdir <- file.path(wd, tools::file_path_sans_ext(basename(x)))

    #upzip file
    unzip(file.path(wd, x), exdir=exdir)

    #find shape file to extract fire name
    zip_files <- list.files(exdir)
    shp_file <- zip_files[stringr::str_detect(zip_files, "burn_bndy.shp")]

    #read file to get fire name
    name <- sf::read_sf(paste(exdir, shp_file, sep="/"))[1,]
    filename <- name$Incid_Name
    state <- substr(name$Event_ID, 1,2)
    year <- lubridate::year(name$Ig_Date)

    if(filename == "UNNAMED"){
      filename <- str_sub(name$Event_ID, 1, 13)
    }

    #format name for file
    filename <- tolower(filename)
    filename <- stringr::str_replace_all(filename, " ", "_")
    filename <- gsub("[][()]|[.]", "", filename)
    filename <- gsub("[/]", "-", filename)

    #rename files
    #get dates of data
    start <- regmatches(name$Pre_ID, regexpr( "(19|20)\\d{6}(?=(_\\d+m)?$)", name$Pre_ID, perl = TRUE))
    end <- regmatches(name$Post_ID, regexpr( "(19|20)\\d{6}(?=(_\\d+m)?$)", name$Post_ID, perl = TRUE))

    dates <- c(start, end)

    #get header, remove files from other fires
    header <- stringr::str_split_1(shp_file, "_")[1]
    zip_files <- zip_files[stringr::str_detect(zip_files, stringr::regex(header, ignore_case=T))]

    #replace long name with fire name
    new_files <- gsub(header, filename, zip_files, ignore.case = TRUE)

    #if both, remove completely
    if(all(!is.na(dates)) & length(dates) == 2){
      new_files <- gsub(paste("_", paste(dates, collapse ="_"), sep=""), "", new_files)
    }

    #if just one use pre or post
    if(length(start) == 1 && !is.na(start)){new_files <- gsub(start, "pre", new_files)}
    if(length(end) == 1 && !is.na(end)){new_files <- gsub(end, "post", new_files)}

    #replace nbr with dnbr
    new_files <- gsub("_nbr", "_dnbr", new_files)

    #rename files
    fire_file <- file.path(save_loc, paste0(state,"_", year, "_", filename))
    dir.create(fire_file, showWarnings = FALSE)
    file.rename(file.path(exdir, zip_files), file.path(fire_file, new_files))

  })
}


#' Interpolate missing MTBS dNBR data
#'
#' Sometimes a raster layer will have areas of "Non-processing" due to data gaps in the Landsat data used to determine the dNBR. While we could just remove
#' these layers, sometimes you don't want to lose that data. This procedure aims to fill those gaps using interpolation using MBA::mba.points() which
#' uses a bivariate scatter of points using multilevel B-splines to interpolate.
#'
#' @param raster Object of class `SpatRaster` to interpolate.
#' @param plots Logical, should plots showing the original raster and the interpolated raster be printed?
#'
#' @returns Object of class `SpatRaster` with the data gaps interpolated.
#' @export
#' @md
#'
#' @examples
fill_mtbs <- function(raster, plots = TRUE){
  stopifnot(inherits(raster, "SpatRaster"))
  library(raster)
  library(pracma)
  library(MBA)
  library(dplyr)

  #interpolate NA values
  # Convert raster to a data frame of coordinates and values
  xyz <- as.data.frame(raster, xy = TRUE)

  #remove NA values (MTBS sets to -32768)
  xy.est <- xyz[xyz[,3] == -32768, 1:2]
  xyz <- xyz[xyz[,3] != -32768,]
  colnames(xyz) <- c("x", "y", "z")

  # Perform interpolation using mba.points
  interp_result <- MBA::mba.points(
    xyz, #existing data
    xy.est, #missing data
  )

  # Convert interpolated results back to a raster
  raster_filled <- rbind(xyz, interp_result$xyz.est)
  raster_filled <- terra::rast(raster_filled, crs=crs(raster))


  if(plots){
    #show comparison
    raster_plot <- raster
    NAflag(raster_plot) <- -32768
    p1 <- ggplot2::ggplot() + tidyterra::geom_spatraster(data=raster_plot) + ggplot2::scale_fill_viridis_c()
    p2 <- ggplot2::ggplot() + tidyterra::geom_spatraster(data=raster_filled) + ggplot2::scale_fill_viridis_c()

    p <- ggpubr::ggarrange(p1, p2, nrow=1, common.legend = TRUE)
    print(p)
  }

  return(raster_filled)
}
